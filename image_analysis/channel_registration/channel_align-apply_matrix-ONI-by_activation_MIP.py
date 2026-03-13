from tifffile import TiffFile, imwrite
from tkinter import filedialog as fd
from rich.progress import Progress
from concurrent.futures import ThreadPoolExecutor, as_completed
import cv2
import pickle
import numpy as np
import psutil
import gc
import time


"""
Activation-Scheme Channel-Aligned TIFF Splitter — MIP variant
==============================================================
Like split_channel_by_activation_MIP.py, but each camera half is also passed
through a perspective-warp alignment matrix (loaded from a .pickle file) and
then border-cropped by 5px on each side — exactly as in
channel_align-apply_matrix-ONI.py.

The final output per activation+camera is a Maximum Intensity Projection
(single 2D frame) of the aligned, cropped stack.

User-editable config (top of this file)
---------------------------------------
  TRANSFORM_CAMERA   : "left" or "right" — which camera side the warp matrix
                       is applied to. The other side is treated as reference
                       (only border-cropped, no warp). This applies uniformly
                       across all activations.

Each activation entry in ACTIVATION_SCHEME specifies:
  - frame_start  : int   -- 1-based, inclusive
  - frame_end    : int   -- 1-based, inclusive
  - cameras      : list  -- subset of ["left", "right"]
  - labels       : list  -- one label per camera (output filename suffix)

Output filenames:  <input_stem>-<label>-MIP.tif
"""


# ─────────────────────────────────────────────────────────────────────────────
# USER-DEFINED CONFIG  ← edit this section before running
# ─────────────────────────────────────────────────────────────────────────────

# Which camera the warp matrix is applied to: "left" or "right"
TRANSFORM_CAMERA = "left"

ACTIVATION_SCHEME = [
    {
        "frame_start": 1,
        "frame_end": 10,
        "cameras": ["right"],
        "labels": ["Nb_JFX650"],
    },
    {
        "frame_start": 11,
        "frame_end": 20,
        "cameras": ["left"],
        "labels": ["Edar-GFP"],
    },
    {
        "frame_start": 21,
        "frame_end": 30,
        "cameras": ["left"],
        "labels": ["DAPI"],
    },
]

# ─────────────────────────────────────────────────────────────────────────────


# ── Alignment helpers (from channel_align-apply_matrix-ONI.py) ───────────────


def _crop_image(img):
    """Crop 5-pixel border from a 2D (H, W) image."""
    h, w = img.shape
    return img[5 : h - 5, 5 : w - 5]


def _warp_frame(frame, warp_matrix):
    """Apply perspective warp to a single 2D frame."""
    h, w = frame.shape
    return cv2.warpPerspective(
        frame,
        warp_matrix,
        (w, h),
        flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP,
    )


# ── Scheme validation ─────────────────────────────────────────────────────────


def _validate_scheme(scheme, transform_camera):
    if transform_camera not in ("left", "right"):
        raise ValueError(
            f"TRANSFORM_CAMERA must be 'left' or 'right', got '{transform_camera}'."
        )
    for i, act in enumerate(scheme):
        required = {"frame_start", "frame_end", "cameras", "labels"}
        missing = required - set(act.keys())
        if missing:
            raise ValueError(f"Activation {i+1} is missing keys: {missing}")

        if act["frame_start"] > act["frame_end"]:
            raise ValueError(
                f"Activation {i+1}: frame_start ({act['frame_start']}) "
                f"> frame_end ({act['frame_end']})"
            )

        for cam in act["cameras"]:
            if cam not in ("left", "right"):
                raise ValueError(
                    f"Activation {i+1}: unknown camera '{cam}'. Must be 'left' or 'right'."
                )

        if len(act["cameras"]) != len(act["labels"]):
            raise ValueError(f"Activation {i+1}: len(cameras) != len(labels).")


# ── Memory helpers ────────────────────────────────────────────────────────────


def _calculate_chunk_size(max_concurrent_chunks=4):
    available_ram_gb = psutil.virtual_memory().available / (1024**3)
    usable_ram_gb = available_ram_gb * 0.8
    ram_per_chunk = usable_ram_gb / (max_concurrent_chunks + 1)
    frame_size_mb = 1.4
    frames_per_chunk = int((ram_per_chunk * 1024) / frame_size_mb)
    return max(500, min(frames_per_chunk, 3000))


# ── Main processor ────────────────────────────────────────────────────────────


class ActivationAlignedMIPProcessor:
    """
    For each activation:
      1. Load the specified frame range in parallel chunks.
      2. Split each frame into left / right halves.
      3. Apply warp alignment to the TRANSFORM_CAMERA side.
      4. Crop 5px border from all selected cameras.
      5. Compute MIP (np.max across frames) and write one 2D TIFF per camera.
    """

    def __init__(
        self, fpath, scheme, transform_camera, warp_matrix, max_concurrent_chunks=4
    ):
        self.fpath = fpath
        self.scheme = scheme
        self.transform_camera = transform_camera
        self.warp_matrix = warp_matrix
        self.max_concurrent_chunks = max_concurrent_chunks
        self.chunk_size = _calculate_chunk_size(max_concurrent_chunks)

    # ── workers ──────────────────────────────────────────────────────────────

    def _load_frames(self, tif, frame_indices):
        return [tif.pages[idx].asarray() for idx in frame_indices]

    def _split_align_crop_frames(self, frames, halfwidth, cameras):
        """
        For each frame:
          - split into left ([:, :halfwidth]) and right ([:, halfwidth:])
          - warp the side matching self.transform_camera
          - crop 5px border from all selected cameras
        Returns {cam: [2D frames]}
        """
        result = {cam: [] for cam in cameras}
        for frame in frames:
            halves = {
                "left": frame[:, :halfwidth],
                "right": frame[:, halfwidth:],
            }
            for cam in cameras:
                half = halves[cam]
                if cam == self.transform_camera:
                    half = _warp_frame(half, self.warp_matrix)
                half = _crop_image(half)
                result[cam].append(half)
        return result

    # ── per-activation ────────────────────────────────────────────────────────

    def _process_activation(
        self, tif, activation, halfwidth, act_index, total_acts, progress, task
    ):
        frame_start = activation["frame_start"] - 1  # 0-based
        frame_end = activation["frame_end"]  # exclusive
        cameras = activation["cameras"]
        labels = activation["labels"]

        num_frames = frame_end - frame_start
        cam_to_label = dict(zip(cameras, labels))

        stem = self.fpath[:-4] if self.fpath.lower().endswith(".tif") else self.fpath
        out_paths = {cam: f"{stem}-{cam_to_label[cam]}-MIP.tif" for cam in cameras}

        print(
            f"  Activation {act_index}/{total_acts}: "
            f"frames {activation['frame_start']}-{activation['frame_end']} | "
            f"cameras: {cameras} | transform: {self.transform_camera}"
        )
        for cam, path in out_paths.items():
            print(f"    → {cam}: {path.split('/')[-1]}")

        # ── chunked parallel I/O ──────────────────────────────────────────────
        all_chunk_results = {}

        frame_indices_all = list(range(frame_start, frame_end))
        chunks_info = []
        for i in range(0, num_frames, self.chunk_size):
            chunk_indices = frame_indices_all[i : i + self.chunk_size]
            chunks_info.append((len(chunks_info), chunk_indices))

        with ThreadPoolExecutor(
            max_workers=self.max_concurrent_chunks,
            thread_name_prefix="IO-Loader",
        ) as io_executor, ThreadPoolExecutor(
            max_workers=min(4, psutil.cpu_count() or 1),
            thread_name_prefix="Processor",
        ) as proc_executor:

            io_futures = {
                io_executor.submit(self._load_frames, tif, ci[1]): ci[0]
                for ci in chunks_info
            }

            proc_futures = {}
            for completed_io in as_completed(io_futures):
                chunk_id = io_futures[completed_io]
                try:
                    frames = completed_io.result()
                    pf = proc_executor.submit(
                        self._split_align_crop_frames,
                        frames,
                        halfwidth,
                        cameras,
                    )
                    proc_futures[pf] = chunk_id
                except Exception as e:
                    print(f"    ✗ Error loading chunk {chunk_id}: {e}")

            for completed_proc in as_completed(proc_futures):
                chunk_id = proc_futures[completed_proc]
                try:
                    all_chunk_results[chunk_id] = completed_proc.result()
                except Exception as e:
                    print(f"    ✗ Error processing chunk {chunk_id}: {e}")

        # ── assemble → MIP → write ────────────────────────────────────────────
        for cam in cameras:
            ordered_frames = []
            for cid in sorted(all_chunk_results.keys()):
                ordered_frames.extend(all_chunk_results[cid][cam])

            arr = np.stack(ordered_frames, axis=0)  # (frames, H, W)
            mip = np.max(arr, axis=0)  # (H, W)
            del arr
            gc.collect()
            imwrite(out_paths[cam], mip, imagej=True)
            del mip
            gc.collect()
            progress.advance(task)

        del all_chunk_results
        gc.collect()

    # ── public entry ──────────────────────────────────────────────────────────

    def process(self, progress, task):
        with TiffFile(self.fpath) as tif:
            num_frames = len(tif.pages)
            first_frame = tif.pages[0].asarray()
            height, width = first_frame.shape
            halfwidth = width // 2

            print(f"  File: {self.fpath.split('/')[-1]}")
            print(
                f"  Frames: {num_frames} | Resolution: {width}×{height} | Halfwidth: {halfwidth}"
            )

            for act in self.scheme:
                if act["frame_end"] > num_frames:
                    raise ValueError(
                        f"Activation frame_end={act['frame_end']} exceeds "
                        f"total frames in file ({num_frames})."
                    )

            total_acts = len(self.scheme)
            for idx, activation in enumerate(self.scheme, start=1):
                self._process_activation(
                    tif, activation, halfwidth, idx, total_acts, progress, task
                )


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────


def process_files():
    _validate_scheme(ACTIVATION_SCHEME, TRANSFORM_CAMERA)

    print("Choose the warp matrix (.pickle) file:")
    matrix_path = fd.askopenfilename(
        title="Choose warp matrix",
        filetypes=[("Pickle files", "*.p *.pickle *.pkl"), ("All files", "*.*")],
    )
    if not matrix_path:
        print("No matrix selected. Exiting.")
        return

    with open(matrix_path, "rb") as f:
        warp_matrix = pickle.load(f)
    print(f"  Loaded matrix: {matrix_path.split('/')[-1]}")

    print("\nChoose the TIF file(s) to process:")
    lst_files = list(fd.askopenfilenames(filetypes=[("TIFF files", "*.tif *.tiff")]))
    if not lst_files:
        print("No files selected. Exiting.")
        return

    total_outputs = sum(len(act["cameras"]) for act in ACTIVATION_SCHEME) * len(
        lst_files
    )

    start_time = time.time()
    with Progress() as progress:
        task = progress.add_task("[green]Writing output files", total=total_outputs)

        for i, fpath in enumerate(lst_files, 1):
            fname = fpath.split("/")[-1]
            print(f"\n[{i}/{len(lst_files)}] {fname}")
            try:
                processor = ActivationAlignedMIPProcessor(
                    fpath,
                    ACTIVATION_SCHEME,
                    TRANSFORM_CAMERA,
                    warp_matrix,
                    max_concurrent_chunks=4,
                )
                processor.process(progress, task)

            except Exception as e:
                print(f"❌ Error processing {fpath}: {e}")

    elapsed = time.time() - start_time
    print(f"\n✅ All done in {elapsed:.1f}s")


if __name__ == "__main__":
    process_files()
