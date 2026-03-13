from tifffile import TiffFile, imwrite
from tkinter import filedialog as fd
from rich.progress import Progress
from concurrent.futures import ThreadPoolExecutor, as_completed
from queue import Queue
import numpy as np
import psutil
import gc
import time


"""
Activation-Scheme TIFF Splitter — MIP variant
==============================================
Identical to split_channel_by_activation.py, except each output file is a
Maximum Intensity Projection (single 2D frame) rather than the full stack.

Each activation entry in the scheme specifies:
  - frame_start : int   -- 1-based, inclusive start frame of this activation
  - frame_end   : int   -- 1-based, inclusive end frame of this activation
  - cameras     : list  -- e.g. ["left"], ["right"], or ["left", "right"]
  - labels      : list  -- one label per camera entry, used as the output filename suffix
                           e.g. ["targetA_dyeA"] or ["targetA_dyeA", "targetB_dyeB"]

Example activation scheme
--------------------------
ACTIVATION_SCHEME = [
    {
        "frame_start": 1,
        "frame_end":   10,
        "cameras":     ["left", "right"],
        "labels":      ["targetA_dyeA", "targetB_dyeB"],
    },
    {
        "frame_start": 11,
        "frame_end":   20,
        "cameras":     ["left"],
        "labels":      ["targetC_dyeC"],
    },
]

Output filenames are:   <input_stem>-<label>-MIP.tif
"""


# ─────────────────────────────────────────────────────────────────────────────
# USER-DEFINED ACTIVATION SCHEME  ← edit this section before running
# ─────────────────────────────────────────────────────────────────────────────

ACTIVATION_SCHEME = [
    {
        "frame_start": 1,
        "frame_end":   10,
        "cameras":     ["right"],
        "labels":      ["Nb_JFX650"],
    },
    {
        "frame_start": 11,
        "frame_end":   20,
        "cameras":     ["left"],
        "labels":      ["Edar-GFP"],
    },
    # {
    #     "frame_start": 21,
    #     "frame_end":   30,
    #     "cameras":     ["left"],
    #     "labels":      ["DAPI"],
    # },
]

# ─────────────────────────────────────────────────────────────────────────────


def _validate_scheme(scheme):
    """Validate the activation scheme for obvious errors."""
    for i, act in enumerate(scheme):
        required_keys = {"frame_start", "frame_end", "cameras", "labels"}
        missing = required_keys - set(act.keys())
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
                    f"Activation {i+1}: unknown camera '{cam}'. "
                    "Must be 'left' or 'right'."
                )

        if len(act["cameras"]) != len(act["labels"]):
            raise ValueError(
                f"Activation {i+1}: len(cameras) != len(labels). "
                f"Got cameras={act['cameras']}, labels={act['labels']}"
            )


def _calculate_chunk_size(max_concurrent_chunks=4):
    """Calculate chunk size based on available system RAM."""
    available_ram_gb = psutil.virtual_memory().available / (1024**3)
    usable_ram_gb = available_ram_gb * 0.8
    ram_per_chunk = usable_ram_gb / (max_concurrent_chunks + 1)
    frame_size_mb = 1.4
    frames_per_chunk = int((ram_per_chunk * 1024) / frame_size_mb)
    return max(500, min(frames_per_chunk, 3000))


class ActivationProcessor:
    """
    Processes a single TIFF file according to an activation scheme.

    For each activation entry the processor:
      1. Reads the specified frame range.
      2. Splits each frame vertically into left / right halves.
      3. Keeps only the camera(s) specified for that activation.
      4. Writes one output TIFF per kept camera, named with the corresponding label.
    """

    def __init__(self, fpath, scheme, max_concurrent_chunks=4):
        self.fpath = fpath
        self.scheme = scheme
        self.max_concurrent_chunks = max_concurrent_chunks
        self.chunk_size = _calculate_chunk_size(max_concurrent_chunks)

    # ------------------------------------------------------------------
    # Low-level workers
    # ------------------------------------------------------------------

    def _load_frames(self, tif, frame_indices):
        """Load a list of frame indices from an open TiffFile."""
        return [tif.pages[idx].asarray() for idx in frame_indices]

    def _split_frames(self, frames, halfwidth, cameras):
        """
        Split frames vertically and return only the requested camera halves.

        Returns a dict: {camera_side: [frames]}
        """
        result = {cam: [] for cam in cameras}
        for frame in frames:
            if "left" in cameras:
                result["left"].append(frame[:, :halfwidth])
            if "right" in cameras:
                result["right"].append(frame[:, halfwidth:])
        return result

    # ------------------------------------------------------------------
    # Per-activation processing
    # ------------------------------------------------------------------

    def _process_activation(self, tif, activation, halfwidth, act_index, total_acts, progress, task):
        """
        Load, split, and save frames for a single activation entry.
        Uses chunked parallel I/O to stay memory-efficient.
        Advances `task` in the shared `progress` bar once per output file written.
        """
        frame_start = activation["frame_start"] - 1   # convert to 0-based
        frame_end   = activation["frame_end"]          # exclusive upper bound (0-based)
        cameras     = activation["cameras"]
        labels      = activation["labels"]

        num_frames = frame_end - frame_start
        cam_to_label = dict(zip(cameras, labels))

        # Build output paths (-MIP suffix)
        stem = self.fpath[:-4] if self.fpath.lower().endswith(".tif") else self.fpath
        out_paths = {
            cam: f"{stem}-{cam_to_label[cam]}-MIP.tif"
            for cam in cameras
        }

        print(
            f"  Activation {act_index}/{total_acts}: "
            f"frames {activation['frame_start']}-{activation['frame_end']} | "
            f"cameras: {cameras}"
        )
        for cam, path in out_paths.items():
            print(f"    → {cam}: {path.split('/')[-1]}")

        # ---- chunked parallel I/O ----------------------------------------
        all_chunk_results = {}   # chunk_id -> {cam: [frames]}

        frame_indices_all = list(range(frame_start, frame_end))
        chunks_info = []
        for i in range(0, num_frames, self.chunk_size):
            chunk_indices = frame_indices_all[i: i + self.chunk_size]
            chunks_info.append((len(chunks_info), chunk_indices))

        with ThreadPoolExecutor(
            max_workers=self.max_concurrent_chunks,
            thread_name_prefix="IO-Loader",
        ) as io_executor, ThreadPoolExecutor(
            max_workers=min(4, psutil.cpu_count() or 1),
            thread_name_prefix="Processor",
        ) as proc_executor:

            # Submit I/O jobs
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
                        self._split_frames, frames, halfwidth, cameras
                    )
                    proc_futures[pf] = chunk_id
                except Exception as e:
                    print(f"    ✗ Error loading chunk {chunk_id}: {e}")

            for completed_proc in as_completed(proc_futures):
                chunk_id = proc_futures[completed_proc]
                try:
                    split_result = completed_proc.result()
                    all_chunk_results[chunk_id] = split_result
                except Exception as e:
                    print(f"    ✗ Error processing chunk {chunk_id}: {e}")

        # ---- assemble, MIP, & write (advances shared progress bar per output file) ----
        for cam in cameras:
            ordered_frames = []
            for cid in sorted(all_chunk_results.keys()):
                ordered_frames.extend(all_chunk_results[cid][cam])

            arr = np.stack(ordered_frames, axis=0)   # (frames, H, W)
            mip = np.max(arr, axis=0)                 # (H, W)
            del arr
            gc.collect()
            imwrite(out_paths[cam], mip, imagej=True)
            del mip
            gc.collect()
            progress.advance(task)   # one tick per output file written

        del all_chunk_results
        gc.collect()

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def process(self, progress, task):
        """Process the file, reporting progress into the shared progress bar."""
        with TiffFile(self.fpath) as tif:
            num_frames = len(tif.pages)
            first_frame = tif.pages[0].asarray()
            height, width = first_frame.shape
            halfwidth = width // 2

            print(f"  File: {self.fpath.split('/')[-1]}")
            print(f"  Frames: {num_frames} | Resolution: {width}×{height} | Halfwidth: {halfwidth}")

            # Validate frame indices against actual file length
            for act in self.scheme:
                if act["frame_end"] > num_frames:
                    raise ValueError(
                        f"Activation frame_end={act['frame_end']} exceeds "
                        f"total frames in file ({num_frames})."
                    )

            total_acts = len(self.scheme)
            for idx, activation in enumerate(self.scheme, start=1):
                self._process_activation(tif, activation, halfwidth, idx, total_acts, progress, task)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def process_files():
    # Validate scheme once before touching any files
    _validate_scheme(ACTIVATION_SCHEME)

    print("Choose the TIF file(s) to process")
    lst_files = list(fd.askopenfilenames(filetypes=[("TIFF files", "*.tif *.tiff")]))

    if not lst_files:
        print("No files selected. Exiting.")
        return

    # Total output files = sum of cameras used across all activations, times all input files
    total_outputs = sum(len(act["cameras"]) for act in ACTIVATION_SCHEME) * len(lst_files)

    start_time = time.time()
    with Progress() as progress:
        task = progress.add_task("[green]Writing output files", total=total_outputs)

        for i, fpath in enumerate(lst_files, 1):
            fname = fpath.split("/")[-1]
            print(f"\n[{i}/{len(lst_files)}] {fname}")
            try:
                processor = ActivationProcessor(
                    fpath,
                    ACTIVATION_SCHEME,
                    max_concurrent_chunks=4,
                )
                processor.process(progress, task)

            except Exception as e:
                print(f"❌ Error processing {fpath}: {e}")

    elapsed = time.time() - start_time
    print(f"\n✅ All done in {elapsed:.1f}s")


if __name__ == "__main__":
    process_files()
