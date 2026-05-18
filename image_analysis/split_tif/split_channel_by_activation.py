from tifffile import TiffFile, TiffWriter
from tkinter import filedialog as fd
from rich.progress import Progress
import numpy as np
import psutil
import gc
import time


"""
Activation-Scheme TIFF Splitter
================================
Splits large TIFF stacks vertically into left and/or right camera halves,
guided by a user-defined activation scheme.

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

Output filenames are:   <input_stem>-<label>.tif
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
    {
        "frame_start": 21,
        "frame_end":   30,
        "cameras":     ["left"],
        "labels":      ["DAPI"],
    },
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


def _calculate_chunk_size(frame_bytes):
    """Calculate chunk size based on available RAM and actual frame size.

    Peak memory per activation = 2 × chunk (full frames + two halves).
    We target at most 50% of available RAM.
    """
    available_ram = psutil.virtual_memory().available
    usable_ram = available_ram * 0.5
    memory_per_frame = frame_bytes * 2   # full + left_half + right_half ≈ 2× full
    frames_per_chunk = int(usable_ram / memory_per_frame)
    return max(100, min(frames_per_chunk, 3000))


class ActivationProcessor:
    """
    Processes a single TIFF file according to an activation scheme.

    For each activation entry the processor:
      1. Reads the specified frame range.
      2. Splits each frame vertically into left / right halves.
      3. Keeps only the camera(s) specified for that activation.
      4. Writes one output TIFF per kept camera, named with the corresponding label.
    """

    def __init__(self, fpath, scheme):
        self.fpath = fpath
        self.scheme = scheme


    def _process_activation(self, activation, halfwidth, frame_bytes, act_index, total_acts, progress, task):
        """
        Load, split, and save frames for a single activation entry using
        sequential chunked I/O and streaming TiffWriter appends.

        Each chunk opens a fresh TiffFile handle to avoid file-handle sharing.
        Frames are written immediately after splitting — no full accumulation.
        Advances `task` in the shared `progress` bar once per output file written.
        """
        frame_start = activation["frame_start"] - 1   # convert to 0-based
        frame_end   = activation["frame_end"]          # exclusive upper bound (0-based)
        cameras     = activation["cameras"]
        labels      = activation["labels"]

        num_frames = frame_end - frame_start
        cam_to_label = dict(zip(cameras, labels))

        # Build output paths
        stem = self.fpath[:-4] if self.fpath.lower().endswith(".tif") else self.fpath
        out_paths = {
            cam: f"{stem}-{cam_to_label[cam]}.tif"
            for cam in cameras
        }

        chunk_size = _calculate_chunk_size(frame_bytes)

        print(
            f"  Activation {act_index}/{total_acts}: "
            f"frames {activation['frame_start']}-{activation['frame_end']} | "
            f"cameras: {cameras}"
        )
        for cam, path in out_paths.items():
            print(f"    → {cam}: {path.split('/')[-1]}")

        # ---- streaming write per camera ----------------------------------
        # Open one TiffWriter per camera output for this activation
        writers = {
            cam: TiffWriter(out_paths[cam], bigtiff=True)
            for cam in cameras
        }
        first_written = {cam: True for cam in cameras}

        try:
            # Sequential chunked I/O
            for chunk_start in range(0, num_frames, chunk_size):
                chunk_end = min(chunk_start + chunk_size, num_frames)
                abs_start = frame_start + chunk_start
                abs_end   = frame_start + chunk_end

                with TiffFile(self.fpath) as tif:
                    for j in range(abs_start, abs_end):
                        frame = tif.pages[j].asarray()

                        for cam in cameras:
                            half = frame[:, :halfwidth] if cam == "left" else frame[:, halfwidth:]

                            if first_written[cam]:
                                # First frame: embed ImageJ-compatible TYX axes metadata
                                writers[cam].write(
                                    half,
                                    contiguous=True,
                                    metadata={'axes': 'TYX'},
                                )
                                first_written[cam] = False
                            else:
                                writers[cam].write(
                                    half,
                                    contiguous=True,
                                    metadata=None,
                                )

                            del half

                        del frame

                gc.collect()

        finally:
            # Close all writers and advance progress bar
            for cam in cameras:
                writers[cam].close()
                progress.advance(task)   # one tick per output file written

        del writers
        gc.collect()

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def process(self, progress, task):
        """Process the file, reporting progress into the shared progress bar."""
        with TiffFile(self.fpath) as tif:
            num_frames = len(tif.pages)
            page0 = tif.pages[0]
            height, width = page0.shape
            halfwidth = width // 2
            dtype = page0.dtype

        frame_bytes = height * width * np.dtype(dtype).itemsize

        print(f"  File: {self.fpath.split('/')[-1]}")
        print(f"  Frames: {num_frames} | Resolution: {width}×{height} | Halfwidth: {halfwidth}")
        print(f"  Frame size: {frame_bytes / (1024**2):.2f} MB | Available RAM: {psutil.virtual_memory().available / (1024**3):.1f} GB")

        # Validate frame indices against actual file length
        for act in self.scheme:
            if act["frame_end"] > num_frames:
                raise ValueError(
                    f"Activation frame_end={act['frame_end']} exceeds "
                    f"total frames in file ({num_frames})."
                )

        total_acts = len(self.scheme)
        for idx, activation in enumerate(self.scheme, start=1):
            self._process_activation(activation, halfwidth, frame_bytes, idx, total_acts, progress, task)


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
                )
                processor.process(progress, task)

            except Exception as e:
                print(f"❌ Error processing {fpath}: {e}")

    elapsed = time.time() - start_time
    print(f"\n✅ All done in {elapsed:.1f}s")


if __name__ == "__main__":
    process_files()
