from tifffile import TiffFile, TiffWriter
from tkinter import filedialog as fd
from rich.progress import Progress
import numpy as np
import psutil
import gc
import time


"""
TIFF Splitter - Splits large TIFF stacks vertically into left and right halves.

Uses sequential chunked I/O with streaming writes to handle arbitrarily large
files while keeping memory usage bounded based on available system RAM.

Architecture:
  For each chunk (sequential):
    1. Open a fresh TiffFile handle (avoids thread-safety issues)
    2. Load chunk_size frames sequentially
    3. Split each frame into left and right halves
    4. Append halves to output TiffWriters (streaming, no accumulation)
    5. Free chunk memory
"""


class ChunkProcessor:
    def __init__(self, fpath):
        self.fpath = fpath

    def _calculate_chunk_size(self, frame_bytes):
        """Calculate chunk size based on available RAM and actual frame size.

        We need memory for:
          - 1 chunk of full frames (loaded from disk)
          - 1 chunk of left halves (before writing)
          - 1 chunk of right halves (before writing)
        So peak memory ≈ 2× chunk of full frames (halves sum to 1 full).
        We target using at most 50% of available RAM.
        """
        available_ram = psutil.virtual_memory().available
        usable_ram = available_ram * 0.5  # Conservative: use 50% of available

        # Each frame held in memory: full + left_half + right_half = 2× full
        memory_per_frame = frame_bytes * 2
        frames_per_chunk = int(usable_ram / memory_per_frame)

        # Clamp to reasonable range: 100-3000 frames per chunk
        chunk_size = max(100, min(frames_per_chunk, 3000))
        return chunk_size

    def process_file(self):
        """Process the TIFF file: split into left and right halves with streaming writes."""

        # --- Phase 1: Read metadata (lightweight, no frame data loaded) ---
        with TiffFile(self.fpath) as tif:
            num_frames = len(tif.pages)
            page0 = tif.pages[0]
            height, width = page0.shape
            dtype = page0.dtype
            halfwidth = width // 2

        frame_bytes = height * width * np.dtype(dtype).itemsize
        chunk_size = self._calculate_chunk_size(frame_bytes)

        # Create output paths
        left_path = self.fpath[:-4] + "-cropped-left.tif"
        right_path = self.fpath[:-4] + "-cropped-right.tif"

        # Calculate chunk boundaries
        chunks = []
        for start in range(0, num_frames, chunk_size):
            end = min(start + chunk_size, num_frames)
            chunks.append((start, end))

        n_chunks = len(chunks)
        print(f"  {num_frames} frames, {height}×{width} {dtype}")
        print(f"  Chunk size: {chunk_size} frames, {n_chunks} chunks")
        print(f"  Available RAM: {psutil.virtual_memory().available / (1024**3):.1f} GB")

        # --- Phase 2: Sequential load → split → stream-write ---
        with TiffWriter(left_path, bigtiff=True) as writer_left, \
             TiffWriter(right_path, bigtiff=True) as writer_right:

            with Progress() as progress:
                task = progress.add_task(
                    "[cyan]Processing", total=num_frames
                )

                first_frame_written = True  # flag to write ImageJ metadata on first frame only

                for chunk_idx, (start, end) in enumerate(chunks):

                    # Load chunk: open a fresh TiffFile handle for each chunk
                    # to avoid any file-handle sharing issues
                    with TiffFile(self.fpath) as tif:
                        for j in range(start, end):
                            frame = tif.pages[j].asarray()
                            left_half = frame[:, :halfwidth]
                            right_half = frame[:, halfwidth:]

                            if first_frame_written:
                                # First frame: embed ImageJ-compatible TYX axes metadata
                                # so Fiji interprets the stack as a time-series
                                writer_left.write(
                                    left_half,
                                    contiguous=True,
                                    metadata={'axes': 'TYX'},
                                )
                                writer_right.write(
                                    right_half,
                                    contiguous=True,
                                    metadata={'axes': 'TYX'},
                                )
                                first_frame_written = False
                            else:
                                # Subsequent frames: contiguous=True reuses the same IFD
                                writer_left.write(
                                    left_half,
                                    contiguous=True,
                                    metadata=None,
                                )
                                writer_right.write(
                                    right_half,
                                    contiguous=True,
                                    metadata=None,
                                )

                            del frame, left_half, right_half
                            progress.advance(task)

                    gc.collect()

        # --- Phase 3: Verify output frame count ---
        with TiffFile(left_path) as tif:
            n_left = len(tif.pages)
        with TiffFile(right_path) as tif:
            n_right = len(tif.pages)

        if n_left != num_frames or n_right != num_frames:
            print(f"  ⚠️  Frame count mismatch! Expected {num_frames}, "
                  f"got left={n_left}, right={n_right}")
        else:
            print(f"  ✅ Verified: {n_left} left + {n_right} right frames written")


# Main processing function
def process_files():
    print("Choose the tif files for crop")
    lst_files = list(fd.askopenfilenames())

    for i, fpath in enumerate(lst_files, 1):
        fname = fpath.split("/")[-1]
        print(f"\n[{i}/{len(lst_files)}] {fname}")
        try:
            start_time = time.time()

            processor = ChunkProcessor(fpath)
            processor.process_file()

            end_time = time.time()
            print(f"  ✅ Done in {end_time - start_time:.1f}s")

        except Exception as e:
            print(f"  ❌ Error processing {fpath}: {str(e)}")


if __name__ == "__main__":
    process_files()
