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
Parallel TIFF Splitter - Splits large TIFF stacks vertically into left and right halves.

Uses chunked parallel I/O and processing pipeline to handle large files efficiently
while managing memory usage based on available system RAM.
"""


class ChunkProcessor:
    def __init__(self, fpath, max_concurrent_chunks=4):
        self.fpath = fpath
        self.max_concurrent_chunks = max_concurrent_chunks
        self.chunk_queue = Queue(maxsize=max_concurrent_chunks + 1)
        self.results_queue = Queue()
        self.chunk_size = self._calculate_chunk_size()

    def _calculate_chunk_size(self):
        """Calculate chunk size allowing for multiple chunks in memory"""
        available_ram_gb = psutil.virtual_memory().available / (1024**3)
        usable_ram_gb = available_ram_gb * 0.8  # Use 80% of available RAM

        # Reserve memory for multiple chunks + processing overhead
        ram_per_chunk = usable_ram_gb / (self.max_concurrent_chunks + 1)
        frame_size_mb = 1.4
        frames_per_chunk = int((ram_per_chunk * 1024) / frame_size_mb)

        return max(500, min(frames_per_chunk, 3000))  # 500-3000 frames per chunk

    def load_chunk_worker(self, tif, chunk_info):
        """Worker function to load a single chunk"""
        chunk_id, start_idx, end_idx = chunk_info

        chunk_frames = []
        for j in range(start_idx, end_idx):
            chunk_frames.append(tif.pages[j].asarray())

        return chunk_id, chunk_frames

    def process_chunk_worker(self, chunk_data, halfwidth):
        """Worker function to process a single chunk"""
        chunk_id, chunk_frames = chunk_data

        processed_left = []
        processed_right = []

        for frame in chunk_frames:
            left_half = frame[:, :halfwidth]
            right_half = frame[:, halfwidth:]
            processed_left.append(left_half)
            processed_right.append(right_half)

        return chunk_id, processed_left, processed_right

    def process_file_parallel(self):
        with TiffFile(self.fpath) as tif:
            num_frames = len(tif.pages)
            first_frame = tif.pages[0].asarray()
            height, width = first_frame.shape
            halfwidth = width // 2

            # Create output paths
            left_path = self.fpath[:-4] + "-cropped-left.tif"
            right_path = self.fpath[:-4] + "-cropped-right.tif"

            # Calculate chunks
            chunks_info = []
            for i in range(0, num_frames, self.chunk_size):
                end_idx = min(i + self.chunk_size, num_frames)
                chunks_info.append((len(chunks_info), i, end_idx))


            # Accumulate all processed frames keyed by chunk_id
            all_chunks = {}  # chunk_id -> (left_frames, right_frames)

            with ThreadPoolExecutor(
                max_workers=self.max_concurrent_chunks,
                thread_name_prefix="IO-Loader",
            ) as io_executor, ThreadPoolExecutor(
                max_workers=min(4, psutil.cpu_count()),
                thread_name_prefix="Processor",
            ) as process_executor:

                with Progress() as progress:
                    load_task = progress.add_task(
                        "[cyan]Loading   ", total=len(chunks_info)
                    )
                    process_task = progress.add_task(
                        "[yellow]Processing", total=len(chunks_info)
                    )

                    # Submit all I/O jobs
                    io_futures = {}
                    for chunk_info in chunks_info:
                        future = io_executor.submit(
                            self.load_chunk_worker, tif, chunk_info
                        )
                        io_futures[future] = chunk_info[0]  # chunk_id

                    # As I/O completes, submit processing jobs
                    process_futures = {}
                    for completed_io in as_completed(io_futures):
                        try:
                            chunk_id, chunk_frames = completed_io.result()
                            progress.advance(load_task)

                            process_future = process_executor.submit(
                                self.process_chunk_worker,
                                (chunk_id, chunk_frames),
                                halfwidth,
                            )
                            process_futures[process_future] = chunk_id

                        except Exception as e:
                            print(f"Error loading chunk {io_futures[completed_io]}: {e}")

                    # Collect all processing results
                    for completed_process in as_completed(process_futures):
                        try:
                            p_chunk_id, left_frames, right_frames = completed_process.result()
                            all_chunks[p_chunk_id] = (left_frames, right_frames)
                            progress.advance(process_task)
                        except Exception as e:
                            print(f"Error in processing: {e}")

            # Assemble frames in order and write as ImageJ-compatible TIFF
            with Progress() as progress:
                write_task = progress.add_task("[green]Writing   ", total=4)

                all_left = []
                all_right = []
                for chunk_id in sorted(all_chunks.keys()):
                    left_frames, right_frames = all_chunks[chunk_id]
                    all_left.extend(left_frames)
                    all_right.extend(right_frames)
                    del all_chunks[chunk_id]
                    gc.collect()

                left_array = np.stack(all_left, axis=0)
                del all_left
                gc.collect()
                progress.advance(write_task)  # 1/4 stacked left

                imwrite(left_path, left_array, imagej=True)
                del left_array
                gc.collect()
                progress.advance(write_task)  # 2/4 wrote left

                right_array = np.stack(all_right, axis=0)
                del all_right
                gc.collect()
                progress.advance(write_task)  # 3/4 stacked right

                imwrite(right_path, right_array, imagej=True)
                del right_array
                gc.collect()
                progress.advance(write_task)  # 4/4 wrote right


# Main processing function
def process_files_parallel():
    print("Choose the tif files for crop")
    lst_files = list(fd.askopenfilenames())

    for i, fpath in enumerate(lst_files, 1):
        fname = fpath.split("/")[-1]
        print(f"\n[{i}/{len(lst_files)}] {fname}")
        try:
            start_time = time.time()

            processor = ChunkProcessor(fpath, max_concurrent_chunks=4)
            processor.process_file_parallel()

            end_time = time.time()
            print(f"✅ Done in {end_time - start_time:.1f}s")

        except Exception as e:
            print(f"❌ Error processing {fpath}: {str(e)}")


if __name__ == "__main__":
    process_files_parallel()
