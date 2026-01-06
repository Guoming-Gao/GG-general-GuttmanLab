from tifffile import TiffFile, TiffWriter
from tkinter import filedialog as fd
from rich.progress import Progress
from concurrent.futures import ThreadPoolExecutor, as_completed
from queue import Queue
import threading
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

    def process_chunk_worker(self, chunk_data, writers, halfwidth):
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

            # Create output files
            left_path = self.fpath[:-4] + "-cropped-left.tif"
            right_path = self.fpath[:-4] + "-cropped-right.tif"

            # Calculate chunks
            chunks_info = []
            for i in range(0, num_frames, self.chunk_size):
                end_idx = min(i + self.chunk_size, num_frames)
                chunks_info.append((len(chunks_info), i, end_idx))

            print(
                f"Processing {len(chunks_info)} chunks with {self.max_concurrent_chunks} parallel I/O operations"
            )

            with TiffWriter(left_path, bigtiff=True) as writer_left, TiffWriter(
                right_path, bigtiff=True
            ) as writer_right:

                # Phase 1: Parallel I/O Loading + Parallel Processing Pipeline
                with ThreadPoolExecutor(
                    max_workers=self.max_concurrent_chunks,
                    thread_name_prefix="IO-Loader",
                ) as io_executor, ThreadPoolExecutor(
                    max_workers=min(4, psutil.cpu_count()),
                    thread_name_prefix="Processor",
                ) as process_executor:

                    with Progress() as progress:
                        load_task = progress.add_task(
                            "Loading chunks", total=len(chunks_info)
                        )
                        process_task = progress.add_task(
                            "Processing chunks", total=len(chunks_info)
                        )
                        write_task = progress.add_task(
                            "Writing chunks", total=len(chunks_info)
                        )

                        # Submit I/O jobs
                        io_futures = {}
                        for chunk_info in chunks_info:
                            future = io_executor.submit(
                                self.load_chunk_worker, tif, chunk_info
                            )
                            io_futures[future] = chunk_info[0]  # chunk_id

                        # Process loaded chunks and maintain order
                        loaded_chunks = {}
                        processed_chunks = {}
                        next_write_chunk = 0
                        process_futures = {}

                        # Pipeline processing
                        for completed_io in as_completed(io_futures):
                            try:
                                chunk_id, chunk_frames = completed_io.result()
                                loaded_chunks[chunk_id] = chunk_frames
                                progress.advance(load_task)

                                # Submit for processing
                                process_future = process_executor.submit(
                                    self.process_chunk_worker,
                                    (chunk_id, chunk_frames),
                                    (writer_left, writer_right),
                                    halfwidth,
                                )
                                process_futures[process_future] = chunk_id

                                # Check for completed processing
                                completed_processes = []
                                for pf in list(process_futures.keys()):
                                    if pf.done():
                                        completed_processes.append(pf)

                                for completed_process in completed_processes:
                                    p_chunk_id, left_frames, right_frames = (
                                        completed_process.result()
                                    )
                                    processed_chunks[p_chunk_id] = (
                                        left_frames,
                                        right_frames,
                                    )
                                    progress.advance(process_task)
                                    del process_futures[completed_process]

                                    # Write chunks in order
                                    while next_write_chunk in processed_chunks:
                                        left_frames, right_frames = processed_chunks[
                                            next_write_chunk
                                        ]

                                        for left_frame, right_frame in zip(
                                            left_frames, right_frames
                                        ):
                                            writer_left.write(left_frame)
                                            writer_right.write(right_frame)

                                        progress.advance(write_task)

                                        # Cleanup memory
                                        del processed_chunks[next_write_chunk]
                                        if next_write_chunk in loaded_chunks:
                                            del loaded_chunks[next_write_chunk]

                                        next_write_chunk += 1
                                        gc.collect()

                            except Exception as e:
                                print(f"Error in chunk {io_futures[completed_io]}: {e}")

                        # Handle remaining processed chunks
                        for remaining_process in process_futures:
                            try:
                                p_chunk_id, left_frames, right_frames = (
                                    remaining_process.result()
                                )
                                processed_chunks[p_chunk_id] = (
                                    left_frames,
                                    right_frames,
                                )
                                progress.advance(process_task)
                            except Exception as e:
                                print(f"Error in remaining processing: {e}")

                        # Write remaining chunks in order
                        while next_write_chunk in processed_chunks:
                            left_frames, right_frames = processed_chunks[
                                next_write_chunk
                            ]

                            for left_frame, right_frame in zip(
                                left_frames, right_frames
                            ):
                                writer_left.write(left_frame)
                                writer_right.write(right_frame)

                            progress.advance(write_task)
                            del processed_chunks[next_write_chunk]
                            next_write_chunk += 1


# Main processing function
def process_files_parallel():
    print("Choose the tif files for crop")
    lst_files = list(fd.askopenfilenames())

    for i, fpath in enumerate(lst_files, 1):
        print(f"\n[{i}/{len(lst_files)}] Processing files...")
        try:
            print(f"\n🚀 Starting parallel processing: {fpath}")
            start_time = time.time()

            processor = ChunkProcessor(fpath, max_concurrent_chunks=4)
            processor.process_file_parallel()

            end_time = time.time()
            print(f"✅ Completed in {end_time - start_time:.1f} seconds: {fpath}")

        except Exception as e:
            print(f"❌ Error processing {fpath}: {str(e)}")


if __name__ == "__main__":
    process_files_parallel()
