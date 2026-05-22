# tif_subtractor_gui.py

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import tifffile
import os
import threading
import queue


class TiffExtractor:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("TIFF Frame Extractor (Batch Mode)")
        self.root.geometry("600x550")
        self.root.minsize(500, 450)

        # Configure styles
        self.style = ttk.Style(self.root)
        self.style.theme_use("clam")

        self.file_paths = []
        self.queue = queue.Queue()
        self.processing = False

        # Main container with padding
        main_frame = ttk.Frame(self.root, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Title/Header
        title_label = ttk.Label(
            main_frame,
            text="TIFF Frame Extractor - Batch Process",
            font=("Helvetica", 14, "bold"),
        )
        title_label.pack(anchor=tk.W, pady=(0, 10))

        # --- File selection frame ---
        file_select_frame = ttk.Frame(main_frame)
        file_select_frame.pack(fill=tk.X, pady=(0, 5))

        self.select_btn = ttk.Button(
            file_select_frame, text="Select TIFF File(s)", command=self.select_files
        )
        self.select_btn.pack(side=tk.LEFT, padx=(0, 5))

        self.clear_btn = ttk.Button(
            file_select_frame, text="Clear Selection", command=self.clear_selection
        )
        self.clear_btn.pack(side=tk.LEFT)

        self.file_count_label = ttk.Label(file_select_frame, text="0 files selected")
        self.file_count_label.pack(side=tk.RIGHT)

        # --- Listbox for showing selected files ---
        list_frame = ttk.Frame(main_frame)
        list_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        list_frame.grid_columnconfigure(0, weight=1)
        list_frame.grid_rowconfigure(0, weight=1)

        self.file_listbox = tk.Listbox(
            list_frame,
            selectmode=tk.EXTENDED,
            font=("Courier", 10),
            bg="#ffffff",
            fg="#333333",
            relief=tk.SUNKEN,
            borderwidth=1,
        )
        self.file_listbox.grid(row=0, column=0, sticky="nsew")

        scrollbar_y = ttk.Scrollbar(
            list_frame, orient=tk.VERTICAL, command=self.file_listbox.yview
        )
        scrollbar_y.grid(row=0, column=1, sticky="ns")

        scrollbar_x = ttk.Scrollbar(
            list_frame, orient=tk.HORIZONTAL, command=self.file_listbox.xview
        )
        scrollbar_x.grid(row=1, column=0, sticky="ew")

        self.file_listbox.configure(
            yscrollcommand=scrollbar_y.set, xscrollcommand=scrollbar_x.set
        )

        # --- Frame range input and Extraction frame ---
        action_frame = ttk.Frame(main_frame)
        action_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(
            action_frame, text="Frame Range (e.g., 1-1000):", font=("Helvetica", 10)
        ).pack(side=tk.LEFT, padx=(0, 5))
        self.range_entry = ttk.Entry(action_frame, width=15)
        self.range_entry.pack(side=tk.LEFT, padx=(0, 10))

        self.extract_btn = ttk.Button(
            action_frame, text="Extract & Save All", command=self.start_extraction
        )
        self.extract_btn.pack(side=tk.RIGHT)

        # --- Log frame ---
        log_frame = ttk.Frame(main_frame)
        log_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        log_frame.grid_columnconfigure(0, weight=1)
        log_frame.grid_rowconfigure(1, weight=1)

        ttk.Label(log_frame, text="Status Log:", font=("Helvetica", 10, "bold")).grid(
            row=0, column=0, sticky="w", pady=(0, 5)
        )

        self.log_text = tk.Text(
            log_frame,
            height=8,
            state=tk.DISABLED,
            font=("Courier", 9),
            bg="#f9f9f9",
            fg="#222222",
            relief=tk.SUNKEN,
            borderwidth=1,
            wrap=tk.WORD,
        )
        self.log_text.grid(row=1, column=0, sticky="nsew")

        log_scroll = ttk.Scrollbar(
            log_frame, orient=tk.VERTICAL, command=self.log_text.yview
        )
        log_scroll.grid(row=1, column=1, sticky="ns")
        self.log_text.configure(yscrollcommand=log_scroll.set)

        # --- Status Bar at the bottom ---
        self.status_label = ttk.Label(
            self.root,
            text="Ready",
            relief=tk.SUNKEN,
            anchor=tk.W,
            padding=(5, 3),
            font=("Helvetica", 10),
        )
        self.status_label.pack(fill=tk.X, side=tk.BOTTOM)

        # Start background polling
        self.poll_queue()

    def select_files(self):
        if self.processing:
            return
        paths = filedialog.askopenfilenames(
            title="Select TIFF File(s)",
            filetypes=[("TIFF files", "*.tif *.tiff")],
        )
        if paths:
            self.file_paths = list(paths)
            self.update_file_listbox()
            self.log_message(f"Selected {len(self.file_paths)} files.")

    def clear_selection(self):
        if self.processing:
            return
        self.file_paths = []
        self.update_file_listbox()
        self.log_message("Selection cleared.")
        self.status_label.config(text="Ready")

    def update_file_listbox(self):
        self.file_listbox.delete(0, tk.END)
        for path in self.file_paths:
            self.file_listbox.insert(tk.END, os.path.basename(path))
        self.file_count_label.config(text=f"{len(self.file_paths)} file(s) selected")

    def log_message(self, message):
        self.queue.put(("log", message))

    def write_log(self, message):
        self.log_text.config(state=tk.NORMAL)
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)
        self.log_text.config(state=tk.DISABLED)

    def poll_queue(self):
        try:
            while True:
                msg_type, data = self.queue.get_nowait()
                if msg_type == "log":
                    self.write_log(data)
                elif msg_type == "status":
                    self.status_label.config(text=data)
                elif msg_type == "finished":
                    self.processing_finished(data)
        except queue.Empty:
            pass
        self.root.after(100, self.poll_queue)

    def start_extraction(self):
        if self.processing:
            return
        if not self.file_paths:
            self.status_label.config(text="Warning: No files selected")
            self.write_log("[Warning] Please select TIFF files first.")
            return

        range_text = self.range_entry.get().strip()
        if not range_text:
            self.status_label.config(text="Warning: Frame range is empty")
            self.write_log(
                "[Warning] Please specify a frame range (e.g., 1-1000)."
            )
            return

        try:
            start, end = map(int, range_text.split("-"))
            if start < 1 or start > end:
                raise ValueError
        except ValueError:
            self.status_label.config(text="Error: Invalid range format")
            self.write_log(
                "[Error] Invalid range format. Please use: start-end (e.g., 1-1000) where start >= 1 and start <= end."
            )
            return

        self.processing = True
        self.select_btn.config(state=tk.DISABLED)
        self.clear_btn.config(state=tk.DISABLED)
        self.extract_btn.config(state=tk.DISABLED)
        self.range_entry.config(state=tk.DISABLED)

        self.status_label.config(text="Processing...")
        self.write_log(f"\n--- Starting extraction for range {start}-{end} ---")

        # Start worker thread
        thread = threading.Thread(
            target=self.extract_worker,
            args=(self.file_paths.copy(), start, end),
            daemon=True,
        )
        thread.start()

    def extract_worker(self, file_paths, start, end):
        success_count = 0
        failure_count = 0
        total_files = len(file_paths)

        for idx, file_path in enumerate(file_paths, 1):
            base_name = os.path.basename(file_path)
            self.queue.put(
                ("status", f"Processing {idx}/{total_files}: {base_name}")
            )
            self.queue.put(
                ("log", f"[{idx}/{total_files}] Processing '{base_name}'...")
            )

            try:
                # Read total frames
                with tifffile.TiffFile(file_path) as tif:
                    total_frames = len(tif.pages)

                # Validate range for this specific file
                if end > total_frames:
                    self.queue.put(
                        (
                            "log",
                            f"  [Error] '{base_name}' has only {total_frames} frames. "
                            f"Requested range {start}-{end} exceeds total frames.",
                        )
                    )
                    failure_count += 1
                    continue

                # Perform extraction
                frames = tifffile.imread(file_path, key=slice(start - 1, end))

                # Output path in the same folder
                dir_name = os.path.dirname(file_path)
                file_base = os.path.splitext(base_name)[0]
                output_name = f"{file_base}_frames_{start}-{end}.tif"
                output_path = os.path.join(dir_name, output_name)

                tifffile.imwrite(output_path, frames)

                self.queue.put(("log", f"  [Success] Saved to '{output_name}'"))
                success_count += 1

            except Exception as e:
                self.queue.put(
                    ("log", f"  [Error] Failed to process '{base_name}': {e}")
                )
                failure_count += 1

        summary = f"Finished: {success_count} succeeded, {failure_count} failed."
        self.queue.put(("finished", (success_count, failure_count, summary)))

    def processing_finished(self, summary_data):
        success_count, failure_count, summary = summary_data
        self.processing = False

        self.select_btn.config(state=tk.NORMAL)
        self.clear_btn.config(state=tk.NORMAL)
        self.extract_btn.config(state=tk.NORMAL)
        self.range_entry.config(state=tk.NORMAL)

        self.status_label.config(text=summary)
        self.write_log(f"\n--- Batch processing finished ---\n{summary}\n")

    def run(self):
        self.root.mainloop()


if __name__ == "__main__":
    app = TiffExtractor()
    app.run()
