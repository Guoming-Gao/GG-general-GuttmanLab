# tif_montage_gui_composite.py

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import tifffile
from PIL import Image, ImageTk
import numpy as np
import os
import math


class TifMontageGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Multi-Channel TIF Montage Creator")
        self.root.geometry("700x650")
        self.selected_files = []
        self.save_folder = ""
        self.save_filename = ""
        self.create_widgets()

    def create_widgets(self):
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Section 1: File Selection
        file_frame = ttk.LabelFrame(
            main_frame, text="1. Select Multi-Channel TIF Files", padding="10"
        )
        file_frame.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=(0, 10))

        ttk.Button(file_frame, text="Choose Files", command=self.select_files).grid(
            row=0, column=0
        )

        # Section 2: File List and Grid Layout
        list_frame = ttk.LabelFrame(
            main_frame, text="2. Selected Files & Layout", padding="10"
        )
        list_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(0, 10))

        # File list
        list_subframe = ttk.Frame(list_frame)
        list_subframe.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        self.file_listbox = tk.Listbox(list_subframe, height=10)
        self.file_listbox.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        scrollbar = ttk.Scrollbar(
            list_subframe, orient="vertical", command=self.file_listbox.yview
        )
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.file_listbox.configure(yscrollcommand=scrollbar.set)

        # Grid info
        grid_info_frame = ttk.Frame(list_frame)
        grid_info_frame.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=(10, 0))

        ttk.Label(grid_info_frame, text="Layout:").grid(row=0, column=0)
        self.layout_label = ttk.Label(
            grid_info_frame, text="No files selected", foreground="blue"
        )
        self.layout_label.grid(row=0, column=1, padx=(10, 0))

        # Channel info
        ttk.Label(grid_info_frame, text="Channels:").grid(row=1, column=0)
        self.channel_label = ttk.Label(grid_info_frame, text="", foreground="green")
        self.channel_label.grid(row=1, column=1, padx=(10, 0))

        # Section 3: Save Configuration
        save_frame = ttk.LabelFrame(
            main_frame, text="3. Save Configuration", padding="10"
        )
        save_frame.grid(row=2, column=0, sticky=(tk.W, tk.E), pady=(0, 10))

        # Folder selection
        folder_frame = ttk.Frame(save_frame)
        folder_frame.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=(0, 10))

        ttk.Label(folder_frame, text="Save Folder:").grid(row=0, column=0, sticky=tk.W)
        self.folder_label = ttk.Label(
            folder_frame, text="No folder selected", foreground="gray", wraplength=400
        )
        self.folder_label.grid(
            row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(5, 0)
        )

        ttk.Button(folder_frame, text="Browse Folder", command=self.select_folder).grid(
            row=0, column=1, padx=(10, 0)
        )

        # Filename input
        filename_frame = ttk.Frame(save_frame)
        filename_frame.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=(10, 0))

        ttk.Label(filename_frame, text="Filename:").grid(row=0, column=0, sticky=tk.W)

        filename_input_frame = ttk.Frame(filename_frame)
        filename_input_frame.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=(5, 0))

        self.filename_entry = ttk.Entry(filename_input_frame, width=40)
        self.filename_entry.grid(row=0, column=0, sticky=(tk.W, tk.E))
        ttk.Label(filename_input_frame, text=".tif").grid(row=0, column=1, padx=(5, 0))

        # Section 4: Execute
        execute_frame = ttk.Frame(main_frame)
        execute_frame.grid(row=3, column=0, pady=(10, 0))

        self.execute_btn = ttk.Button(
            execute_frame,
            text="Create Multi-Channel Montage",
            command=self.create_montage,
            state="disabled",
        )
        self.execute_btn.grid(row=0, column=0)

        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(1, weight=1)
        list_frame.columnconfigure(0, weight=1)
        list_frame.rowconfigure(0, weight=1)
        list_subframe.columnconfigure(0, weight=1)
        list_subframe.rowconfigure(0, weight=1)
        folder_frame.columnconfigure(0, weight=1)
        filename_frame.columnconfigure(0, weight=1)
        filename_input_frame.columnconfigure(0, weight=1)

    def select_files(self):
        filetypes = [("TIF files", "*.tif *.tiff"), ("All files", "*.*")]
        files = filedialog.askopenfilenames(
            title="Select Multi-Channel TIF files", filetypes=filetypes
        )
        if files:
            self.selected_files = list(files)
            self.update_file_list()
            self.update_layout_info()
            self.update_channel_info()
            self.update_execute_button()

    def update_file_list(self):
        self.file_listbox.delete(0, tk.END)
        for file_path in self.selected_files:
            filename = os.path.basename(file_path)
            self.file_listbox.insert(tk.END, filename)

    def calculate_optimal_grid(self, n_files):
        if n_files == 0:
            return 0, 0
        if n_files == 1:
            return 1, 1

        best_cols, best_rows = 1, n_files
        best_score = float("inf")

        max_cols = min(n_files, int(math.sqrt(n_files)) + 3)

        for cols in range(1, max_cols + 1):
            rows = math.ceil(n_files / cols)
            aspect_ratio = max(cols, rows) / min(cols, rows)

            if aspect_ratio > 4:
                continue

            empty_slots = rows * cols - n_files
            square_penalty = abs(cols - rows) * 0.5
            empty_penalty = empty_slots * 0.3
            aspect_penalty = aspect_ratio * 0.2

            score = square_penalty + empty_penalty + aspect_penalty

            if score < best_score:
                best_cols, best_rows = cols, rows
                best_score = score

        return best_cols, best_rows

    def update_layout_info(self):
        if not self.selected_files:
            self.layout_label.config(text="No files selected")
            return

        n_files = len(self.selected_files)
        cols, rows = self.calculate_optimal_grid(n_files)
        empty_slots = rows * cols - n_files

        layout_text = f"{cols} columns × {rows} rows"
        if empty_slots > 0:
            layout_text += (
                f" ({empty_slots} empty slot{'s' if empty_slots > 1 else ''})"
            )

        self.layout_label.config(text=layout_text)

    def update_channel_info(self):
        if not self.selected_files:
            self.channel_label.config(text="")
            return

        try:
            channel_counts = []
            for file_path in self.selected_files:
                with tifffile.TiffFile(file_path) as tif:
                    # Read metadata to get shape
                    page = tif.pages[0]
                    shape = page.shape
                    # For CYX format, channels are first dimension
                    if len(shape) == 3:  # Multi-channel
                        channels = shape[0]
                    else:  # Single channel
                        channels = 1
                    channel_counts.append(channels)

            unique_counts = list(set(channel_counts))
            if len(unique_counts) == 1:
                self.channel_label.config(
                    text=f"All files: {unique_counts[0]} channels"
                )
            else:
                min_ch, max_ch = min(channel_counts), max(channel_counts)
                self.channel_label.config(text=f"Mixed: {min_ch}-{max_ch} channels")

        except Exception as e:
            self.channel_label.config(text="Error reading channel info")

    def select_folder(self):
        folder = filedialog.askdirectory(title="Select save folder")
        if folder:
            self.save_folder = folder
            display_path = folder if len(folder) <= 50 else "..." + folder[-47:]
            self.folder_label.config(text=display_path, foreground="black")
            self.update_execute_button()

    def update_execute_button(self):
        if (
            self.selected_files
            and self.save_folder
            and self.filename_entry.get().strip()
        ):
            self.execute_btn.config(state="normal")
        else:
            self.execute_btn.config(state="disabled")

    def create_montage(self):
        if not self.selected_files or not self.save_folder:
            messagebox.showerror("Error", "Please select files and save folder")
            return

        filename = self.filename_entry.get().strip()
        if not filename:
            messagebox.showerror("Error", "Please enter a filename")
            return

        try:
            # Load all images and determine maximum dimensions and channel count
            all_images = []
            max_height, max_width = 0, 0
            max_channels = 0

            for file_path in self.selected_files:
                img_array = tifffile.imread(file_path)

                # Handle different formats - ensure CYX format
                if len(img_array.shape) == 2:  # Single channel (YX)
                    height, width = img_array.shape
                    channels = 1
                    # Convert to CYX format
                    img_array = img_array[np.newaxis, :, :]
                elif len(img_array.shape) == 3:  # Multi-channel (CYX)
                    channels, height, width = img_array.shape
                else:
                    raise ValueError(f"Unsupported image format in {file_path}")

                max_height = max(max_height, height)
                max_width = max(max_width, width)
                max_channels = max(max_channels, channels)
                all_images.append(img_array)

            # Calculate grid layout
            cols, rows = self.calculate_optimal_grid(len(all_images))

            # Create montage for each channel (Option A)
            montage_height = max_height * rows
            montage_width = max_width * cols

            # Initialize the final montage with maximum channels
            final_montage = np.zeros(
                (max_channels, montage_height, montage_width), dtype=all_images[0].dtype
            )

            # Process each channel separately
            for channel_idx in range(max_channels):
                for img_idx, img in enumerate(all_images):
                    # Calculate position in grid
                    row = img_idx // cols
                    col = img_idx % cols
                    y_start = row * max_height
                    x_start = col * max_width

                    # Get image dimensions
                    if img.shape[0] > channel_idx:  # This image has this channel
                        channel_img = img[channel_idx]
                        h, w = channel_img.shape

                        # Center the image in the allocated space
                        y_offset = (max_height - h) // 2
                        x_offset = (max_width - w) // 2
                        y_end = y_start + y_offset + h
                        x_end = x_start + x_offset + w

                        final_montage[
                            channel_idx,
                            y_start + y_offset : y_end,
                            x_start + x_offset : x_end,
                        ] = channel_img
                    # If image doesn't have this channel, leave as zeros (black)

            # Save the montage with proper metadata
            save_path = os.path.join(self.save_folder, filename + ".tif")
            tifffile.imwrite(save_path, final_montage, metadata={"axes": "CYX"})

            messagebox.showinfo(
                "Success", f"Multi-channel montage saved to:\n{save_path}"
            )

            # Clear selected files but keep save settings
            self.selected_files = []
            self.update_file_list()
            self.update_layout_info()
            self.update_channel_info()
            self.update_execute_button()

        except Exception as e:
            messagebox.showerror("Error", f"Failed to create montage:\n{str(e)}")


def main():
    root = tk.Tk()
    app = TifMontageGUI(root)
    app.filename_entry.bind("<KeyRelease>", lambda e: app.update_execute_button())
    root.mainloop()


if __name__ == "__main__":
    main()
