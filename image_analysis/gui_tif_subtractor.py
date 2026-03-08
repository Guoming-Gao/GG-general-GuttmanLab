# tif_subtractor_gui.py

import tkinter as tk
from tkinter import filedialog, messagebox
import tifffile
import os


class TiffExtractor:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("TIFF Frame Extractor")
        self.root.geometry("400x200")

        self.file_path = None
        self.total_frames = 0

        # File selection
        tk.Button(self.root, text="Select TIFF File", command=self.select_file).pack(
            pady=10
        )
        self.file_label = tk.Label(self.root, text="No file selected")
        self.file_label.pack(pady=5)

        # Frame range input
        tk.Label(self.root, text="Frame Range (e.g., 1-1000):").pack(pady=5)
        self.range_entry = tk.Entry(self.root, width=20)
        self.range_entry.pack(pady=5)

        # Extract button
        tk.Button(self.root, text="Extract & Save", command=self.extract_frames).pack(
            pady=10
        )

    def select_file(self):
        self.file_path = filedialog.askopenfilename(
            filetypes=[("TIFF files", "*.tif *.tiff")]
        )
        if self.file_path:
            try:
                # Fast frame count - only reads TIFF directory, not image data
                with tifffile.TiffFile(self.file_path) as tif:
                    self.total_frames = len(tif.pages)
                self.file_label.config(
                    text=f"Selected: {os.path.basename(self.file_path)}\nTotal frames: {self.total_frames}"
                )
            except Exception as e:
                messagebox.showerror("Error", f"Cannot read TIFF: {e}")

    def extract_frames(self):
        if not self.file_path:
            messagebox.showwarning("Warning", "Please select a TIFF file first")
            return

        try:
            # Parse range
            range_text = self.range_entry.get().strip()
            start, end = map(int, range_text.split("-"))

            if start < 1 or end > self.total_frames or start > end:
                messagebox.showerror(
                    "Error", f"Invalid range. Use 1-{self.total_frames}"
                )
                return

            # Create output filename
            base_name = os.path.splitext(self.file_path)[0]
            output_path = f"{base_name}_frames_{start}-{end}.tif"

            # Extract frames efficiently - only reads specified range
            frames = tifffile.imread(
                self.file_path, key=slice(start - 1, end)
            )  # Convert to 0-based

            # Save extracted frames
            tifffile.imwrite(output_path, frames)

            messagebox.showinfo(
                "Success",
                f"Extracted {end-start+1} frames to:\n{os.path.basename(output_path)}",
            )

        except ValueError:
            messagebox.showerror("Error", "Invalid range format. Use: start-end")
        except Exception as e:
            messagebox.showerror("Error", f"Extraction failed: {e}")

    def run(self):
        self.root.mainloop()


if __name__ == "__main__":
    app = TiffExtractor()
    app.run()
