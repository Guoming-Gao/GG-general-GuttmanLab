#!/usr/bin/env python3
"""
CellposeSAM_result_viewer_gui.py

GUI application for viewing and analyzing CellposeSAM segmentation results.
Displays nucleus statistics and plots cropped images with mask overlays.

Usage:
    python CellposeSAM_result_viewer_gui.py

Author: Generated for Cellpose-SAM result visualization
Date: 2025-11-05
"""

import numpy as np
import pandas as pd
from pathlib import Path
import pickle
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import threading
import queue
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.lines import Line2D
import tifffile
import matplotlib

matplotlib.use("TkAgg")


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def calculate_channel_stats(image_channel, mask):
    """
    Calculate intensity statistics for a single channel within a mask.

    Parameters:
    -----------
    image_channel : 2D array
        Single channel image
    mask : 2D boolean array
        Binary mask

    Returns:
    --------
    dict : Statistics (mean, max)
    """
    if np.sum(mask) == 0:
        return {"mean": 0.0, "max": 0}

    pixels = image_channel[mask]
    return {"mean": float(np.mean(pixels)), "max": int(np.max(pixels))}


# =============================================================================
# MAIN GUI CLASS
# =============================================================================


class CellposeSAMViewer:
    def __init__(self, root):
        self.root = root
        self.root.title("CellposeSAM Result Viewer")
        self.root.geometry("1500x720")

        # Data storage
        self.selected_pkl_files = []
        self.current_pkl_path = None
        self.current_pkl_data = None
        self.nucleus_stats_df = None
        self.nucleus_stats_df_original = None  # Keep original unsorted copy
        self.nucleus_data_dict = None
        self.n_channels = 0

        # Progress tracking
        self.analysis_queue = queue.Queue()
        self.analyzing = False

        # Sorting state
        self.sort_columns = []  # List of (column, ascending) tuples

        # Create GUI
        self.create_widgets()

        # Start queue polling
        self.poll_queue()

    def create_widgets(self):
        """Build all GUI components."""

        # Main container
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        row = 0

        # ===== FILE SELECTION =====
        ttk.Label(main_frame, text="Pickle Files", font=("Arial", 12, "bold")).grid(
            row=row, column=0, sticky=tk.W, pady=(0, 5)
        )

        ttk.Button(main_frame, text="Select Files...", command=self.select_files).grid(
            row=row, column=1, sticky=tk.E, pady=(0, 5)
        )
        row += 1

        # File listbox
        list_frame = ttk.Frame(main_frame)
        list_frame.grid(
            row=row,
            column=0,
            columnspan=2,
            sticky=(tk.W, tk.E, tk.N, tk.S),
            pady=(0, 10),
        )

        self.file_listbox = tk.Listbox(list_frame, height=5, selectmode=tk.SINGLE)
        scrollbar = ttk.Scrollbar(
            list_frame, orient=tk.VERTICAL, command=self.file_listbox.yview
        )
        self.file_listbox.configure(yscrollcommand=scrollbar.set)
        self.file_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        row += 1

        # Selection info and buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(
            row=row, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10)
        )

        self.selected_file_label = ttk.Label(
            button_frame, text="No file selected", foreground="gray"
        )
        self.selected_file_label.pack(side=tk.LEFT)

        ttk.Button(button_frame, text="View", command=self.view_file, width=12).pack(
            side=tk.RIGHT, padx=5
        )
        ttk.Button(button_frame, text="Clear", command=self.clear_files, width=12).pack(
            side=tk.RIGHT
        )
        row += 1

        # ===== PROGRESS BAR =====
        self.progress_frame = ttk.Frame(main_frame)
        self.progress_frame.grid(
            row=row, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10)
        )

        self.progress_label = ttk.Label(self.progress_frame, text="")
        self.progress_label.pack(side=tk.LEFT)

        self.progress_bar = ttk.Progressbar(
            self.progress_frame, mode="determinate", length=600
        )
        self.progress_bar.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(10, 0))

        self.progress_frame.grid_remove()  # Hidden by default
        row += 1

        # ===== NUCLEUS TABLE =====
        table_header_frame = ttk.Frame(main_frame)
        table_header_frame.grid(
            row=row, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(10, 5)
        )

        ttk.Label(
            table_header_frame, text="Nucleus Table", font=("Arial", 12, "bold")
        ).pack(side=tk.LEFT)

        ttk.Button(
            table_header_frame, text="Export to CSV", command=self.export_csv
        ).pack(side=tk.RIGHT)
        row += 1

        # Table with scrollbars
        table_frame = ttk.Frame(main_frame)
        table_frame.grid(
            row=row,
            column=0,
            columnspan=2,
            sticky=(tk.W, tk.E, tk.N, tk.S),
            pady=(0, 5),
        )

        # Create treeview with scrollbars
        self.tree_scroll_y = ttk.Scrollbar(table_frame, orient=tk.VERTICAL)
        self.tree_scroll_x = ttk.Scrollbar(table_frame, orient=tk.HORIZONTAL)

        self.nucleus_tree = ttk.Treeview(
            table_frame,
            yscrollcommand=self.tree_scroll_y.set,
            xscrollcommand=self.tree_scroll_x.set,
            selectmode="browse",  # Single select
            height=20,  # Changed from 15 to 20
        )

        self.tree_scroll_y.config(command=self.nucleus_tree.yview)
        self.tree_scroll_x.config(command=self.nucleus_tree.xview)

        self.nucleus_tree.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.tree_scroll_y.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.tree_scroll_x.grid(row=1, column=0, sticky=(tk.W, tk.E))

        table_frame.rowconfigure(0, weight=1)
        table_frame.columnconfigure(0, weight=1)

        # Bind selection event
        self.nucleus_tree.bind("<<TreeviewSelect>>", self.on_nucleus_select)
        row += 1

        # ===== SORTING CONTROLS =====
        sort_frame = ttk.Frame(main_frame)
        sort_frame.grid(
            row=row, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(5, 10)
        )

        ttk.Label(sort_frame, text="Sort by:", font=("Arial", 10)).pack(
            side=tk.LEFT, padx=(0, 5)
        )

        # First sort column
        self.sort_col1 = ttk.Combobox(sort_frame, width=18, state="readonly")
        self.sort_col1.pack(side=tk.LEFT, padx=2)

        self.sort_dir1 = ttk.Combobox(
            sort_frame,
            values=["High to Low", "Low to High"],
            width=12,
            state="readonly",
        )
        self.sort_dir1.set("High to Low")
        self.sort_dir1.pack(side=tk.LEFT, padx=2)

        ttk.Label(sort_frame, text="then", font=("Arial", 9, "italic")).pack(
            side=tk.LEFT, padx=5
        )

        # Second sort column
        self.sort_col2 = ttk.Combobox(sort_frame, width=18, state="readonly")
        self.sort_col2.pack(side=tk.LEFT, padx=2)

        self.sort_dir2 = ttk.Combobox(
            sort_frame,
            values=["High to Low", "Low to High"],
            width=12,
            state="readonly",
        )
        self.sort_dir2.set("High to Low")
        self.sort_dir2.pack(side=tk.LEFT, padx=2)

        ttk.Label(sort_frame, text="then", font=("Arial", 9, "italic")).pack(
            side=tk.LEFT, padx=5
        )

        # Third sort column
        self.sort_col3 = ttk.Combobox(sort_frame, width=18, state="readonly")
        self.sort_col3.pack(side=tk.LEFT, padx=2)

        self.sort_dir3 = ttk.Combobox(
            sort_frame,
            values=["High to Low", "Low to High"],
            width=12,
            state="readonly",
        )
        self.sort_dir3.set("High to Low")
        self.sort_dir3.pack(side=tk.LEFT, padx=2)

        # Apply sort button
        ttk.Button(
            sort_frame, text="Apply Sort", command=self.apply_sort, width=12
        ).pack(side=tk.LEFT, padx=10)
        ttk.Button(sort_frame, text="Reset", command=self.reset_sort, width=8).pack(
            side=tk.LEFT, padx=2
        )

        row += 1

        # Selection info and plot button
        plot_frame = ttk.Frame(main_frame)
        plot_frame.grid(
            row=row, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10)
        )

        self.selection_label = ttk.Label(
            plot_frame, text="No nucleus selected", foreground="gray"
        )
        self.selection_label.pack(side=tk.LEFT)

        self.plot_button = ttk.Button(
            plot_frame,
            text="Plot Nucleus",
            command=self.plot_nucleus,
            state=tk.DISABLED,
            width=15,
        )
        self.plot_button.pack(side=tk.RIGHT)
        row += 1

        # Configure grid weights
        main_frame.rowconfigure(2, weight=0)  # Listbox - small
        main_frame.rowconfigure(6, weight=1)  # Table - expandable
        main_frame.columnconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=0)

        self.root.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)

    def select_files(self):
        """Open file selection dialog for pickle files."""
        files = filedialog.askopenfilenames(
            title="Select pickle files",
            filetypes=[("Pickle files", "*.pkl"), ("All files", "*.*")],
        )

        if files:
            self.selected_pkl_files = list(files)
            self.file_listbox.delete(0, tk.END)
            for f in self.selected_pkl_files:
                self.file_listbox.insert(tk.END, Path(f).name)

    def clear_files(self):
        """Clear file list."""
        self.selected_pkl_files = []
        self.file_listbox.delete(0, tk.END)
        self.selected_file_label.configure(text="No file selected", foreground="gray")
        self.clear_table()
        self.clear_sort_controls()

    def clear_sort_controls(self):
        """Clear sorting dropdowns."""
        self.sort_col1.set("")
        self.sort_col2.set("")
        self.sort_col3.set("")
        self.sort_col1["values"] = []
        self.sort_col2["values"] = []
        self.sort_col3["values"] = []

    def view_file(self):
        """View selected pickle file."""
        selection = self.file_listbox.curselection()
        if not selection:
            messagebox.showerror("Error", "Please select a pickle file to view")
            return

        idx = selection[0]
        pkl_path = Path(self.selected_pkl_files[idx])

        self.selected_file_label.configure(text=pkl_path.name, foreground="black")

        # Load and analyze in thread
        self.analyzing = True
        thread = threading.Thread(
            target=self.analyze_pkl_worker, args=(pkl_path,), daemon=True
        )
        thread.start()

        # Show progress
        self.progress_frame.grid()
        self.progress_label.configure(text="Loading pickle file...")
        self.progress_bar["mode"] = "indeterminate"
        self.progress_bar.start()

    def analyze_pkl_worker(self, pkl_path):
        """Analyze pickle file in worker thread."""
        try:
            # Load pickle
            with open(pkl_path, "rb") as f:
                pkl_data = pickle.load(f)

            self.analysis_queue.put(
                {"type": "loaded", "pkl_data": pkl_data, "pkl_path": pkl_path}
            )

            # Extract data
            nucleus_data = pkl_data["nucleus_data"]
            n_nuclei = len(nucleus_data)

            # Determine number of channels from first nucleus
            first_nuc = list(nucleus_data.values())[0]
            n_channels = first_nuc["image_crop"].shape[0]

            self.analysis_queue.put(
                {
                    "type": "start_analysis",
                    "n_nuclei": n_nuclei,
                    "n_channels": n_channels,
                }
            )

            # Calculate statistics for each nucleus
            stats_list = []

            for i, (nuc_id, nuc_data) in enumerate(nucleus_data.items()):
                # Extract data
                image_crop = nuc_data["image_crop"]
                nuc_mask = nuc_data["nucleus_mask_crop"]
                cell_mask = nuc_data["cell_mask_crop"]

                # Calculate cytoplasm mask
                cyto_mask = cell_mask & ~nuc_mask

                # Basic info - following exact column order requested
                row_data = {
                    "Nucleus_ID": nuc_id,
                    "Matched": "✓" if nuc_data["is_matched"] else "✗",
                    "Nucleus_Area": nuc_data["nucleus_area"],
                    "Cyto_Area": int(np.sum(cyto_mask)),
                }

                # Calculate stats for each channel in specific order:
                # Ch{N}_Nuc_Mean, Ch{N}_Cyto_Mean, Ch{N}_Nuc_Max, Ch{N}_Cyto_Max
                for ch in range(n_channels):
                    ch_num = ch + 1
                    ch_img = image_crop[ch]

                    # Nucleus stats
                    nuc_stats = calculate_channel_stats(ch_img, nuc_mask)

                    # Cyto stats
                    cyto_stats = calculate_channel_stats(ch_img, cyto_mask)

                    # Add in specified order
                    row_data[f"Ch{ch_num}_Nuc_Mean"] = nuc_stats["mean"]
                    row_data[f"Ch{ch_num}_Cyto_Mean"] = cyto_stats["mean"]
                    row_data[f"Ch{ch_num}_Nuc_Max"] = nuc_stats["max"]
                    row_data[f"Ch{ch_num}_Cyto_Max"] = cyto_stats["max"]

                stats_list.append(row_data)

                # Update progress
                self.analysis_queue.put(
                    {"type": "progress", "current": i + 1, "total": n_nuclei}
                )

            # Create DataFrame with explicit column order
            columns = ["Nucleus_ID", "Matched", "Nucleus_Area", "Cyto_Area"]
            for ch in range(n_channels):
                ch_num = ch + 1
                columns.extend(
                    [
                        f"Ch{ch_num}_Nuc_Mean",
                        f"Ch{ch_num}_Cyto_Mean",
                        f"Ch{ch_num}_Nuc_Max",
                        f"Ch{ch_num}_Cyto_Max",
                    ]
                )

            df = pd.DataFrame(stats_list, columns=columns)

            self.analysis_queue.put(
                {"type": "complete", "df": df, "nucleus_data": nucleus_data}
            )

        except Exception as e:
            import traceback

            self.analysis_queue.put(
                {"type": "error", "error": str(e), "traceback": traceback.format_exc()}
            )

    def poll_queue(self):
        """Poll analysis queue for updates."""
        try:
            while True:
                item = self.analysis_queue.get_nowait()

                if item["type"] == "loaded":
                    self.current_pkl_data = item["pkl_data"]
                    self.current_pkl_path = item["pkl_path"]

                elif item["type"] == "start_analysis":
                    self.n_channels = item["n_channels"]
                    self.progress_bar.stop()
                    self.progress_bar["mode"] = "determinate"
                    self.progress_bar["maximum"] = item["n_nuclei"]
                    self.progress_bar["value"] = 0

                elif item["type"] == "progress":
                    self.progress_bar["value"] = item["current"]
                    self.progress_label.configure(
                        text=f"Analyzing nucleus {item['current']}/{item['total']}..."
                    )

                elif item["type"] == "complete":
                    self.nucleus_stats_df = item["df"].copy()
                    self.nucleus_stats_df_original = item["df"].copy()  # Keep original
                    self.nucleus_data_dict = item["nucleus_data"]
                    self.populate_table()
                    self.update_sort_dropdowns()
                    self.progress_frame.grid_remove()
                    self.analyzing = False

                elif item["type"] == "error":
                    messagebox.showerror("Error", f"Analysis failed:\n{item['error']}")
                    self.progress_frame.grid_remove()
                    self.analyzing = False

        except queue.Empty:
            pass

        self.root.after(100, self.poll_queue)

    def update_sort_dropdowns(self):
        """Update sort dropdown menus with column names."""
        if self.nucleus_stats_df is None:
            return

        columns = list(self.nucleus_stats_df.columns)

        self.sort_col1["values"] = [""] + columns
        self.sort_col2["values"] = [""] + columns
        self.sort_col3["values"] = [""] + columns

    def apply_sort(self):
        """Apply multi-column sorting."""
        if self.nucleus_stats_df_original is None:
            return

        # Build list of sort criteria
        sort_criteria = []
        ascending_list = []

        col1 = self.sort_col1.get()
        if col1 and col1 in self.nucleus_stats_df_original.columns:
            sort_criteria.append(col1)
            ascending_list.append(self.sort_dir1.get() == "Low to High")

        col2 = self.sort_col2.get()
        if col2 and col2 in self.nucleus_stats_df_original.columns:
            sort_criteria.append(col2)
            ascending_list.append(self.sort_dir2.get() == "Low to High")

        col3 = self.sort_col3.get()
        if col3 and col3 in self.nucleus_stats_df_original.columns:
            sort_criteria.append(col3)
            ascending_list.append(self.sort_dir3.get() == "Low to High")

        if not sort_criteria:
            messagebox.showwarning(
                "No Sort", "Please select at least one column to sort by"
            )
            return

        # Sort DataFrame
        self.nucleus_stats_df = self.nucleus_stats_df_original.sort_values(
            by=sort_criteria, ascending=ascending_list
        )

        # Repopulate table
        self.populate_table()

    def reset_sort(self):
        """Reset to original unsorted order."""
        if self.nucleus_stats_df_original is None:
            return

        self.nucleus_stats_df = self.nucleus_stats_df_original.copy()
        self.populate_table()

        # Clear dropdowns
        self.sort_col1.set("")
        self.sort_col2.set("")
        self.sort_col3.set("")

    def populate_table(self):
        """Populate table with nucleus statistics."""
        # Clear existing data
        self.clear_table()

        if self.nucleus_stats_df is None:
            return

        # Define columns (already in correct order from DataFrame)
        columns = list(self.nucleus_stats_df.columns)

        self.nucleus_tree["columns"] = columns
        self.nucleus_tree["show"] = "headings"

        # Configure column headers and widths
        for col in columns:
            self.nucleus_tree.heading(col, text=col)

            # Set column widths
            if col == "Nucleus_ID":
                width = 90
            elif col == "Matched":
                width = 70
            elif col == "Nucleus_Area":
                width = 100
            elif col == "Cyto_Area":
                width = 90
            elif "Mean" in col:
                width = 120
            elif "Max" in col:
                width = 110
            else:
                width = 100

            self.nucleus_tree.column(col, width=width, anchor=tk.CENTER)

        # Insert data
        for idx, row in self.nucleus_stats_df.iterrows():
            values = [row[col] for col in columns]
            # Format floats
            formatted_values = []
            for val in values:
                if isinstance(val, float):
                    formatted_values.append(f"{val:.1f}")
                else:
                    formatted_values.append(val)

            self.nucleus_tree.insert("", tk.END, values=formatted_values, tags=(idx,))

    def clear_table(self):
        """Clear table contents."""
        for item in self.nucleus_tree.get_children():
            self.nucleus_tree.delete(item)

    def on_nucleus_select(self, event):
        """Handle nucleus selection."""
        selection = self.nucleus_tree.selection()
        if selection:
            item = selection[0]
            values = self.nucleus_tree.item(item)["values"]
            nuc_id = values[0]
            self.selection_label.configure(
                text=f"Selected: Nucleus ID {nuc_id}", foreground="black"
            )
            self.plot_button.configure(state=tk.NORMAL)
        else:
            self.selection_label.configure(
                text="No nucleus selected", foreground="gray"
            )
            self.plot_button.configure(state=tk.DISABLED)

    def plot_nucleus(self):
        """Plot selected nucleus in popup window."""
        selection = self.nucleus_tree.selection()
        if not selection:
            messagebox.showwarning("Warning", "Please select a nucleus to plot")
            return

        item = selection[0]
        values = self.nucleus_tree.item(item)["values"]
        nuc_id = int(values[0])

        # Get nucleus data
        nuc_data = self.nucleus_data_dict[nuc_id]
        image_crop = nuc_data["image_crop"]
        nuc_mask = nuc_data["nucleus_mask_crop"]
        cell_mask = nuc_data["cell_mask_crop"]

        # Create popup window
        plot_window = tk.Toplevel(self.root)
        plot_window.title(f"Nucleus {nuc_id} - {Path(self.current_pkl_path).name}")

        # Create figure
        n_channels = image_crop.shape[0]
        fig_width = n_channels * 4
        fig_height = 5

        fig, axes = plt.subplots(1, n_channels, figsize=(fig_width, fig_height))

        if n_channels == 1:
            axes = [axes]

        # Plot each channel
        for ch in range(n_channels):
            ax = axes[ch]
            ch_img = image_crop[ch]

            # Calculate vmax (99th percentile)
            vmax = np.quantile(ch_img, 0.99)

            # Show image
            ax.imshow(ch_img, cmap="gray", vmin=0, vmax=vmax)

            # Overlay contours
            # Nucleus: firebrick
            ax.contour(
                nuc_mask.astype(float), levels=[0.5], colors="firebrick", linewidths=2
            )

            # Cell: orange (if exists)
            if np.sum(cell_mask) > 0:
                ax.contour(
                    cell_mask.astype(float), levels=[0.5], colors="orange", linewidths=2
                )

            ax.set_title(f"Channel {ch+1}", fontsize=10)
            ax.axis("off")

        # Add legend to the LAST (rightmost) subplot - bottom right, no frame
        legend_elements = [
            Line2D([0], [0], color="firebrick", linewidth=2, label="Nucleus"),
            Line2D([0], [0], color="orange", linewidth=2, label="Cell"),
        ]

        axes[-1].legend(
            handles=legend_elements,
            loc="lower right",
            frameon=False,
            fontsize=9,
            bbox_to_anchor=(1.0, -0.2),
        )

        plt.tight_layout()

        # Embed in tkinter
        canvas_frame = ttk.Frame(plot_window)
        canvas_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=(10, 5))

        canvas = FigureCanvasTkAgg(fig, master=canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Save button
        button_frame = ttk.Frame(plot_window)
        button_frame.pack(fill=tk.X, padx=10, pady=(0, 10))

        def save_plot():
            pkl_dir = Path(self.current_pkl_path).parent
            pkl_stem = Path(self.current_pkl_path).stem

            # Save PNG (plot with masks)
            png_filename = f"{pkl_stem}-nucleus{nuc_id}-allchwithmask.png"
            png_path = pkl_dir / png_filename
            fig.savefig(png_path, dpi=300, bbox_inches="tight")

            # Save TIFF (raw data: channels + masks)
            # Stack: [Ch1, Ch2, ..., ChN, Nucleus_Mask, Cell_Mask]
            tif_filename = f"{pkl_stem}-nucleus{nuc_id}.tif"
            tif_path = pkl_dir / tif_filename

            # Create stack with all channels and masks
            tif_stack = []

            # Add all image channels
            for ch in range(n_channels):
                tif_stack.append(image_crop[ch])

            # Add nucleus mask (convert bool to uint8 or uint16)
            tif_stack.append(nuc_mask.astype(np.uint16) * 65535)

            # Add cell mask (convert bool to uint8 or uint16)
            tif_stack.append(cell_mask.astype(np.uint16) * 65535)

            # Stack into single array (C, H, W)
            tif_array = np.stack(tif_stack, axis=0)

            # Save as TIFF
            tifffile.imwrite(tif_path, tif_array, imagej=True)

            messagebox.showinfo(
                "Saved",
                f"Files saved to:\n\n"
                f"Plot: {png_path.name}\n"
                f"TIFF: {tif_path.name}",
            )

        ttk.Button(button_frame, text="Save Plot", command=save_plot, width=15).pack(
            side=tk.RIGHT
        )
        ttk.Button(
            button_frame, text="Close", command=plot_window.destroy, width=15
        ).pack(side=tk.RIGHT, padx=5)

    def export_csv(self):
        """Export nucleus statistics to CSV."""
        if self.nucleus_stats_df is None:
            messagebox.showwarning("Warning", "No data to export")
            return

        pkl_dir = Path(self.current_pkl_path).parent
        pkl_stem = Path(self.current_pkl_path).stem
        filename = f"{pkl_stem}_nucleus_stats.csv"
        save_path = pkl_dir / filename

        self.nucleus_stats_df.to_csv(save_path, index=False)
        messagebox.showinfo("Exported", f"Table exported to:\n{save_path}")


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================


def main():
    root = tk.Tk()
    app = CellposeSAMViewer(root)
    root.mainloop()


if __name__ == "__main__":
    main()
