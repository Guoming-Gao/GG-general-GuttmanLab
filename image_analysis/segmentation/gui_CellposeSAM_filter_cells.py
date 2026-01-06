#!/usr/bin/env python3
"""
CellposeSAM_filter_cells_gui.py

GUI application for filtering CellposeSAM segmentation results across multiple files.
Shows 24 histograms (3 masks × 4 channels × 2 stats) in fixed grid layout.

Usage:
    python CellposeSAM_filter_cells_gui.py

Author: Generated for Cellpose-SAM result filtering
Date: 2025-11-08
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
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.lines import Line2D
import tifffile
from skimage.filters import threshold_otsu
from datetime import datetime
import time
import matplotlib

matplotlib.use("TkAgg")


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def calculate_channel_stats(image_channel, mask):
    """Calculate intensity statistics for a single channel within a mask."""
    if np.sum(mask) == 0:
        return {"mean": 0.0, "max": 0}

    pixels = image_channel[mask]
    return {"mean": float(np.mean(pixels)), "max": int(np.max(pixels))}


def calculate_otsu_threshold(data):
    """Calculate Otsu threshold with fallback for edge cases."""
    try:
        if len(data) == 0:
            return 0.0
        if len(set(data)) == 1 or np.min(data) == np.max(data):
            return float(data[0])
        return float(threshold_otsu(data))
    except:
        return float(np.median(data)) if len(data) > 0 else 0.0


def format_time(seconds):
    """Format seconds into readable time string."""
    if seconds < 60:
        return f"{int(seconds)}s"
    elif seconds < 3600:
        mins = int(seconds / 60)
        secs = int(seconds % 60)
        return f"{mins}m {secs}s"
    else:
        hours = int(seconds / 3600)
        mins = int((seconds % 3600) / 60)
        return f"{hours}h {mins}m"


# =============================================================================
# MAIN GUI CLASS
# =============================================================================


class CellposeSAMFilter:
    def __init__(self, root):
        self.root = root
        self.root.title("CellposeSAM Filter Cells")
        self.root.geometry("1600x1130")  # Slightly increased for progress section

        # Data storage
        self.selected_pkl_files = []
        self.all_nuclei_data = []
        self.filtered_nuclei_df = None
        self.pickle_cache = {}
        self.n_channels = 0
        self.metric_names = []

        # Histogram and threshold data
        self.otsu_thresholds = {}
        self.threshold_entries = {}
        self.threshold_modes = {}
        self.threshold_lines = {}
        self.histogram_figure = None
        self.histogram_canvas = None

        # Progress tracking
        self.analysis_queue = queue.Queue()
        self.export_queue = queue.Queue()
        self.analyzing = False
        self.exporting = False
        self.analysis_start_time = 0
        self.export_start_time = 0

        # Create GUI
        self.create_widgets()

        # Start queue polling
        self.poll_queues()

    def create_widgets(self):
        """Build all GUI components with fixed layout."""

        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # ===== FILE SELECTION =====
        file_header = ttk.Frame(main_frame)
        file_header.pack(fill=tk.X, pady=(0, 5))

        ttk.Label(
            file_header, text="Pickle Files (Multi-Select)", font=("Arial", 12, "bold")
        ).pack(side=tk.LEFT)
        ttk.Button(file_header, text="Select Files...", command=self.select_files).pack(
            side=tk.RIGHT, padx=5
        )
        ttk.Button(file_header, text="Clear", command=self.clear_files).pack(
            side=tk.RIGHT
        )
        ttk.Button(file_header, text="Analyze", command=self.analyze_files).pack(
            side=tk.RIGHT, padx=5
        )

        # File listbox
        list_frame = ttk.Frame(main_frame)
        list_frame.pack(fill=tk.X, pady=(0, 5))

        self.file_listbox = tk.Listbox(list_frame, height=3, selectmode=tk.MULTIPLE)
        list_scroll = ttk.Scrollbar(
            list_frame, orient=tk.VERTICAL, command=self.file_listbox.yview
        )
        self.file_listbox.configure(yscrollcommand=list_scroll.set)
        self.file_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        list_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        self.selected_file_label = ttk.Label(
            main_frame, text="No files selected", foreground="gray"
        )
        self.selected_file_label.pack(anchor=tk.W, pady=(0, 10))

        # ===== HISTOGRAM SECTION (FIXED 450px) =====
        self.histogram_frame = ttk.Frame(main_frame, height=450)
        self.histogram_frame.pack(fill=tk.X, pady=(0, 10))
        self.histogram_frame.pack_propagate(False)

        self.histogram_placeholder = ttk.Label(
            self.histogram_frame,
            text="Histograms will appear after clicking [Analyze]\n24 plots: 3 mask types × 4 channels × 2 statistics",
            font=("Arial", 10),
            foreground="gray",
        )
        self.histogram_placeholder.pack(expand=True)

        # ===== THRESHOLD CONTROLS SECTION (280px) =====
        self.threshold_frame = ttk.Frame(main_frame, height=280)
        self.threshold_frame.pack(fill=tk.X, pady=(0, 10))
        self.threshold_frame.pack_propagate(False)

        # ===== FILTER BUTTON =====
        self.filter_button_frame = ttk.Frame(main_frame)
        self.filter_button_frame.pack(fill=tk.X, pady=(0, 10))

        self.apply_filter_button = ttk.Button(
            self.filter_button_frame,
            text="Apply Filter",
            command=self.apply_filter,
            width=15,
            state=tk.DISABLED,
        )
        self.apply_filter_button.pack(side=tk.RIGHT)

        # ===== TABLE SECTION =====
        table_header = ttk.Frame(main_frame)
        table_header.pack(fill=tk.X, pady=(0, 5))

        self.table_label = ttk.Label(
            table_header, text="Filtered Cells: 0 / 0", font=("Arial", 11, "bold")
        )
        self.table_label.pack(side=tk.LEFT)

        ttk.Button(
            table_header, text="Export All as TIF", command=self.export_all_tif
        ).pack(side=tk.RIGHT, padx=5)
        ttk.Button(table_header, text="Export CSV", command=self.export_csv).pack(
            side=tk.RIGHT
        )

        # Table
        table_frame = ttk.Frame(main_frame)
        table_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 5))

        self.tree_scroll_y = ttk.Scrollbar(table_frame, orient=tk.VERTICAL)
        self.tree_scroll_x = ttk.Scrollbar(table_frame, orient=tk.HORIZONTAL)

        self.nucleus_tree = ttk.Treeview(
            table_frame,
            yscrollcommand=self.tree_scroll_y.set,
            xscrollcommand=self.tree_scroll_x.set,
            selectmode="browse",
            height=3,
        )

        self.tree_scroll_y.config(command=self.nucleus_tree.yview)
        self.tree_scroll_x.config(command=self.nucleus_tree.xview)

        self.nucleus_tree.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.tree_scroll_y.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.tree_scroll_x.grid(row=1, column=0, sticky=(tk.W, tk.E))

        table_frame.rowconfigure(0, weight=1)
        table_frame.columnconfigure(0, weight=1)

        self.nucleus_tree.bind("<<TreeviewSelect>>", self.on_nucleus_select)

        # Sort and plot controls
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(fill=tk.X, pady=(5, 10))

        ttk.Label(control_frame, text="Sort by:").pack(side=tk.LEFT, padx=(0, 5))

        self.sort_col1 = ttk.Combobox(control_frame, width=18, state="readonly")
        self.sort_col1.pack(side=tk.LEFT, padx=2)

        self.sort_dir1 = ttk.Combobox(
            control_frame,
            values=["High to Low", "Low to High"],
            width=12,
            state="readonly",
        )
        self.sort_dir1.set("High to Low")
        self.sort_dir1.pack(side=tk.LEFT, padx=2)

        ttk.Button(control_frame, text="Sort", command=self.apply_sort, width=8).pack(
            side=tk.LEFT, padx=10
        )

        self.selection_label = ttk.Label(
            control_frame, text="No nucleus selected", foreground="gray"
        )
        self.selection_label.pack(side=tk.LEFT, padx=20)

        self.plot_button = ttk.Button(
            control_frame,
            text="Plot Nucleus",
            command=self.plot_nucleus,
            state=tk.DISABLED,
            width=12,
        )
        self.plot_button.pack(side=tk.RIGHT)

        # ===== DEDICATED PROGRESS SECTION (Always at bottom) =====
        separator = ttk.Separator(main_frame, orient="horizontal")
        separator.pack(fill=tk.X, pady=(5, 5))

        progress_container = ttk.Frame(main_frame)
        progress_container.pack(fill=tk.X, side=tk.BOTTOM)

        # Progress label and bar
        self.unified_progress_label = ttk.Label(
            progress_container, text="Status: Ready", font=("Arial", 9)
        )
        self.unified_progress_label.pack(side=tk.LEFT, padx=(0, 10))

        self.unified_progress_bar = ttk.Progressbar(
            progress_container, mode="determinate", length=400
        )
        self.unified_progress_bar.pack(side=tk.LEFT, padx=(0, 10))

        self.unified_eta_label = ttk.Label(
            progress_container, text="", font=("Arial", 9), width=15
        )
        self.unified_eta_label.pack(side=tk.LEFT)

    def select_files(self):
        """Open file selection dialog."""
        files = filedialog.askopenfilenames(
            title="Select pickle files",
            filetypes=[("Pickle files", "*.pkl"), ("All files", "*.*")],
        )

        if files:
            self.selected_pkl_files = list(files)
            self.file_listbox.delete(0, tk.END)
            for f in self.selected_pkl_files:
                self.file_listbox.insert(tk.END, Path(f).name)

            self.selected_file_label.configure(
                text=f"Selected: {len(self.selected_pkl_files)} files",
                foreground="black",
            )

    def clear_files(self):
        """Clear file list and reset."""
        self.selected_pkl_files = []
        self.file_listbox.delete(0, tk.END)
        self.selected_file_label.configure(text="No files selected", foreground="gray")
        self.all_nuclei_data = []
        self.filtered_nuclei_df = None
        self.pickle_cache = {}

        if self.histogram_canvas:
            self.histogram_canvas.get_tk_widget().destroy()
            self.histogram_canvas = None

        self.histogram_placeholder.pack(expand=True)

        for widget in self.threshold_frame.winfo_children():
            widget.destroy()

        self.clear_table()
        self.table_label.configure(text="Filtered Cells: 0 / 0")

        self.apply_filter_button.configure(state=tk.DISABLED)
        self.plot_button.configure(state=tk.DISABLED)

        # Reset progress
        self.unified_progress_label.configure(text="Status: Ready")
        self.unified_progress_bar["value"] = 0
        self.unified_eta_label.configure(text="")

    def analyze_files(self):
        """Analyze all selected pickle files."""
        if not self.selected_pkl_files:
            messagebox.showerror("Error", "Please select pickle files to analyze")
            return

        self.analyzing = True
        self.analysis_start_time = time.time()
        thread = threading.Thread(target=self.analyze_worker, daemon=True)
        thread.start()

        # Update unified progress bar
        self.unified_progress_label.configure(text="Analysis: Starting...")
        self.unified_progress_bar["mode"] = "determinate"
        self.unified_progress_bar["maximum"] = len(self.selected_pkl_files)
        self.unified_progress_bar["value"] = 0
        self.unified_eta_label.configure(text="")

    def analyze_worker(self):
        """Worker thread for analysis."""
        try:
            all_nuclei = []
            n_channels = None

            for file_idx, pkl_path in enumerate(self.selected_pkl_files):
                if file_idx > 0:
                    elapsed = time.time() - self.analysis_start_time
                    avg_time = elapsed / file_idx
                    remaining = (len(self.selected_pkl_files) - file_idx) * avg_time
                    eta = f"ETA: {format_time(remaining)}"
                else:
                    eta = ""

                self.analysis_queue.put(
                    {
                        "type": "progress",
                        "current": file_idx,
                        "total": len(self.selected_pkl_files),
                        "message": f"Analysis: File {file_idx + 1}/{len(self.selected_pkl_files)}",
                        "eta": eta,
                    }
                )

                with open(pkl_path, "rb") as f:
                    pkl_data = pickle.load(f)

                self.pickle_cache[pkl_path] = pkl_data
                nucleus_data = pkl_data["nucleus_data"]

                if n_channels is None:
                    first_nuc = list(nucleus_data.values())[0]
                    n_channels = first_nuc["image_crop"].shape[0]

                for nuc_id, nuc_data in nucleus_data.items():
                    image_crop = nuc_data["image_crop"]
                    nuc_mask = nuc_data["nucleus_mask_crop"]
                    cell_mask = nuc_data["cell_mask_crop"]
                    cyto_mask = cell_mask & ~nuc_mask

                    nucleus_info = {
                        "source_file": pkl_path,
                        "source_filename": Path(pkl_path).name,
                        "nucleus_id": nuc_id,
                        "matched": nuc_data["is_matched"],
                        "nucleus_area": nuc_data["nucleus_area"],
                        "cyto_area": int(np.sum(cyto_mask)),
                    }

                    for ch in range(n_channels):
                        ch_num = ch + 1
                        ch_img = image_crop[ch]

                        cell_stats = calculate_channel_stats(ch_img, cell_mask)
                        nucleus_info[f"Cell_Ch{ch_num}_Mean"] = cell_stats["mean"]
                        nucleus_info[f"Cell_Ch{ch_num}_Max"] = cell_stats["max"]

                        nuc_stats = calculate_channel_stats(ch_img, nuc_mask)
                        nucleus_info[f"Nuc_Ch{ch_num}_Mean"] = nuc_stats["mean"]
                        nucleus_info[f"Nuc_Ch{ch_num}_Max"] = nuc_stats["max"]

                        cyto_stats = calculate_channel_stats(ch_img, cyto_mask)
                        nucleus_info[f"Cyto_Ch{ch_num}_Mean"] = cyto_stats["mean"]
                        nucleus_info[f"Cyto_Ch{ch_num}_Max"] = cyto_stats["max"]

                    all_nuclei.append(nucleus_info)

            self.analysis_queue.put(
                {"type": "complete", "all_nuclei": all_nuclei, "n_channels": n_channels}
            )

        except Exception as e:
            import traceback

            self.analysis_queue.put(
                {"type": "error", "error": str(e), "traceback": traceback.format_exc()}
            )

    def poll_queues(self):
        """Poll analysis and export queues."""
        try:
            while True:
                item = self.analysis_queue.get_nowait()

                if item["type"] == "progress":
                    self.unified_progress_bar["value"] = item["current"]
                    self.unified_progress_label.configure(text=item["message"])
                    self.unified_eta_label.configure(text=item.get("eta", ""))

                elif item["type"] == "complete":
                    self.all_nuclei_data = item["all_nuclei"]
                    self.n_channels = item["n_channels"]
                    self.analyzing = False

                    # Reset progress bar
                    self.unified_progress_label.configure(
                        text="Status: Analysis complete"
                    )
                    self.unified_progress_bar["value"] = 0
                    self.unified_eta_label.configure(text="")

                    self.create_histograms()
                    self.create_threshold_controls()

                    self.apply_filter_button.configure(state=tk.NORMAL)
                    self.table_label.configure(
                        text=f"Filtered Cells: 0 / {len(self.all_nuclei_data)} (Click Apply Filter)"
                    )

                elif item["type"] == "error":
                    messagebox.showerror("Error", f"Analysis failed:\n{item['error']}")
                    self.unified_progress_label.configure(
                        text="Status: Analysis failed"
                    )
                    self.unified_progress_bar["value"] = 0
                    self.unified_eta_label.configure(text="")
                    self.analyzing = False
        except queue.Empty:
            pass

        try:
            while True:
                item = self.export_queue.get_nowait()

                if item["type"] == "progress":
                    self.unified_progress_bar["value"] = item["current"]
                    self.unified_progress_label.configure(text=item["message"])
                    self.unified_eta_label.configure(text=item.get("eta", ""))

                elif item["type"] == "complete":
                    self.exporting = False
                    self.unified_progress_label.configure(
                        text="Status: Export complete"
                    )
                    self.unified_progress_bar["value"] = 0
                    self.unified_eta_label.configure(text="")
                    messagebox.showinfo("Success", item["message"])

                elif item["type"] == "error":
                    messagebox.showerror("Error", f"Export failed:\n{item['error']}")
                    self.unified_progress_label.configure(text="Status: Export failed")
                    self.unified_progress_bar["value"] = 0
                    self.unified_eta_label.configure(text="")
                    self.exporting = False
        except queue.Empty:
            pass

        self.root.after(100, self.poll_queues)

    def create_histograms(self):
        """Create 3x8 histogram grid."""
        self.histogram_placeholder.pack_forget()

        self.metric_names = []
        mask_types = ["Cell", "Nuc", "Cyto"]
        for ch in range(1, self.n_channels + 1):
            for mask in mask_types:
                for stat in ["Mean", "Max"]:
                    self.metric_names.append(f"{mask}_Ch{ch}_{stat}")

        fig, axes = plt.subplots(3, 8, figsize=(15.5, 4.5))
        fig.subplots_adjust(
            left=0.05, right=0.98, top=0.95, bottom=0.08, hspace=0.4, wspace=0.3
        )

        self.histogram_figure = fig
        self.threshold_lines = {}
        self.otsu_thresholds = {}

        fig.text(
            0.02, 0.83, "Cell", rotation=90, va="center", fontsize=11, weight="bold"
        )
        fig.text(
            0.02, 0.50, "Nucleus", rotation=90, va="center", fontsize=11, weight="bold"
        )
        fig.text(
            0.02,
            0.17,
            "Cytoplasm",
            rotation=90,
            va="center",
            fontsize=11,
            weight="bold",
        )

        stat_types = ["Mean", "Max"]

        for row, mask in enumerate(mask_types):
            for stat_idx, stat in enumerate(stat_types):
                for ch in range(1, self.n_channels + 1):
                    col = (ch - 1) + (stat_idx * self.n_channels)
                    ax = axes[row, col]

                    metric = f"{mask}_Ch{ch}_{stat}"

                    data = [
                        nuc[metric]
                        for nuc in self.all_nuclei_data
                        if nuc.get(metric, 0) > 0
                    ]

                    if len(data) == 0:
                        ax.text(
                            0.5,
                            0.5,
                            "No Data",
                            ha="center",
                            va="center",
                            fontsize=8,
                            transform=ax.transAxes,
                        )
                        ax.set_xticks([])
                        ax.set_yticks([])
                        self.otsu_thresholds[metric] = 0.0
                    else:
                        otsu_val = calculate_otsu_threshold(data)
                        self.otsu_thresholds[metric] = otsu_val

                        x_min = np.quantile(data, 0.005)
                        x_max = np.quantile(data, 0.995)

                        ax.hist(
                            data,
                            bins=30,
                            range=(x_min, x_max),
                            color="skyblue",
                            edgecolor="black",
                            alpha=0.7,
                        )

                        ax.axvline(
                            otsu_val,
                            color="blue",
                            linestyle="--",
                            linewidth=1,
                            label=f"Otsu",
                        )

                        user_line = ax.axvline(
                            otsu_val,
                            color="red",
                            linestyle="-",
                            linewidth=1,
                            label=f"User",
                        )
                        self.threshold_lines[metric] = user_line

                        ax.legend(fontsize=5, loc="upper right")
                        ax.tick_params(labelsize=5)

                    if row == 0:
                        title = f"Ch{ch} {stat}"
                        ax.set_title(title, fontsize=7)

                    if col == 0:
                        ax.set_ylabel("Count", fontsize=6)
                    if row == 2:
                        ax.set_xlabel("Intensity", fontsize=6)

        canvas = FigureCanvasTkAgg(fig, master=self.histogram_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.histogram_frame)
        toolbar.update()

        self.histogram_canvas = canvas

    def create_threshold_controls(self):
        """Create 6 rows × 4 columns threshold control grid."""
        for widget in self.threshold_frame.winfo_children():
            widget.destroy()

        ttk.Label(
            self.threshold_frame, text="Threshold Controls:", font=("Arial", 10, "bold")
        ).pack(anchor=tk.W, pady=(5, 5))

        grid_frame = ttk.Frame(self.threshold_frame)
        grid_frame.pack(fill=tk.BOTH, expand=True)

        headers = [
            "Cell Mean",
            "Cell Max",
            "Nuc Mean",
            "Nuc Max",
            "Cyto Mean",
            "Cyto Max",
        ]

        for col_idx in range(6):
            if col_idx < len(headers):
                ttk.Label(
                    grid_frame, text=headers[col_idx], font=("Arial", 9, "bold")
                ).grid(row=0, column=col_idx, padx=8, pady=(0, 5), sticky=tk.W)

        self.threshold_entries = {}
        self.threshold_modes = {}

        row_idx = 1
        for ch in range(1, self.n_channels + 1):
            self._create_control(grid_frame, row_idx, 0, f"Cell_Ch{ch}_Mean", f"Ch{ch}")
            self._create_control(grid_frame, row_idx, 1, f"Cell_Ch{ch}_Max", f"Ch{ch}")
            self._create_control(grid_frame, row_idx, 2, f"Nuc_Ch{ch}_Mean", f"Ch{ch}")
            self._create_control(grid_frame, row_idx, 3, f"Nuc_Ch{ch}_Max", f"Ch{ch}")
            self._create_control(grid_frame, row_idx, 4, f"Cyto_Ch{ch}_Mean", f"Ch{ch}")
            self._create_control(grid_frame, row_idx, 5, f"Cyto_Ch{ch}_Max", f"Ch{ch}")

            row_idx += 1

    def _create_control(self, parent, row, col, metric, label):
        """Helper to create a single threshold control."""
        control_frame = ttk.Frame(parent)
        control_frame.grid(row=row, column=col, sticky=tk.W, padx=5, pady=1)

        ttk.Label(control_frame, text=f"{label}:", width=4).pack(side=tk.LEFT)

        entry = ttk.Entry(control_frame, width=6)
        entry.insert(0, f"{self.otsu_thresholds.get(metric, 0):.1f}")
        entry.pack(side=tk.LEFT, padx=2)
        self.threshold_entries[metric] = entry

        entry.bind("<KeyRelease>", lambda e, m=metric: self.update_threshold_line(m))

        mode = ttk.Combobox(
            control_frame,
            values=["Not Used", "Above", "Below"],
            width=8,
            state="readonly",
        )
        mode.set("Not Used")
        mode.pack(side=tk.LEFT)
        self.threshold_modes[metric] = mode

    def update_threshold_line(self, metric):
        """Update threshold line when user changes value."""
        try:
            new_val = float(self.threshold_entries[metric].get())
            line = self.threshold_lines.get(metric)
            if line:
                line.set_xdata([new_val, new_val])
                self.histogram_canvas.draw_idle()
        except ValueError:
            pass

    def apply_filter(self):
        """Apply threshold filters."""
        if not self.all_nuclei_data:
            return

        filtered = []

        for nucleus in self.all_nuclei_data:
            passes = True

            for metric in self.metric_names:
                try:
                    threshold = float(self.threshold_entries[metric].get())
                    mode = self.threshold_modes[metric].get()
                    value = nucleus.get(metric, 0)

                    if mode == "Above" and value <= threshold:
                        passes = False
                        break
                    elif mode == "Below" and value >= threshold:
                        passes = False
                        break
                except:
                    continue

            if passes:
                filtered.append(nucleus)

        self.filtered_nuclei_df = pd.DataFrame(filtered)
        self.populate_table()

        total = len(self.all_nuclei_data)
        filtered_count = len(filtered)
        self.table_label.configure(
            text=f"Filtered Cells: {filtered_count} / {total} total"
        )

    def populate_table(self):
        """Populate table with filtered results."""
        self.clear_table()

        if self.filtered_nuclei_df is None or len(self.filtered_nuclei_df) == 0:
            return

        columns = [
            "source_filename",
            "nucleus_id",
            "matched",
            "nucleus_area",
            "cyto_area",
        ] + self.metric_names
        columns = [c for c in columns if c in self.filtered_nuclei_df.columns]

        self.nucleus_tree["columns"] = columns
        self.nucleus_tree["show"] = "headings"

        for col in columns:
            display = col.replace("source_filename", "File").replace("_", " ").title()
            self.nucleus_tree.heading(col, text=display)

            if "file" in col.lower():
                width = 150
            elif col in ["nucleus_id", "matched"]:
                width = 60
            elif "area" in col.lower():
                width = 80
            else:
                width = 90

            self.nucleus_tree.column(col, width=width, anchor=tk.CENTER)

        for idx, row in self.filtered_nuclei_df.iterrows():
            values = []
            for col in columns:
                val = row[col]
                if col == "matched":
                    val = "✓" if val else "✗"
                elif isinstance(val, float):
                    val = f"{val:.1f}"
                values.append(val)

            self.nucleus_tree.insert("", tk.END, values=values, tags=(idx,))

        self.update_sort_dropdown()

    def clear_table(self):
        """Clear table contents."""
        for item in self.nucleus_tree.get_children():
            self.nucleus_tree.delete(item)

    def update_sort_dropdown(self):
        """Update sort dropdown."""
        if self.filtered_nuclei_df is not None:
            columns = list(self.filtered_nuclei_df.columns)
            self.sort_col1["values"] = [""] + columns

    def apply_sort(self):
        """Apply sorting."""
        if self.filtered_nuclei_df is None or len(self.filtered_nuclei_df) == 0:
            return

        col = self.sort_col1.get()
        if col and col in self.filtered_nuclei_df.columns:
            ascending = self.sort_dir1.get() == "Low to High"
            self.filtered_nuclei_df = self.filtered_nuclei_df.sort_values(
                by=col, ascending=ascending
            )
            self.populate_table()

    def on_nucleus_select(self, event):
        """Handle nucleus selection."""
        selection = self.nucleus_tree.selection()
        if selection:
            item = selection[0]
            values = self.nucleus_tree.item(item)["values"]
            filename = values[0]
            nuc_id = values[1]
            self.selection_label.configure(
                text=f"Selected: {filename}, Nucleus {nuc_id}", foreground="black"
            )
            self.plot_button.configure(state=tk.NORMAL)
        else:
            self.selection_label.configure(
                text="No nucleus selected", foreground="gray"
            )
            self.plot_button.configure(state=tk.DISABLED)

    def plot_nucleus(self):
        """Plot selected nucleus."""
        selection = self.nucleus_tree.selection()
        if not selection:
            return

        item = selection[0]
        values = self.nucleus_tree.item(item)["values"]
        filename = values[0]
        nuc_id = int(values[1])

        source_path = None
        for path in self.selected_pkl_files:
            if Path(path).name == filename:
                source_path = path
                break

        if not source_path:
            messagebox.showerror("Error", "Source file not found")
            return

        pkl_data = self.pickle_cache.get(source_path)
        if not pkl_data:
            messagebox.showerror("Error", "Pickle data not loaded")
            return

        nuc_data = pkl_data["nucleus_data"].get(nuc_id)
        if not nuc_data:
            messagebox.showerror("Error", "Nucleus data not found")
            return

        image_crop = nuc_data["image_crop"]
        nuc_mask = nuc_data["nucleus_mask_crop"]
        cell_mask = nuc_data["cell_mask_crop"]

        plot_window = tk.Toplevel(self.root)
        plot_window.title(f"Nucleus {nuc_id} - {filename}")

        n_channels = image_crop.shape[0]
        fig, axes = plt.subplots(1, n_channels, figsize=(n_channels * 4, 5))

        if n_channels == 1:
            axes = [axes]

        for ch in range(n_channels):
            ax = axes[ch]
            ch_img = image_crop[ch]

            vmax = np.quantile(ch_img, 0.99)
            ax.imshow(ch_img, cmap="gray", vmin=0, vmax=vmax)

            ax.contour(
                nuc_mask.astype(float), levels=[0.5], colors="firebrick", linewidths=2
            )

            if np.sum(cell_mask) > 0:
                ax.contour(
                    cell_mask.astype(float), levels=[0.5], colors="orange", linewidths=2
                )

            ax.set_title(f"Channel {ch+1}", fontsize=10)
            ax.axis("off")

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

        canvas_frame = ttk.Frame(plot_window)
        canvas_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=(10, 5))

        canvas = FigureCanvasTkAgg(fig, master=canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        button_frame = ttk.Frame(plot_window)
        button_frame.pack(fill=tk.X, padx=10, pady=(0, 10))

        ttk.Button(
            button_frame, text="Close", command=plot_window.destroy, width=15
        ).pack(side=tk.RIGHT)

    def export_csv(self):
        """Export filtered results to CSV."""
        if self.filtered_nuclei_df is None or len(self.filtered_nuclei_df) == 0:
            messagebox.showwarning("Warning", "No data to export")
            return

        parent_dir = Path(self.selected_pkl_files[0]).parent
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"filtered_cells_{timestamp}.csv"
        save_path = parent_dir / filename

        self.filtered_nuclei_df.to_csv(save_path, index=False)
        messagebox.showinfo("Exported", f"Table exported to:\n{save_path}")

    def export_all_tif(self):
        """Export all filtered nuclei as TIFF files."""
        if self.filtered_nuclei_df is None or len(self.filtered_nuclei_df) == 0:
            messagebox.showwarning("Warning", "No filtered cells to export")
            return

        n_cells = len(self.filtered_nuclei_df)
        response = messagebox.askyesno(
            "Confirm Export",
            f"Export {n_cells} nuclei as TIFF files?\nThis may take several minutes.",
        )

        if not response:
            return

        self.exporting = True
        self.export_start_time = time.time()
        thread = threading.Thread(target=self.export_tif_worker, daemon=True)
        thread.start()

        # Update unified progress bar
        self.unified_progress_label.configure(text="Export: Starting...")
        self.unified_progress_bar["mode"] = "determinate"
        self.unified_progress_bar["maximum"] = n_cells
        self.unified_progress_bar["value"] = 0
        self.unified_eta_label.configure(text="")

    def export_tif_worker(self):
        """Worker thread for TIF export."""
        try:
            parent_dir = Path(self.selected_pkl_files[0]).parent
            output_dir = parent_dir / "filtered_cells"
            output_dir.mkdir(exist_ok=True)

            for idx, row in self.filtered_nuclei_df.iterrows():
                if idx > 0:
                    elapsed = time.time() - self.export_start_time
                    avg_time = elapsed / (idx + 1)
                    remaining = (len(self.filtered_nuclei_df) - idx - 1) * avg_time
                    eta = f"ETA: {format_time(remaining)}"
                else:
                    eta = ""

                source_file = row["source_file"]
                nuc_id = row["nucleus_id"]

                self.export_queue.put(
                    {
                        "type": "progress",
                        "current": idx + 1,
                        "total": len(self.filtered_nuclei_df),
                        "message": f"Export: Nucleus {idx + 1}/{len(self.filtered_nuclei_df)}",
                        "eta": eta,
                    }
                )

                pkl_data = self.pickle_cache.get(source_file)
                nuc_data = pkl_data["nucleus_data"][nuc_id]

                image_crop = nuc_data["image_crop"]
                nuc_mask = nuc_data["nucleus_mask_crop"]
                cell_mask = nuc_data["cell_mask_crop"]

                tif_stack = []
                for ch in range(image_crop.shape[0]):
                    tif_stack.append(image_crop[ch])

                tif_stack.append(nuc_mask.astype(np.uint16) * 65535)
                tif_stack.append(cell_mask.astype(np.uint16) * 65535)

                tif_array = np.stack(tif_stack, axis=0)

                stem = Path(source_file).stem
                filename = f"{stem}-nucleus{nuc_id}.tif"
                save_path = output_dir / filename

                tifffile.imwrite(save_path, tif_array, imagej=True)

            self.export_queue.put(
                {
                    "type": "complete",
                    "message": f"Exported {len(self.filtered_nuclei_df)} nuclei to:\n{output_dir}",
                }
            )

        except Exception as e:
            import traceback

            self.export_queue.put(
                {"type": "error", "error": str(e), "traceback": traceback.format_exc()}
            )


def main():
    root = tk.Tk()
    app = CellposeSAMFilter(root)
    root.mainloop()


if __name__ == "__main__":
    main()
