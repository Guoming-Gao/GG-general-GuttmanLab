#!/usr/bin/env python3
"""
run_CellposeSAM_GUI.py

GUI application for batch Cellpose-SAM segmentation with optimized storage.
Creates compact pickle files by storing cropped bounding boxes per nucleus.

Usage:
    python run_CellposeSAM_GUI.py

Author: Generated for Cellpose-SAM batch processing
Date: 2025-11-04
"""

import numpy as np
from cellpose import models, core, io
from pathlib import Path
import pickle
from datetime import datetime
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import threading
import queue
import time


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def normalize_img(img_in):
    """Normalize image to [0, 1] range."""
    img_float = img_in.astype(np.float32)
    max_val = np.max(img_float)
    if max_val > 0:
        img_float /= max_val
    return img_float


def get_bbox_with_padding(mask, padding, fov_shape):
    """
    Get bounding box coordinates for a mask with padding.

    Parameters:
    -----------
    mask : 2D boolean array
        Binary mask for the object
    padding : int
        Number of pixels to pad around the mask
    fov_shape : tuple
        (height, width) of the full field of view

    Returns:
    --------
    bbox : tuple
        (y_min, y_max, x_min, x_max) clipped to FOV boundaries
    """
    coords = np.where(mask)
    if len(coords[0]) == 0:
        return None

    y_min = max(0, coords[0].min() - padding)
    y_max = min(fov_shape[0], coords[0].max() + padding + 1)
    x_min = max(0, coords[1].min() - padding)
    x_max = min(fov_shape[1], coords[1].max() + padding + 1)

    return (y_min, y_max, x_min, x_max)


def calculate_padding(mask):
    """
    Calculate proportional padding based on object size.
    Uses 10% of object size with minimum of 7 pixels.

    Parameters:
    -----------
    mask : 2D boolean array
        Binary mask for the object

    Returns:
    --------
    padding : int
        Number of pixels to use as padding
    """
    coords = np.where(mask)
    if len(coords[0]) == 0:
        return 7

    height = coords[0].max() - coords[0].min()
    width = coords[1].max() - coords[1].min()

    # 10% of the larger dimension, minimum 7 pixels
    padding = max(7, int(0.1 * max(height, width)))

    return padding


def create_nucleus_data_dict(
    masks_cells, masks_nuclei, original_image, has_cyto_channel
):
    """
    Create optimized data structure with cropped bounding boxes per nucleus.

    Strategy:
    - For MATCHED nuclei: Crop based on cell boundary + padding
    - For UNMATCHED nuclei: Crop based on nucleus boundary + padding
    - Store cropped masks and original image for each nucleus

    Parameters:
    -----------
    masks_cells : 2D array (H, W)
        Cell segmentation masks
    masks_nuclei : 2D array (H, W)
        Nucleus segmentation masks
    original_image : array (C, H, W)
        Original multi-channel image
    has_cyto_channel : bool
        Whether cytoplasm channel was provided

    Returns:
    --------
    nucleus_data : dict
        Optimized per-nucleus data structure
    cell_to_nucleus_map : dict
        Mapping from cell_id to nucleus_id
    nucleus_to_cell_map : dict
        Mapping from nucleus_id to cell_id (or None)
    """
    nucleus_data = {}
    cell_to_nucleus_map = {}
    nucleus_to_cell_map = {}

    fov_shape = masks_cells.shape
    nucleus_ids = np.unique(masks_nuclei)[1:]

    if not has_cyto_channel:
        # Fast path: All nuclei are unmatched (no cell segmentation)
        # Skip the expensive overlap calculations
        for nucleus_id in nucleus_ids:
            nucleus_to_cell_map[int(nucleus_id)] = None

            nucleus_mask = masks_nuclei == nucleus_id

            # Use NUCLEUS mask for bounding box
            padding = calculate_padding(nucleus_mask)
            bbox = get_bbox_with_padding(nucleus_mask, padding, fov_shape)

            if bbox is None:
                continue

            y_min, y_max, x_min, x_max = bbox

            # Crop masks
            nucleus_mask_crop = nucleus_mask[y_min:y_max, x_min:x_max]
            cell_mask_crop = np.zeros_like(nucleus_mask_crop, dtype=bool)

            # Crop original image (all channels)
            image_crop = original_image[:, y_min:y_max, x_min:x_max]

            # Calculate statistics
            nucleus_area = int(np.sum(nucleus_mask))

            # Get centroid
            coords = np.where(nucleus_mask)
            centroid = (int(np.mean(coords[0])), int(np.mean(coords[1])))

            # Store data
            nucleus_data[int(nucleus_id)] = {
                "cell_id": None,
                "is_matched": False,
                "bbox": bbox,
                "bbox_size": (y_max - y_min, x_max - x_min),
                "padding_used": padding,
                "nucleus_mask_crop": nucleus_mask_crop,
                "cell_mask_crop": cell_mask_crop,
                "image_crop": image_crop,
                "nucleus_area": nucleus_area,
                "cell_area": 0,
                "centroid": centroid,
            }

        return nucleus_data, cell_to_nucleus_map, nucleus_to_cell_map

    # Normal path: Match nuclei to cells
    cell_ids = np.unique(masks_cells)[1:]
    matched_nuclei = set()

    # First pass: Match nuclei to cells
    for cell_id in cell_ids:
        cell_mask = masks_cells == cell_id

        # Find overlapping nuclei
        overlaps = {}
        for nucleus_id in nucleus_ids:
            if nucleus_id in matched_nuclei:
                continue
            nucleus_mask = masks_nuclei == nucleus_id
            overlap_area = np.sum(cell_mask & nucleus_mask)
            if overlap_area > 0:
                overlaps[nucleus_id] = overlap_area

        if overlaps:
            # Get nucleus with maximum overlap
            best_nucleus_id = max(overlaps, key=overlaps.get)
            matched_nuclei.add(best_nucleus_id)
            cell_to_nucleus_map[int(cell_id)] = int(best_nucleus_id)
            nucleus_to_cell_map[int(best_nucleus_id)] = int(cell_id)

    # Mark unmatched nuclei
    for nucleus_id in nucleus_ids:
        if nucleus_id not in matched_nuclei:
            nucleus_to_cell_map[int(nucleus_id)] = None

    # Second pass: Create cropped data for each nucleus
    for nucleus_id in nucleus_ids:
        nucleus_mask = masks_nuclei == nucleus_id
        cell_id = nucleus_to_cell_map[int(nucleus_id)]
        is_matched = cell_id is not None

        # Determine which mask to use for bounding box
        if is_matched:
            # Use CELL mask for padding calculation and bounding box
            cell_mask = masks_cells == cell_id
            reference_mask = cell_mask
            padding = calculate_padding(cell_mask)
        else:
            # Use NUCLEUS mask for padding calculation and bounding box
            reference_mask = nucleus_mask
            padding = calculate_padding(nucleus_mask)

        # Get bounding box
        bbox = get_bbox_with_padding(reference_mask, padding, fov_shape)

        if bbox is None:
            continue

        y_min, y_max, x_min, x_max = bbox

        # Crop masks
        nucleus_mask_crop = nucleus_mask[y_min:y_max, x_min:x_max]

        if is_matched:
            cell_mask = masks_cells == cell_id
            cell_mask_crop = cell_mask[y_min:y_max, x_min:x_max]
        else:
            cell_mask_crop = np.zeros_like(nucleus_mask_crop, dtype=bool)

        # Crop original image (all channels)
        image_crop = original_image[:, y_min:y_max, x_min:x_max]

        # Calculate statistics
        nucleus_area = int(np.sum(nucleus_mask))
        cell_area = int(np.sum(masks_cells == cell_id)) if is_matched else 0

        # Get centroid
        coords = np.where(nucleus_mask)
        centroid = (int(np.mean(coords[0])), int(np.mean(coords[1])))

        # Store data
        nucleus_data[int(nucleus_id)] = {
            "cell_id": cell_id,
            "is_matched": is_matched,
            "bbox": bbox,
            "bbox_size": (y_max - y_min, x_max - x_min),
            "padding_used": padding,
            "nucleus_mask_crop": nucleus_mask_crop,
            "cell_mask_crop": cell_mask_crop,
            "image_crop": image_crop,
            "nucleus_area": nucleus_area,
            "cell_area": cell_area,
            "centroid": centroid,
        }

    return nucleus_data, cell_to_nucleus_map, nucleus_to_cell_map


# =============================================================================
# MAIN GUI CLASS
# =============================================================================


class CellposeSAMGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Cellpose-SAM Batch Segmentation")
        self.root.geometry("700x1050")  # Reduced width from 1000px to 700px

        # Data storage
        self.selected_files = []
        self.output_dir = None
        self.processing = False
        self.cancel_requested = False
        self.result_queue = queue.Queue()
        self.log_messages = []

        # Create GUI
        self.create_widgets()

        # Start queue polling
        self.poll_queue()

    def create_widgets(self):
        """Build all GUI components."""

        # Main container with padding
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        row = 0

        # ===== FILE SELECTION =====
        ttk.Label(main_frame, text="Input Files", font=("Arial", 12, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        ttk.Button(main_frame, text="Select Files...", command=self.select_files).grid(
            row=row, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 5)
        )
        row += 1

        # File list
        list_frame = ttk.Frame(main_frame)
        list_frame.grid(
            row=row,
            column=0,
            columnspan=3,
            sticky=(tk.W, tk.E, tk.N, tk.S),
            pady=(0, 10),
        )

        self.file_listbox = tk.Listbox(list_frame, height=6)
        scrollbar = ttk.Scrollbar(
            list_frame, orient=tk.VERTICAL, command=self.file_listbox.yview
        )
        self.file_listbox.configure(yscrollcommand=scrollbar.set)
        self.file_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        row += 1

        # ===== OUTPUT DIRECTORY =====
        ttk.Label(main_frame, text="Output Directory", font=("Arial", 12, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(10, 5)
        )
        row += 1

        # Add explanation label
        explanation_label = ttk.Label(
            main_frame,
            text="A 'cellposeSAM_segmentation' folder will be created in the selected directory",
            font=("Arial", 9, "italic"),
            foreground="gray",
        )
        explanation_label.grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        ttk.Button(
            main_frame, text="Select Output Folder...", command=self.select_output
        ).grid(row=row, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 5))
        row += 1

        self.output_label = ttk.Label(
            main_frame, text="No folder selected", foreground="gray"
        )
        self.output_label.grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 10)
        )
        row += 1

        # ===== CHANNEL CONFIGURATION =====
        ttk.Label(
            main_frame, text="Channel Configuration", font=("Arial", 12, "bold")
        ).grid(row=row, column=0, columnspan=3, sticky=tk.W, pady=(10, 5))
        row += 1

        ttk.Label(main_frame, text="Cytoplasm Channel:").grid(
            row=row, column=0, sticky=tk.W, pady=2
        )
        self.cyto_channel = ttk.Combobox(
            main_frame, values=["None", "1", "2", "3", "4"], state="readonly", width=10
        )
        self.cyto_channel.set("2")
        self.cyto_channel.grid(row=row, column=1, sticky=tk.W, pady=2, padx=(5, 10))
        row += 1

        ttk.Label(main_frame, text="Nucleus Channel:").grid(
            row=row, column=0, sticky=tk.W, pady=2
        )
        self.nuc_channel = ttk.Combobox(
            main_frame, values=["1", "2", "3", "4"], state="readonly", width=10
        )
        self.nuc_channel.set("3")
        self.nuc_channel.grid(
            row=row, column=1, sticky=tk.W, pady=(2, 10), padx=(5, 10)
        )
        row += 1

        # ===== SEGMENTATION PARAMETERS (3-COLUMN LAYOUT) =====
        ttk.Label(
            main_frame, text="Segmentation Parameters", font=("Arial", 12, "bold")
        ).grid(row=row, column=0, columnspan=3, sticky=tk.W, pady=(10, 5))
        row += 1

        # Cell Diameter
        ttk.Label(main_frame, text="Cell Diameter:").grid(
            row=row, column=0, sticky=tk.W, pady=2
        )
        self.cell_diameter = ttk.Entry(main_frame, width=10)
        self.cell_diameter.insert(0, "70")
        self.cell_diameter.grid(row=row, column=1, sticky=tk.W, pady=2, padx=(5, 10))
        ttk.Label(
            main_frame,
            text="Expected cell diameter in pixels. Measure typical cells in ImageJ to estimate.",
            font=("Arial", 9, "italic"),
            foreground="gray",
            wraplength=400,
            justify=tk.LEFT,
        ).grid(row=row, column=2, sticky=tk.W, pady=2)
        row += 1

        # Cell Flow Threshold
        ttk.Label(main_frame, text="Cell Flow Threshold:").grid(
            row=row, column=0, sticky=tk.W, pady=2
        )
        self.cell_flow = ttk.Entry(main_frame, width=10)
        self.cell_flow.insert(0, "0.5")
        self.cell_flow.grid(row=row, column=1, sticky=tk.W, pady=2, padx=(5, 10))
        ttk.Label(
            main_frame,
            text="Max flow error (0-1). Default 0.4-0.5 works well. Increase to 0.6-0.8 if too few masks detected, decrease to 0.2-0.3 if getting ill-shaped masks.",
            font=("Arial", 9, "italic"),
            foreground="gray",
            wraplength=400,
            justify=tk.LEFT,
        ).grid(row=row, column=2, sticky=tk.W, pady=2)
        row += 1

        # Nucleus Diameter
        ttk.Label(main_frame, text="Nucleus Diameter:").grid(
            row=row, column=0, sticky=tk.W, pady=2
        )
        self.nuc_diameter = ttk.Entry(main_frame, width=10)
        self.nuc_diameter.insert(0, "50")
        self.nuc_diameter.grid(row=row, column=1, sticky=tk.W, pady=2, padx=(5, 10))
        ttk.Label(
            main_frame,
            text="Expected nucleus diameter in pixels. Measure typical nuclei in ImageJ to estimate.",
            font=("Arial", 9, "italic"),
            foreground="gray",
            wraplength=400,
            justify=tk.LEFT,
        ).grid(row=row, column=2, sticky=tk.W, pady=2)
        row += 1

        # Nucleus Flow Threshold
        ttk.Label(main_frame, text="Nucleus Flow Threshold:").grid(
            row=row, column=0, sticky=tk.W, pady=2
        )
        self.nuc_flow = ttk.Entry(main_frame, width=10)
        self.nuc_flow.insert(0, "0.5")
        self.nuc_flow.grid(row=row, column=1, sticky=tk.W, pady=(2, 10), padx=(5, 10))
        ttk.Label(
            main_frame,
            text="Max flow error (0-1). Default 0.4-0.5 works well. Increase to 0.6-0.8 if too few masks detected, decrease to 0.2-0.3 if getting ill-shaped masks.",
            font=("Arial", 9, "italic"),
            foreground="gray",
            wraplength=400,
            justify=tk.LEFT,
        ).grid(row=row, column=2, sticky=tk.W, pady=(2, 10))
        row += 1

        # ===== CONTROL BUTTONS =====
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=row, column=0, columnspan=3, pady=(10, 10))

        self.start_button = ttk.Button(
            button_frame,
            text="Start Processing",
            command=self.start_processing,
            width=20,
        )
        self.start_button.pack(side=tk.LEFT, padx=5)

        self.cancel_button = ttk.Button(
            button_frame,
            text="Cancel",
            command=self.cancel_processing,
            state=tk.DISABLED,
            width=15,
        )
        self.cancel_button.pack(side=tk.LEFT, padx=5)
        row += 1

        # ===== PROGRESS =====
        ttk.Label(main_frame, text="Progress", font=("Arial", 12, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(10, 5)
        )
        row += 1

        self.progress_label = ttk.Label(main_frame, text="Ready to process")
        self.progress_label.grid(row=row, column=0, columnspan=3, sticky=tk.W, pady=2)
        row += 1

        self.progress_bar = ttk.Progressbar(main_frame, mode="determinate", length=400)
        self.progress_bar.grid(
            row=row, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=2
        )
        row += 1

        self.eta_label = ttk.Label(main_frame, text="")
        self.eta_label.grid(row=row, column=0, columnspan=3, sticky=tk.W, pady=(2, 10))
        row += 1

        # ===== RESULTS TABLE =====
        ttk.Label(main_frame, text="Results Summary", font=("Arial", 12, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(10, 5)
        )
        row += 1

        # Create treeview for results
        tree_frame = ttk.Frame(main_frame)
        tree_frame.grid(
            row=row,
            column=0,
            columnspan=3,
            sticky=(tk.W, tk.E, tk.N, tk.S),
            pady=(0, 10),
        )

        self.results_tree = ttk.Treeview(
            tree_frame,
            columns=("File", "Nuclei", "Matched", "Unmatched"),
            show="headings",
            height=6,
        )
        self.results_tree.heading("File", text="Image File")
        self.results_tree.heading("Nuclei", text="Nuclei")
        self.results_tree.heading("Matched", text="Matched")
        self.results_tree.heading("Unmatched", text="Unmatched")

        self.results_tree.column("File", width=300)
        self.results_tree.column("Nuclei", width=80, anchor=tk.CENTER)
        self.results_tree.column("Matched", width=80, anchor=tk.CENTER)
        self.results_tree.column("Unmatched", width=100, anchor=tk.CENTER)

        tree_scroll = ttk.Scrollbar(
            tree_frame, orient=tk.VERTICAL, command=self.results_tree.yview
        )
        self.results_tree.configure(yscrollcommand=tree_scroll.set)

        self.results_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        tree_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        row += 1

        # Total row
        self.total_label = ttk.Label(
            main_frame,
            text="Total: 0 nuclei, 0 matched, 0 unmatched",
            font=("Arial", 10, "bold"),
        )
        self.total_label.grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(5, 10)
        )
        row += 1

        # ===== LOG WINDOW =====
        ttk.Label(main_frame, text="Log", font=("Arial", 12, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(10, 5)
        )
        row += 1

        self.log_text = scrolledtext.ScrolledText(
            main_frame, height=8, width=70, wrap=tk.WORD, state=tk.DISABLED
        )
        self.log_text.grid(
            row=row, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S)
        )

        # Configure grid weights for resizing
        main_frame.columnconfigure(0, weight=0)  # Labels - fixed width
        main_frame.columnconfigure(1, weight=0)  # Input boxes - fixed width
        main_frame.columnconfigure(2, weight=1)  # Descriptions - expandable
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

    def log(self, message):
        """Add message to log window and list."""
        timestamp = datetime.now().strftime("%H:%M:%S")
        log_entry = f"[{timestamp}] {message}"
        self.log_messages.append(log_entry)

        self.log_text.configure(state=tk.NORMAL)
        self.log_text.insert(tk.END, log_entry + "\n")
        self.log_text.see(tk.END)
        self.log_text.configure(state=tk.DISABLED)

    def select_files(self):
        """Open file selection dialog."""
        files = filedialog.askopenfilenames(
            title="Select image files",
            filetypes=[("TIFF files", "*.tif *.tiff"), ("All files", "*.*")],
        )

        if files:
            self.selected_files = list(files)
            self.file_listbox.delete(0, tk.END)
            for f in self.selected_files:
                self.file_listbox.insert(tk.END, Path(f).name)
            self.log(f"Selected {len(self.selected_files)} files")

    def select_output(self):
        """Open output directory selection dialog."""
        folder = filedialog.askdirectory(title="Select output folder")

        if folder:
            self.output_dir = Path(folder) / "cellposeSAM_segmentation"
            self.output_label.configure(text=str(self.output_dir), foreground="black")
            self.log(f"Output directory: {self.output_dir}")

    def validate_inputs(self):
        """Validate all user inputs."""
        if not self.selected_files:
            messagebox.showerror("Error", "Please select input files")
            return False

        if self.output_dir is None:
            messagebox.showerror("Error", "Please select output folder")
            return False

        if self.nuc_channel.get() == "":
            messagebox.showerror("Error", "Nucleus channel must be selected")
            return False

        # Validate numeric parameters
        try:
            if self.cyto_channel.get() != "None":
                float(self.cell_diameter.get())
                float(self.cell_flow.get())
            float(self.nuc_diameter.get())
            float(self.nuc_flow.get())
        except ValueError:
            messagebox.showerror("Error", "All parameter values must be numeric")
            return False

        return True

    def start_processing(self):
        """Start processing in a separate thread."""
        if not self.validate_inputs():
            return

        # Clear previous results
        for item in self.results_tree.get_children():
            self.results_tree.delete(item)
        self.total_label.configure(text="Total: 0 nuclei, 0 matched, 0 unmatched")
        self.progress_bar["value"] = 0

        # Disable controls
        self.start_button.configure(state=tk.DISABLED)
        self.cancel_button.configure(state=tk.NORMAL)
        self.processing = True
        self.cancel_requested = False

        # Start worker thread
        thread = threading.Thread(target=self.process_files_worker, daemon=True)
        thread.start()

        self.log("=" * 60)
        self.log("Processing started")

    def cancel_processing(self):
        """Request cancellation of processing."""
        self.cancel_requested = True
        self.log("Cancellation requested...")
        self.cancel_button.configure(state=tk.DISABLED)

    def process_files_worker(self):
        """Main processing loop (runs in worker thread)."""
        try:
            # Create output directory
            self.output_dir.mkdir(exist_ok=True, parents=True)

            # Get parameters - Convert from 1-based to 0-based indexing
            cyto_ch_str = self.cyto_channel.get()
            cyto_ch = None if cyto_ch_str == "None" else int(cyto_ch_str) - 1
            nuc_ch = int(self.nuc_channel.get()) - 1

            # Get segmentation parameters
            if cyto_ch is not None:
                cell_diam = float(self.cell_diameter.get())
                cell_flow_th = float(self.cell_flow.get())
            nuc_diam = float(self.nuc_diameter.get())
            nuc_flow_th = float(self.nuc_flow.get())

            # Initialize model
            self.result_queue.put(
                {"type": "log", "message": "Loading Cellpose model..."}
            )
            model = models.CellposeModel(gpu=core.use_gpu())
            self.result_queue.put(
                {"type": "log", "message": "Model loaded successfully"}
            )

            # Log channel configuration
            if cyto_ch is not None:
                self.result_queue.put(
                    {
                        "type": "log",
                        "message": f"Channel config: Cyto=Ch{cyto_ch_str} (index {cyto_ch}), Nuc=Ch{self.nuc_channel.get()} (index {nuc_ch})",
                    }
                )
            else:
                self.result_queue.put(
                    {
                        "type": "log",
                        "message": f"Channel config: Cyto=None, Nuc=Ch{self.nuc_channel.get()} (index {nuc_ch})",
                    }
                )
                self.result_queue.put(
                    {
                        "type": "log",
                        "message": "No cytoplasm channel - cell parameters will be ignored",
                    }
                )

            # Process files
            total_files = len(self.selected_files)
            start_time = time.time()

            for i, file_path in enumerate(self.selected_files):
                if self.cancel_requested:
                    self.result_queue.put(
                        {"type": "log", "message": "Processing cancelled by user"}
                    )
                    break

                file_path = Path(file_path)

                # Update progress at START
                self.result_queue.put(
                    {
                        "type": "progress_start",
                        "current": i + 1,
                        "total": total_files,
                        "file": file_path.name,
                    }
                )

                try:
                    # Load image
                    img = io.imread(str(file_path))
                    self.result_queue.put(
                        {
                            "type": "log",
                            "message": f"Loaded {file_path.name}: shape={img.shape}",
                        }
                    )

                    # Segment based on whether cytoplasm channel is provided
                    if cyto_ch is not None:
                        # === CASE 1: Both cytoplasm and nucleus channels ===

                        # Segment cells
                        self.result_queue.put(
                            {"type": "log", "message": "Starting cell segmentation..."}
                        )
                        seg_start = time.time()

                        selected_channels = [cyto_ch, nuc_ch]
                        img_selected = np.zeros_like(img)
                        img_selected[:2, :, :] = normalize_img(
                            img[selected_channels, :, :]
                        )

                        masks_cells, _, _ = model.eval(
                            img_selected,
                            batch_size=32,
                            diameter=cell_diam,
                            flow_threshold=cell_flow_th,
                            cellprob_threshold=0.0,
                            normalize={"tile_norm_blocksize": 0},
                        )

                        seg_elapsed = time.time() - seg_start
                        self.result_queue.put(
                            {
                                "type": "log",
                                "message": f"Cell segmentation completed in {seg_elapsed:.1f} sec",
                            }
                        )

                        # Segment nuclei
                        self.result_queue.put(
                            {
                                "type": "log",
                                "message": "Starting nucleus segmentation...",
                            }
                        )
                        seg_start = time.time()

                        img_nuc_only = np.zeros_like(img)
                        img_nuc_only[0, :, :] = normalize_img(img[nuc_ch, :, :])

                        masks_nuclei, _, _ = model.eval(
                            img_nuc_only,
                            batch_size=32,
                            diameter=nuc_diam,
                            flow_threshold=nuc_flow_th,
                            cellprob_threshold=0.0,
                            normalize={"tile_norm_blocksize": 0},
                        )

                        seg_elapsed = time.time() - seg_start
                        self.result_queue.put(
                            {
                                "type": "log",
                                "message": f"Nucleus segmentation completed in {seg_elapsed:.1f} sec",
                            }
                        )

                        has_cyto = True

                    else:
                        # === CASE 2: Only nucleus channel (no cytoplasm) ===

                        # Segment nuclei only
                        self.result_queue.put(
                            {
                                "type": "log",
                                "message": "Starting nucleus segmentation...",
                            }
                        )
                        seg_start = time.time()

                        img_nuc_only = np.zeros_like(img)
                        img_nuc_only[0, :, :] = normalize_img(img[nuc_ch, :, :])

                        masks_nuclei, _, _ = model.eval(
                            img_nuc_only,
                            batch_size=32,
                            diameter=nuc_diam,
                            flow_threshold=nuc_flow_th,
                            cellprob_threshold=0.0,
                            normalize={"tile_norm_blocksize": 0},
                        )

                        seg_elapsed = time.time() - seg_start
                        self.result_queue.put(
                            {
                                "type": "log",
                                "message": f"Nucleus segmentation completed in {seg_elapsed:.1f} sec",
                            }
                        )

                        # Use nucleus masks as cell masks
                        masks_cells = masks_nuclei.copy()

                        has_cyto = False

                    # Process nucleus data
                    self.result_queue.put(
                        {"type": "log", "message": "Processing nucleus data..."}
                    )
                    process_start = time.time()

                    nucleus_data, cell_to_nucleus_map, nucleus_to_cell_map = (
                        create_nucleus_data_dict(
                            masks_cells, masks_nuclei, img, has_cyto
                        )
                    )

                    process_elapsed = time.time() - process_start
                    self.result_queue.put(
                        {
                            "type": "log",
                            "message": f"Data processing completed in {process_elapsed:.1f} sec",
                        }
                    )

                    # Save results
                    self.result_queue.put(
                        {"type": "log", "message": "Saving results..."}
                    )
                    save_start = time.time()

                    # Save masks as TIFF
                    if has_cyto:
                        # Save both cell and nucleus masks
                        io.imsave(
                            str(
                                self.output_dir / (file_path.stem + "_masks_cells.tif")
                            ),
                            masks_cells,
                        )
                        io.imsave(
                            str(
                                self.output_dir / (file_path.stem + "_masks_nuclei.tif")
                            ),
                            masks_nuclei,
                        )
                    else:
                        # Only save nucleus masks (no duplicate cell mask)
                        io.imsave(
                            str(
                                self.output_dir / (file_path.stem + "_masks_nuclei.tif")
                            ),
                            masks_nuclei,
                        )
                        self.result_queue.put(
                            {
                                "type": "log",
                                "message": "Skipping cell mask save (no cytoplasm channel)",
                            }
                        )

                    # Calculate statistics
                    n_nuclei_total = len(nucleus_data)
                    n_nuclei_matched = sum(
                        1 for data in nucleus_data.values() if data["is_matched"]
                    )
                    n_nuclei_unmatched = n_nuclei_total - n_nuclei_matched
                    n_cells = len(cell_to_nucleus_map)

                    # Save optimized pickle
                    pickle_data = {
                        # Metadata
                        "image_name": file_path.name,
                        "image_stem": file_path.stem,
                        "processing_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        "image_shape": img.shape,
                        # Parameters
                        "parameters": {
                            "cyto_channel": cyto_ch_str,
                            "nuc_channel": self.nuc_channel.get(),
                            "cyto_channel_0based": cyto_ch,
                            "nuc_channel_0based": nuc_ch,
                            "cell_diameter": cell_diam if cyto_ch is not None else None,
                            "cell_flow_threshold": (
                                cell_flow_th if cyto_ch is not None else None
                            ),
                            "nucleus_diameter": nuc_diam,
                            "nucleus_flow_threshold": nuc_flow_th,
                            "has_cytoplasm_channel": cyto_ch is not None,
                        },
                        # Full FOV masks (small, for visualization)
                        "masks_cells": masks_cells,
                        "masks_nuclei": masks_nuclei,
                        # Optimized per-nucleus cropped data (MAIN DATA)
                        "nucleus_data": nucleus_data,
                        # Quick lookup maps
                        "cell_to_nucleus_map": cell_to_nucleus_map,
                        "nucleus_to_cell_map": nucleus_to_cell_map,
                        # Summary statistics
                        "summary": {
                            "n_nuclei_total": n_nuclei_total,
                            "n_nuclei_matched": n_nuclei_matched,
                            "n_nuclei_unmatched": n_nuclei_unmatched,
                            "n_cells": n_cells,
                            "has_cytoplasm_channel": cyto_ch is not None,
                        },
                        # Documentation
                        "readme": {
                            "description": "Optimized Cellpose-SAM results with cropped bounding boxes per nucleus",
                            "storage_strategy": (
                                "Each nucleus stored with: "
                                "cropped masks (nucleus + cell), "
                                "cropped original image (all channels), "
                                "bounding box coordinates. "
                                "Padding based on 10% of object size (min 7px). "
                                "Matched nuclei use cell boundary, unmatched use nucleus boundary."
                            ),
                            "note_no_cyto": (
                                "If cyto_channel='None', only nucleus segmentation was performed. "
                                "masks_cells and masks_nuclei are identical in pickle. "
                                "Only nucleus mask TIFF was saved (no duplicate cell mask). "
                                "All nuclei are marked as 'unmatched'."
                            ),
                        },
                    }

                    with open(
                        self.output_dir / (file_path.stem + "_segmentation.pkl"), "wb"
                    ) as f:
                        pickle.dump(pickle_data, f)

                    save_elapsed = time.time() - save_start
                    self.result_queue.put(
                        {
                            "type": "log",
                            "message": f"Save completed in {save_elapsed:.1f} sec",
                        }
                    )

                    # Send result
                    self.result_queue.put(
                        {
                            "type": "result",
                            "file": file_path.name,
                            "n_nuclei": n_nuclei_total,
                            "n_matched": n_nuclei_matched,
                            "n_unmatched": n_nuclei_unmatched,
                        }
                    )

                    # Update progress bar AFTER completion
                    percent = ((i + 1) / total_files) * 100
                    self.result_queue.put(
                        {"type": "progress_complete", "percent": percent}
                    )

                    # Calculate ETA
                    elapsed = time.time() - start_time
                    avg_time = elapsed / (i + 1)
                    remaining = (total_files - i - 1) * avg_time
                    eta_min = int(remaining // 60)
                    eta_sec = int(remaining % 60)

                    self.result_queue.put(
                        {
                            "type": "eta",
                            "eta": (
                                f"ETA: {eta_min} min {eta_sec} sec"
                                if remaining > 0
                                else "Complete"
                            ),
                        }
                    )

                except Exception as e:
                    import traceback

                    error_msg = f"{str(e)}\n{traceback.format_exc()}"
                    self.result_queue.put(
                        {"type": "error", "file": file_path.name, "error": error_msg}
                    )

            # Save log file
            log_path = self.output_dir / "processing_log.txt"
            with open(log_path, "w") as f:
                f.write("\n".join(self.log_messages))

            self.result_queue.put({"type": "done"})

        except Exception as e:
            import traceback

            error_msg = f"Fatal error: {str(e)}\n{traceback.format_exc()}"
            self.result_queue.put({"type": "error", "error": error_msg})
            self.result_queue.put({"type": "done"})

    def poll_queue(self):
        """Poll the result queue and update GUI."""
        try:
            while True:
                item = self.result_queue.get_nowait()

                if item["type"] == "log":
                    self.log(item["message"])

                elif item["type"] == "progress_start":
                    self.progress_label.configure(
                        text=f"Processing: {item['file']} ({item['current']}/{item['total']})"
                    )

                elif item["type"] == "progress_complete":
                    self.progress_bar["value"] = item["percent"]

                elif item["type"] == "eta":
                    self.eta_label.configure(text=item["eta"])

                elif item["type"] == "result":
                    self.results_tree.insert(
                        "",
                        tk.END,
                        values=(
                            item["file"],
                            item["n_nuclei"],
                            item["n_matched"],
                            item["n_unmatched"],
                        ),
                    )
                    self.results_tree.yview_moveto(1)
                    self.update_totals()
                    self.log(
                        f"Completed: {item['file']} - {item['n_nuclei']} nuclei, {item['n_matched']} matched"
                    )

                elif item["type"] == "error":
                    self.log(f"ERROR: {item.get('file', 'Unknown')}: {item['error']}")

                elif item["type"] == "done":
                    self.processing_complete()

        except queue.Empty:
            pass

        # Schedule next poll
        self.root.after(100, self.poll_queue)

    def update_totals(self):
        """Update total statistics."""
        total_nuclei = 0
        total_matched = 0
        total_unmatched = 0

        for item in self.results_tree.get_children():
            values = self.results_tree.item(item)["values"]
            total_nuclei += int(values[1])
            total_matched += int(values[2])
            total_unmatched += int(values[3])

        self.total_label.configure(
            text=f"Total: {total_nuclei} nuclei, {total_matched} matched, {total_unmatched} unmatched"
        )

    def processing_complete(self):
        """Handle completion of processing."""
        self.processing = False
        self.start_button.configure(state=tk.NORMAL)
        self.cancel_button.configure(state=tk.DISABLED)
        self.progress_label.configure(text="Processing complete")
        self.eta_label.configure(text="")

        self.log("=" * 60)
        self.log("Processing complete")
        self.log(f"Results saved to: {self.output_dir}")

        messagebox.showinfo("Complete", "Processing finished successfully!")


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================


def main():
    root = tk.Tk()
    app = CellposeSAMGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
