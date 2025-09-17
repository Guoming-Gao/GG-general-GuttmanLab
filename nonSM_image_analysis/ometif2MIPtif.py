# ometif2MIPtif.py

import numpy as np
import tifffile
from pathlib import Path
from tkinter import filedialog as fd
from rich.progress import track
import shutil

print("Choose the mother folder containing subfolders with OME-TIFF files:")

# Select mother folder instead of individual files
mother_folder = fd.askdirectory(title="Select Mother Folder (e.g., 20250916)")

if not mother_folder:
    exit()

mother_path = Path(mother_folder)

# Find all .ome.tif files in immediate subfolders
ome_tif_files = []
subfolders_with_files = set()

for subfolder in mother_path.iterdir():
    if subfolder.is_dir():
        # Look for .ome.tif files directly in this subfolder
        ome_files = list(subfolder.glob("*.ome.tif"))
        if ome_files:
            ome_tif_files.extend(ome_files)
            subfolders_with_files.add(subfolder)

if not ome_tif_files:
    print("No OME-TIFF files found in subfolders. Exiting.")
    exit()

print(f"Found {len(ome_tif_files)} OME-TIFF files to process")

# Fixed pixel size: 103 nm = 0.103 μm per pixel
pixel_size_um = 0.103
resolution = (1 / pixel_size_um, 1 / pixel_size_um)

for fpath in track(ome_tif_files, description="Processing files"):
    try:
        # Read OME-TIFF file
        img_stack = tifffile.imread(fpath)

        # CRITICAL FIX: Reverse channel order to match CZI convention
        # OME-TIFF: wavelength small to large (DAPI first)
        # CZI convention: wavelength large to small (DAPI last)
        if img_stack.ndim == 4:  # (C, Z, Y, X)
            img_stack = img_stack[::-1]  # Reverse channel order
            num_channels = img_stack.shape[0]
        elif img_stack.ndim == 3:
            if img_stack.shape[0] <= 10:  # Likely (C, Y, X) multi-channel
                img_stack = img_stack[::-1]  # Reverse channel order
                num_channels = img_stack.shape[0]
            else:  # Likely (Z, Y, X) single channel
                num_channels = 1
        else:
            print(f"Unexpected image dimensions in {fpath.name}. Skipping.")
            continue

        mip_list = []

        for channel_idx in range(num_channels):
            # Extract channel data
            if img_stack.ndim == 3:
                if num_channels == 1:
                    squeezed_stack = img_stack
                else:
                    squeezed_stack = img_stack[channel_idx]
            else:  # ndim == 4
                squeezed_stack = img_stack[channel_idx]

            # Ensure we have a proper stack for MIP
            if squeezed_stack.ndim == 2:
                # Single slice, just use as is
                max_proj = squeezed_stack
            else:
                # Multiple Z slices, create MIP
                max_proj = np.max(squeezed_stack, axis=0)

            # Create channel folder
            channel_folder = mother_path / f"ch{channel_idx}"
            channel_folder.mkdir(exist_ok=True)

            # Save single MIP
            tifffile.imwrite(
                channel_folder / f"{fpath.stem}_ch{channel_idx}_MIP.tif",
                max_proj,
                imagej=True,
                resolution=resolution,
            )

            mip_list.append(max_proj)

        # Create composite and average projection folders
        composite_folder = mother_path / "composite"
        composite_folder.mkdir(exist_ok=True)
        ave_proj_folder = mother_path / "ave_proj"
        ave_proj_folder.mkdir(exist_ok=True)

        # Create composite stack (reverse order to match original script)
        composite_stack = np.stack(mip_list[::-1], axis=0)

        tifffile.imwrite(
            composite_folder / f"{fpath.stem}_composite_stack.tif",
            composite_stack.astype(np.float32),
            imagej=True,
            metadata={"axes": "CYX"},
            resolution=resolution,
        )

        # Create average projection (exclude DAPI channel - first channel)
        if len(mip_list) > 1:
            ave_proj = np.mean(composite_stack[1:], axis=0)  # Skip first (DAPI) channel
        else:
            ave_proj = composite_stack[0]  # If only one channel, use it

        tifffile.imwrite(
            ave_proj_folder / f"{fpath.stem}_ave_proj.tif",
            ave_proj.astype(np.float32),
            imagej=True,
            resolution=resolution,
        )

    except Exception as e:
        print(f"Error processing {fpath.name}: {e}")
        continue

# Move original subfolders to rawdata folder
rawdata_folder = mother_path / "rawdata"
rawdata_folder.mkdir(exist_ok=True)

# Define result folder names to exclude from moving
result_folders = {"ch0", "ch1", "ch2", "ch3", "composite", "ave_proj", "rawdata"}

print("Moving original subfolders to rawdata...")
for subfolder in subfolders_with_files:
    if subfolder.name not in result_folders:
        try:
            shutil.move(str(subfolder), str(rawdata_folder / subfolder.name))
        except Exception as e:
            print(f"Error moving {subfolder.name}: {e}")

print("Processing complete!")
