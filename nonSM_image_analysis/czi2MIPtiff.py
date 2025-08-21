# czi2MIPtiff.py
import numpy as np
import tifffile
from pathlib import Path
from tkinter import filedialog as fd
from rich.progress import track
import shutil
import re

# czifile package has problem reading the image, aicspylibczi has problem reading the meatadata, and that's why we need both
from czifile import CziFile
from aicspylibczi import CziFile as czi_ai

N_threshold = 45


def get_czi_resolution(metadata_string, num_xpixels, num_ypixels):
    """Extract pixel size using simple string methods"""
    # Find the ScaledImageRectangleSize element using regex pattern
    pattern = r'<ScaledImageRectangleSize Status="Valid">([\d\.,]+)</ScaledImageRectangleSize>'
    match = re.search(pattern, metadata_string)

    if match:
        # Extract values and convert to floats
        fov_sizes = match.group(1).split(",")
        fov_x = float(fov_sizes[0])
        fov_y = float(fov_sizes[1])

        # Calculate pixel size
        pixel_size_x = fov_x / num_xpixels  # in microns
        pixel_size_y = fov_y / num_ypixels  # in microns

        return (pixel_size_x, pixel_size_y)
    else:
        # Fallback to simple string search if regex fails
        start_tag = '<ScaledImageRectangleSize Status="Valid">'
        end_tag = "</ScaledImageRectangleSize>"

        start_idx = metadata_string.find(start_tag)
        if start_idx != -1:
            start_idx += len(start_tag)
            end_idx = metadata_string.find(end_tag, start_idx)
            if end_idx != -1:
                content = metadata_string[start_idx:end_idx]
                fov_x, fov_y = map(float, content.split(","))

                # Calculate pixel size
                pixel_size_x = fov_x / num_xpixels  # in microns
                pixel_size_y = fov_y / num_ypixels  # in microns

                return (pixel_size_x, pixel_size_y)

    raise ValueError("Could not find ScaledImageRectangleSize in metadata")


print("Choose all CZI files to be batch processed:")
file_paths = fd.askopenfilenames(
    title="Select CZI Files", filetypes=[("CZI Files", "*.czi")]
)

if not file_paths:
    exit()

parent_dir = Path(file_paths[0]).parent

for fpath in track(file_paths, description="Processing files"):
    czi_for_img = czi_ai(fpath)

    with CziFile(fpath) as czi:
        # Get metadata as a string and extract resolution
        metadata_string = czi.metadata()

        # Determine number of channels
        axes = czi.axes
        num_channels = czi.shape[axes.index("C")] if "C" in axes else 1
        num_xpixels = czi.shape[axes.index("X")]
        num_ypixels = czi.shape[axes.index("Y")]

        # Extract resolution using the simple string approach
        try:
            x_micron_per_pxl, y_micron_per_pxl = get_czi_resolution(
                metadata_string, num_xpixels, num_ypixels
            )
        except Exception as e:
            print(f"Critical error in {Path(fpath).name}: {e}")
            print("Skipping file due to missing essential metadata")
            continue

        # Check if any channel needs sectioned processing
        has_large_z_stack = False
        for channel_idx in range(num_channels):
            single_color_stack, _ = czi_for_img.read_image(C=channel_idx)
            squeezed_stack = single_color_stack.squeeze()
            if squeezed_stack.shape[0] > N_threshold:
                has_large_z_stack = True
                break

        mip_list_top = []
        mip_list_middle = []
        mip_list_bottom = []

        for channel_idx in range(num_channels):
            # Read and process channel data
            single_color_stack, _ = czi_for_img.read_image(C=channel_idx)
            squeezed_stack = single_color_stack.squeeze()
            num_z_slices = squeezed_stack.shape[0]

            channel_folder = parent_dir / f"ch{channel_idx}"
            channel_folder.mkdir(exist_ok=True)

            if num_z_slices > N_threshold:
                # Create three sectioned MIPs
                # Top 25 slices
                top_mip = np.max(squeezed_stack[0:25], axis=0)
                tifffile.imwrite(
                    channel_folder / f"{Path(fpath).stem}_ch{channel_idx}_MIP_top.tif",
                    top_mip,
                    imagej=True,
                    resolution=(1 / x_micron_per_pxl, 1 / y_micron_per_pxl),
                )

                # Bottom 25 slices
                bottom_mip = np.max(squeezed_stack[-25:], axis=0)
                tifffile.imwrite(
                    channel_folder
                    / f"{Path(fpath).stem}_ch{channel_idx}_MIP_bottom.tif",
                    bottom_mip,
                    imagej=True,
                    resolution=(1 / x_micron_per_pxl, 1 / y_micron_per_pxl),
                )

                # Middle 25 slices
                middle_start = (num_z_slices - 25) // 2
                middle_end = middle_start + 25
                middle_mip = np.max(squeezed_stack[middle_start:middle_end], axis=0)
                tifffile.imwrite(
                    channel_folder
                    / f"{Path(fpath).stem}_ch{channel_idx}_MIP_middle.tif",
                    middle_mip,
                    imagej=True,
                    resolution=(1 / x_micron_per_pxl, 1 / y_micron_per_pxl),
                )

                mip_list_top.append(top_mip)
                mip_list_bottom.append(bottom_mip)
                mip_list_middle.append(middle_mip)
            else:
                max_proj = np.max(squeezed_stack, axis=0)
                # Save single MIP
                tifffile.imwrite(
                    channel_folder / f"{Path(fpath).stem}_ch{channel_idx}_MIP.tif",
                    max_proj,
                    imagej=True,
                    resolution=(1 / x_micron_per_pxl, 1 / y_micron_per_pxl),
                )

                # Append same max_proj to all mip_lists so it's consistent for composite and average
                mip_list_top.append(max_proj)
                mip_list_bottom.append(max_proj)
                mip_list_middle.append(max_proj)

        # Use original composite and ave_proj folders
        composite_folder = parent_dir / "composite"
        composite_folder.mkdir(exist_ok=True)
        ave_proj_folder = parent_dir / "ave_proj"
        ave_proj_folder.mkdir(exist_ok=True)

        if has_large_z_stack:
            # Create sectioned composite and average projections
            if mip_list_top:
                composite_stack_top = np.stack(mip_list_top[::-1], axis=0)
                tifffile.imwrite(
                    composite_folder / f"{Path(fpath).stem}_composite_stack_top.tif",
                    composite_stack_top.astype(np.float32),
                    imagej=True,
                    metadata={"axes": "CYX"},
                    resolution=(1 / x_micron_per_pxl, 1 / y_micron_per_pxl),
                )

                ave_proj_top = np.mean(
                    composite_stack_top[1:], axis=0
                )  # only non DAPI channel
                tifffile.imwrite(
                    ave_proj_folder / f"{Path(fpath).stem}_ave_proj_top.tif",
                    ave_proj_top.astype(np.float32),
                    imagej=True,
                    resolution=(1 / x_micron_per_pxl, 1 / y_micron_per_pxl),
                )

            if mip_list_middle:
                composite_stack_middle = np.stack(mip_list_middle[::-1], axis=0)
                tifffile.imwrite(
                    composite_folder / f"{Path(fpath).stem}_composite_stack_middle.tif",
                    composite_stack_middle.astype(np.float32),
                    imagej=True,
                    metadata={"axes": "CYX"},
                    resolution=(1 / x_micron_per_pxl, 1 / y_micron_per_pxl),
                )

                ave_proj_middle = np.mean(
                    composite_stack_middle[1:], axis=0
                )  # only non DAPI channel
                tifffile.imwrite(
                    ave_proj_folder / f"{Path(fpath).stem}_ave_proj_middle.tif",
                    ave_proj_middle.astype(np.float32),
                    imagej=True,
                    resolution=(1 / x_micron_per_pxl, 1 / y_micron_per_pxl),
                )

            if mip_list_bottom:
                composite_stack_bottom = np.stack(mip_list_bottom[::-1], axis=0)
                tifffile.imwrite(
                    composite_folder / f"{Path(fpath).stem}_composite_stack_bottom.tif",
                    composite_stack_bottom.astype(np.float32),
                    imagej=True,
                    metadata={"axes": "CYX"},
                    resolution=(1 / x_micron_per_pxl, 1 / y_micron_per_pxl),
                )

                ave_proj_bottom = np.mean(
                    composite_stack_bottom[1:], axis=0
                )  # only non DAPI channel
                tifffile.imwrite(
                    ave_proj_folder / f"{Path(fpath).stem}_ave_proj_bottom.tif",
                    ave_proj_bottom.astype(np.float32),
                    imagej=True,
                    resolution=(1 / x_micron_per_pxl, 1 / y_micron_per_pxl),
                )
        else:
            # Create single composite and average projection (original behavior)
            if mip_list_top:  # All lists contain the same data for ≤N_threshold slices
                composite_stack = np.stack(mip_list_top[::-1], axis=0)
                tifffile.imwrite(
                    composite_folder / f"{Path(fpath).stem}_composite_stack.tif",
                    composite_stack.astype(np.float32),
                    imagej=True,
                    metadata={"axes": "CYX"},
                    resolution=(1 / x_micron_per_pxl, 1 / y_micron_per_pxl),
                )

                ave_proj = np.mean(composite_stack[1:], axis=0)  # only non DAPI channel
                tifffile.imwrite(
                    ave_proj_folder / f"{Path(fpath).stem}_ave_proj.tif",
                    ave_proj.astype(np.float32),
                    imagej=True,
                    resolution=(1 / x_micron_per_pxl, 1 / y_micron_per_pxl),
                )


# Move original files to rawdata folder
rawdata_folder = parent_dir / "rawdata"
rawdata_folder.mkdir(exist_ok=True)

for fpath in file_paths:
    shutil.move(fpath, rawdata_folder / Path(fpath).name)
