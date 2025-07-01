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

        mip_list = []
        for channel_idx in range(num_channels):
            # Read and process channel data
            single_color_stack, _ = czi_for_img.read_image(C=channel_idx)
            max_proj = np.max(single_color_stack.squeeze(), axis=0)

            # Save individual MIP
            channel_folder = parent_dir / f"ch{channel_idx}"
            channel_folder.mkdir(exist_ok=True)
            tifffile.imwrite(
                channel_folder / f"{Path(fpath).stem}_ch{channel_idx}_MIP.tif",
                max_proj,
                imagej=True,
                resolution=(1 / x_micron_per_pxl, 1 / y_micron_per_pxl),
            )
            mip_list.append(max_proj)

        # Create composite stack
        if mip_list:
            composite_stack = np.stack(mip_list[::-1], axis=0)
            tifffile.imwrite(
                parent_dir / f"{Path(fpath).stem}_composite_stack.tif",
                composite_stack.astype(np.float32),
                imagej=True,
                metadata={"axes": "CYX"},
                resolution=(1 / x_micron_per_pxl, 1 / y_micron_per_pxl),
            )

# Move original files to rawdata folder
rawdata_folder = parent_dir / "rawdata"
rawdata_folder.mkdir(exist_ok=True)
for fpath in file_paths:
    shutil.move(fpath, rawdata_folder / Path(fpath).name)
