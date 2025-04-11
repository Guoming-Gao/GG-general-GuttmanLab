import numpy as np
from aicspylibczi import CziFile
import tifffile
from pathlib import Path
from tkinter import filedialog as fd
from rich.progress import track
import shutil


print("Choose all CZI files to be batch processed:")
file_paths = list(
    fd.askopenfilenames(title="Select CZI Files", filetypes=[("CZI Files", "*.czi")])
)

if not file_paths:
    exit()

parent_dir = Path(file_paths[0]).parent

for fpath in track(file_paths, description="Processing files"):
    czi = CziFile(fpath)
    dimensions = czi.get_dims_shape()
    num_channels = dimensions[0]["C"][-1]

    # List to collect MIPs for composite stack
    mip_list = []

    for channel_idx in range(num_channels):
        # Process individual channel
        channel_folder = parent_dir / f"ch{channel_idx}"
        channel_folder.mkdir(exist_ok=True)

        single_color_stack, _ = czi.read_image(C=channel_idx)
        max_proj = np.max(single_color_stack.squeeze(), axis=0)

        # Save individual channel MIP
        output_name = f"{Path(fpath).stem}_ch{channel_idx}_MIP.tif"
        output_path = channel_folder / output_name
        tifffile.imwrite(output_path, max_proj, imagej=True)

        # Collect for composite stack
        mip_list.append(max_proj)

    # Reverse the order of channels for blue-to-red stacking
    # The Zessi microscope saves channels from red (long wavelength) to blue (short wavelength)
    mip_list.reverse()  # Reverse order: blue → green → red

    # Save composite stack (channels now ordered blue-to-red)
    if mip_list:
        composite_stack = np.stack(mip_list, axis=0)  # Shape: (C, H, W)
        composite_name = f"{Path(fpath).stem}_composite_stack.tif"
        composite_path = parent_dir / composite_name
        tifffile.imwrite(
            composite_path,
            composite_stack.astype(np.float32),
            imagej=True,
            metadata={"axes": "CYX"},  # Explicit channel axis for ImageJ
            resolution=(
                1.0 / 0.108,
                1.0 / 0.108,
            ),  # Example 0.108 µm/pixel, adjust as needed
        )

# Move original CZI files to rawdata folder
rawdata_folder = parent_dir / "rawdata"
rawdata_folder.mkdir(exist_ok=True)

for fpath in file_paths:
    destination = rawdata_folder / Path(fpath).name
    shutil.move(fpath, destination)
