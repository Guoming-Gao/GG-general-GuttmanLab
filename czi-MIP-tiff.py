import numpy as np
from aicspylibczi import CziFile
import tifffile
from pathlib import Path
from tkinter import filedialog as fd
from rich.progress import track


print("Choose all tif files to be batch proccessed:")
file_paths = list(
    fd.askopenfilenames(title="Select CZI Files", filetypes=[("CZI Files", "*.czi")])
)

# Process each selected file
for fpath in track(file_paths):
    czi = CziFile(fpath)
    # Get dimensions and shape
    dimensions = czi.get_dims_shape()  # Returns dimension names and sizes
    num_channels = dimensions[0]["C"][-1]
    # num_zslices = dimensions[0]['Z'][-1]

    # Read image data for a specific channel (C=0) and Z slice (Z=0)
    parent_dir = Path(fpath).parent

    for channel_idx in range(num_channels):
        # Read image data for a specific channel (C=0) and Z slice (Z=0)
        single_color_stack, _ = czi.read_image(C=channel_idx)
        max_proj = np.max(single_color_stack.squeeze(), axis=0)

        output_name = f"{Path(fpath).stem}_ch{channel_idx}_MIP.tif"
        output_path = parent_dir / output_name

        tifffile.imwrite(output_path, max_proj, imagej=True)
