from tifffile import imread, imwrite
from scipy.ndimage import white_tophat
from skimage.filters import difference_of_gaussians
from skimage.util import img_as_uint
from tkinter import filedialog as fd
import numpy as np
from rich.progress import track
from os.path import dirname, basename, join
import os
import shutil


def create_output_folder(parent_dir, sigma_low, sigma_high):
    """Create output folder with sigma parameters in name"""
    folder_name = (
        f"{basename(parent_dir)}-bandpass-sigma1-{sigma_low}-sigma2-{sigma_high}"
    )
    output_path = join(parent_dir, folder_name)

    try:
        os.mkdir(output_path)
    except FileExistsError:
        shutil.rmtree(output_path)
        os.mkdir(output_path)
    return output_path


def process_video(file_path, sigma_low, sigma_high, output_folder):
    """Process individual video file with given DoG parameters"""
    video = imread(file_path)
    video_out = []

    for img in video:
        img_filtered = difference_of_gaussians(img, sigma_low, sigma_high)
        video_out.append(img_filtered)

    video_out = np.stack(video_out)
    base_name = f"{basename(file_path)[:-4]}-bandpass-sigma1-{sigma_low}-sigma2-{sigma_high}.tif"
    save_path = join(output_folder, base_name)

    # Handle different dimensional cases
    if len(video.shape) < 3:  # Single image
        imwrite(save_path, img_as_uint(video_out), imagej=True)
    else:  # Video stack
        imwrite(
            save_path, img_as_uint(video_out), imagej=True, metadata={"axes": "TYX"}
        )


def batch_process(param_sets):
    """Main batch processing function handling multiple parameter sets"""
    print("Choose all tif files to be batch processed:")
    file_list = list(fd.askopenfilenames())

    if not file_list:
        return  # Exit if no files selected

    parent_dir = dirname(file_list[0])

    for sigma_low, sigma_high in param_sets:
        output_folder = create_output_folder(parent_dir, sigma_low, sigma_high)

        for file_path in track(
            file_list, description=f"Processing σ={sigma_low},{sigma_high}"
        ):
            process_video(file_path, sigma_low, sigma_high, output_folder)


if __name__ == "__main__":
    # Define your parameter sets here
    PARAMETER_SETS = [
        (2, 5),
        # (1, 3),
    ]

    batch_process(PARAMETER_SETS)
