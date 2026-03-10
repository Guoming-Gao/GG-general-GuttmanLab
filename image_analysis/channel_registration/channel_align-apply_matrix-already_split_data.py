from tifffile import imread, imwrite
import cv2
import pickle
import numpy as np
from tkinter import filedialog as fd
from rich.progress import track
import os


def find_common_pairs(file_names):
    """Find pairs of red/green files with matching base names."""
    red_files = {f for f in file_names if f.endswith("red.tif")}
    green_files = {f for f in file_names if f.endswith("green.tif")}

    pairs = []
    for red in red_files:
        base = red.rsplit("red.tif", 1)[0]
        green = f"{base}green.tif"
        if green in green_files:
            pairs.append((red, green))
    return pairs


def crop_image(img):
    """Handle both 2D and 3D images (t, h, w) or (h, w)"""
    if img.ndim == 3:
        t, h, w = img.shape
        return img[:, 5 : h - 5, 5 : w - 5]
    elif img.ndim == 2:
        h, w = img.shape
        return img[5 : h - 5, 5 : w - 5]
    else:
        raise ValueError("Unsupported image dimensions")


def transform(img, warp_matrix):
    """Handle both 2D and 3D images"""
    if img.ndim == 3:
        return np.stack(
            [
                cv2.warpPerspective(
                    frame,
                    warp_matrix,
                    (img.shape[2], img.shape[1]),
                    flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP,
                )
                for frame in img
            ],
            axis=0,
        )
    elif img.ndim == 2:
        return cv2.warpPerspective(
            img,
            warp_matrix,
            (img.shape[1], img.shape[0]),
            flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP,
        )
    else:
        raise ValueError("Unsupported image dimensions")


#########################################
# Main
print("Pick the channel to perform transformation, 1 for green, 2 for red:")
selector = input()
if selector == "1":
    print("Choose the matrix for green")
elif selector == "2":
    print("Choose the matrix for red")
else:
    print("Please enter either 1 or 2")
    exit()

path_matrix = fd.askopenfilename()
warp_matrix = pickle.load(open(path_matrix, "rb"))

print("Choose the tif files for channel alignment:")
lst_files = list(fd.askopenfilenames())
pairs = find_common_pairs(lst_files)

for red_path, green_path in track(pairs):
    img_green = imread(green_path)
    img_red = imread(red_path)

    if selector == "1":
        img_green_aligned = transform(img_green, warp_matrix)
        img_red_aligned = img_red
    elif selector == "2":
        img_red_aligned = transform(img_red, warp_matrix)
        img_green_aligned = img_green

    # Process both channels with dimension-agnostic functions
    img_green_cropped = crop_image(img_green_aligned)
    img_red_cropped = crop_image(img_red_aligned)

    # Save with appropriate metadata
    fsave_green = green_path.replace(".tif", "-aligned.tif")
    fsave_red = red_path.replace(".tif", "-aligned.tif")

    imwrite(
        fsave_green,
        img_green_cropped,
        imagej=True,
        metadata={"axes": "TYX" if img_green_cropped.ndim == 3 else "YX"},
    )

    imwrite(
        fsave_red,
        img_red_cropped,
        imagej=True,
        metadata={"axes": "TYX" if img_red_cropped.ndim == 3 else "YX"},
    )
