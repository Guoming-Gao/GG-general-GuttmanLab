from tifffile import imread, imwrite
import numpy as np
from tkinter import filedialog as fd
from rich.progress import track


#########################################
# Load and organize files
print("Choose the tif files for crop")
lst_files = list(fd.askopenfilenames())


#########################################
# Apply registration and Crop
for fpath in track(lst_files):
    # load the tiff file
    video = imread(fpath)
    halfwidth = int(video.shape[2] / 2)

    print(video.shape, fpath)
    img1 = video[0, :, 0:halfwidth]
    img2 = video[3, :, 0:halfwidth]
    img3 = video[6, :, halfwidth:]

    imwrite(
        fpath[:-4] + "-composite.tif",
        np.stack([img1, img2, img3], axis=0),
        imagej=True,
        metadata={"axes": "CYX"},
    )
