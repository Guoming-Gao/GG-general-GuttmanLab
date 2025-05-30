from tifffile import imread, imwrite
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

    imwrite(
        fpath[:-4] + "-cropped-left.tif",
        video[:, :, 0:halfwidth],
        imagej=True,
        metadata={"axes": "TYX"},
    )
    imwrite(
        fpath[:-4] + "-cropped-right.tif",
        video[:, :, halfwidth:],
        imagej=True,
        metadata={"axes": "TYX"},
    )
