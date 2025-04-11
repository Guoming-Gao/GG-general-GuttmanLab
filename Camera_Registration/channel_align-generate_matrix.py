import os
import cv2
import pickle
from tifffile import imread, imwrite
from tkinter import filedialog as fd
import numpy as np


# For ONI:
# path = fd.askopenfilename(title="Choose a combined bead image:")
# img = imread(path)
# img_left = img[:, 0:428]
# img_right = img[:, 428:]
# imwrite(path[:-4] + "-left.tif", img_left, imagej=True)
# imwrite(path[:-4] + "-right.tif", img_right, imagej=True)
# img_left = cv2.imread(path[:-4] + "-left.tif", cv2.IMREAD_GRAYSCALE)
# img_right = cv2.imread(path[:-4] + "-right.tif", cv2.IMREAD_GRAYSCALE)


# For mannual chosen images
path_red = fd.askopenfilename(title="Choose the red bead image:")
img_red = cv2.imread(path_red, cv2.IMREAD_GRAYSCALE)
path_green = fd.askopenfilename(title="Choose the green bead image:")
img_green = cv2.imread(path_green, cv2.IMREAD_GRAYSCALE)

# Find size of image1
sz = img_red.shape

# Define the motion model
warp_mode = cv2.MOTION_HOMOGRAPHY

# Define a 3x3 matrices and initialize the matrix to identity
warp_matrix_green = np.eye(3, 3, dtype=np.float32)
warp_matrix_red = np.eye(3, 3, dtype=np.float32)

# Specify the number of iterations.
number_of_iterations = 100

# Specify the threshold of the increment
# in the correlation coefficient between two iterations
termination_eps = 1e-5

# Define termination criteria
criteria = (
    cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT,
    number_of_iterations,
    termination_eps,
)

# Run the ECC algorithm. The results are stored in warp_matrix.
(cc, warp_matrix_green) = cv2.findTransformECC(
    img_red, img_green, warp_matrix_green, warp_mode, criteria
)
(cc, warp_matrix_red) = cv2.findTransformECC(
    img_green, img_red, warp_matrix_red, warp_mode, criteria
)

# Save the matrix
folderpath = os.path.dirname(path_red)
pickle.dump(
    warp_matrix_green, open(os.path.join(folderpath, "ONI_warp_matrix_green.p"), "wb")
)
pickle.dump(
    warp_matrix_red, open(os.path.join(folderpath, "ONI_warp_matrix_red.p"), "wb")
)

# Test and save alignment results on beads
img_green_unit16 = cv2.imread(path_green, -1)
img_green_aligned = cv2.warpPerspective(
    img_green_unit16,
    warp_matrix_green,
    (sz[1], sz[0]),
    flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP,
)
imwrite(path_green[:-4] + "-aligned.tif", img_green_aligned, imagej=True)

img_red_unit16 = cv2.imread(path_red, -1)
img_red_aligned = cv2.warpPerspective(
    img_red_unit16,
    warp_matrix_red,
    (sz[1], sz[0]),
    flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP,
)
imwrite(path_red[:-4] + "-aligned.tif", img_red_aligned, imagej=True)
