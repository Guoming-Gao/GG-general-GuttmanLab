# nd2MIPtif.py

import numpy as np
import tifffile
from pathlib import Path
from tkinter import filedialog as fd
from rich.progress import track
import shutil
import nd2


def get_nd2_resolution(nd2_file):
    """Extract pixel size from ND2 metadata"""
    try:
        # Get pixel calibration from metadata
        # ND2 files store this in metadata attributes
        pixel_size_x = nd2_file.metadata.channels[0].volume.axesCalibration[
            0
        ]  # X axis in microns
        pixel_size_y = nd2_file.metadata.channels[0].volume.axesCalibration[
            1
        ]  # Y axis in microns
        return (pixel_size_x, pixel_size_y)
    except (AttributeError, IndexError, KeyError):
        print(
            f"Warning: Could not extract pixel size from {nd2_file.path.name}, using default 1.0 µm"
        )
        return (1.0, 1.0)


def main():
    # Ask user to select ND2 files
    file_names = fd.askopenfilenames(
        title="Select ND2 files for batch processing",
        filetypes=[("ND2 files", "*.nd2")],
    )

    if not file_names:
        print("No files selected. Exiting.")
        return

    # Convert to Path objects
    file_paths = [Path(f) for f in file_names]

    # Get the parent directory of the first file
    parent_dir = file_paths[0].parent

    # Create output directories
    composite_dir = parent_dir / "composite"
    rawdata_dir = parent_dir / "rawdata"

    composite_dir.mkdir(exist_ok=True)
    rawdata_dir.mkdir(exist_ok=True)

    print(f"\nProcessing {len(file_paths)} ND2 file(s)...")
    print(f"Output directory: {parent_dir}")

    # Process each file
    for file_path in track(file_paths, description="Processing files"):
        try:
            # Open ND2 file
            with nd2.ND2File(file_path) as nd2_file:
                # Get resolution
                pixel_size_x, pixel_size_y = get_nd2_resolution(nd2_file)

                # Get image dimensions
                # ND2 file shape is typically (T, Z, C, Y, X) or subset thereof
                img_array = nd2_file.asarray()

                # Determine the shape and handle different dimension orders
                if img_array.ndim == 5:  # (T, Z, C, Y, X)
                    # Take first timepoint if multiple exist
                    img_array = img_array[0]  # Now (Z, C, Y, X)
                elif img_array.ndim == 4:  # Could be (Z, C, Y, X) or (T, C, Y, X)
                    # Check metadata to determine if first dim is T or Z
                    if (
                        hasattr(nd2_file, "sizes")
                        and "T" in nd2_file.sizes
                        and nd2_file.sizes["T"] > 1
                    ):
                        img_array = img_array[0]  # Take first timepoint
                    # Otherwise assume it's already (Z, C, Y, X)
                elif img_array.ndim == 3:  # (C, Y, X) - single Z plane
                    img_array = np.expand_dims(
                        img_array, axis=0
                    )  # Add Z dimension: (1, C, Y, X)

                # At this point, img_array should be (Z, C, Y, X)
                num_z, num_channels, height, width = img_array.shape

                # Create MIP for each channel
                channel_mips = []
                for c in range(num_channels):
                    channel_data = img_array[:, c, :, :]  # (Z, Y, X)
                    mip = np.max(channel_data, axis=0)  # Max projection along Z
                    channel_mips.append(mip)

                # Stack channels in reverse order (to match original script behavior)
                composite_stack = np.array(channel_mips[::-1])  # (C, Y, X)

                # Prepare metadata for TIFF
                resolution_dpi = (10000 / pixel_size_x, 10000 / pixel_size_y)
                metadata = {
                    "axes": "CYX",
                    "Channel": {"Name": [f"Channel{i}" for i in range(num_channels)]},
                }

                # Save composite TIFF
                output_filename = file_path.stem + "_composite.tif"
                output_path = composite_dir / output_filename

                tifffile.imwrite(
                    output_path,
                    composite_stack,
                    resolution=resolution_dpi,
                    metadata=metadata,
                    imagej=True,
                )

            # Move original ND2 file to rawdata folder
            destination = rawdata_dir / file_path.name
            shutil.move(str(file_path), str(destination))

        except Exception as e:
            print(f"\nError processing {file_path.name}: {str(e)}")
            continue

    print(f"\n✓ Processing complete!")
    print(f"  Composite images saved to: {composite_dir}")
    print(f"  Original ND2 files moved to: {rawdata_dir}")


if __name__ == "__main__":
    main()
