import os
import numpy as np
import tifffile
from pathlib import Path
from tkinter import filedialog as fd
from tkinter import Tk
from os.path import commonprefix


def ask_for_first_channel_files():
    """Open a file dialog to select TIF files for the first channel."""
    print("Select TIF files for the first channel:")
    root = Tk()
    root.withdraw()  # Hide the main window
    first_channel_files = list(
        fd.askopenfilenames(
            title="Select First Channel TIF Files", filetypes=[("TIF Files", "*.tif")]
        )
    )
    root.destroy()

    if not first_channel_files:
        print("No files selected for the first channel.")
        return None
    return first_channel_files


def ask_for_second_channel_files():
    """Open a file dialog to select TIF files for the second channel."""
    print("Select TIF files for the second channel:")
    root = Tk()
    root.withdraw()  # Hide the main window
    second_channel_files = list(
        fd.askopenfilenames(
            title="Select Second Channel TIF Files", filetypes=[("TIF Files", "*.tif")]
        )
    )
    root.destroy()

    if not second_channel_files:
        print("No files selected for the second channel.")
        return None
    return second_channel_files


def find_matching_pairs(first_channel_files, second_channel_files):
    """Find matching file pairs based on common filename prefixes."""
    pairs = []

    # Get file stems (filenames without extension)
    first_file_stems = {Path(f).stem: f for f in first_channel_files}
    second_file_stems = {Path(f).stem: f for f in second_channel_files}

    # Identify potential pairs
    for first_stem, first_path in first_file_stems.items():
        best_match = None
        best_match_len = 0

        for second_stem, second_path in second_file_stems.items():
            # Find common prefix between stems
            prefix = commonprefix([first_stem, second_stem])

            # If prefix is significant (more than 3 characters) and better than previous matches
            if len(prefix) > 3 and len(prefix) > best_match_len:
                best_match = second_path
                best_match_len = len(prefix)

        # If a good match was found, add to pairs
        if best_match:
            pairs.append((first_path, best_match))
            # Remove the matched second channel file from consideration for other pairs
            second_file_stems.pop(Path(best_match).stem)

    # Print pairing information
    if pairs:
        print(f"\nFound {len(pairs)} matching file pairs:")
        for i, (first, second) in enumerate(pairs):
            print(f"Pair {i+1}: {Path(first).name} with {Path(second).name}")
    else:
        print("No matching file pairs found based on filename prefixes.")

    return pairs


def create_composite_tif(file_pairs, output_postfix="-pleaserename.tif"):
    """Create composite TIF files from pairs of channel files."""
    if not file_pairs:
        print("No file pairs to process.")
        return

    for first_file, second_file in file_pairs:
        try:
            # Read the TIF files
            print(f"Reading {Path(first_file).name} and {Path(second_file).name}...")
            first_image = tifffile.imread(first_file)
            second_image = tifffile.imread(second_file)

            # Check if images have compatible dimensions
            if first_image.shape != second_image.shape:
                print(
                    f"Warning: Image dimensions don't match for {Path(first_file).name} and {Path(second_file).name}"
                )
                continue

            # Create a composite image by stacking
            composite_image = np.stack([first_image, second_image], axis=0)

            # Define output path - parent directory of the first channel file
            first_file_path = Path(first_file)
            output_dir = first_file_path.parent.parent
            output_name = first_file_path.stem + output_postfix
            output_path = output_dir / output_name

            # Try to get resolution from first image
            try:
                with tifffile.TiffFile(first_file) as tif:
                    tags = tif.pages[0].tags
                    if "XResolution" in tags and "YResolution" in tags:
                        x_res = tags["XResolution"].value
                        y_res = tags["YResolution"].value
                        resolution = (x_res, y_res)
                        resolution_unit = tags.get("ResolutionUnit", None)

                        # Save with resolution metadata
                        tifffile.imwrite(
                            str(output_path),
                            composite_image,
                            imagej=True,
                            resolution=resolution,
                            resolutionunit=(
                                resolution_unit.value if resolution_unit else None
                            ),
                            metadata={"axes": "CYX"},
                        )
                    else:
                        # Save without resolution metadata
                        tifffile.imwrite(
                            str(output_path),
                            composite_image,
                            imagej=True,
                            metadata={"axes": "CYX"},
                        )
            except Exception as e:
                print(f"Could not read resolution metadata: {e}")
                # Save without resolution metadata
                tifffile.imwrite(
                    str(output_path),
                    composite_image,
                    imagej=True,
                    metadata={"axes": "CYX"},
                )

            print(f"Saved composite TIF: {output_path}")

        except Exception as e:
            print(
                f"Error processing {Path(first_file).name} and {Path(second_file).name}: {e}"
            )


def main():
    """Main script execution."""
    # Step 1: Ask for first channel files
    first_channel_files = ask_for_first_channel_files()
    if not first_channel_files:
        return

    # Step 2: Ask for second channel files
    second_channel_files = ask_for_second_channel_files()
    if not second_channel_files:
        return

    # Step 3: Find matching pairs based on common prefixes
    file_pairs = find_matching_pairs(first_channel_files, second_channel_files)

    # Step 4: Create composite TIF files
    create_composite_tif(file_pairs)

    print("Processing complete!")


if __name__ == "__main__":
    main()
