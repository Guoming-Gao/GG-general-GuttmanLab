import pandas as pd
import os
import tkinter as tk
from tkinter import filedialog, messagebox


def file_to_fasta(input_file, output_file):
    """
    Reads a file (Excel or CSV), extracts columns "ProbesNames" and "Seq",
    and writes them to a FASTA file.
    """
    # Read the file based on extension
    if input_file.endswith(".xlsx"):
        df = pd.read_excel(input_file)
    elif input_file.endswith(".csv"):
        df = pd.read_csv(input_file)
    else:
        raise ValueError("File must be .xlsx or .csv")

    # Ensure required columns are present
    required_columns = ["ProbesNames", "Seq"]
    if not all(col in df.columns for col in required_columns):
        raise ValueError(f"Input file must contain columns: {required_columns}")

    # Write FASTA file
    with open(output_file, "w") as fasta_file:
        for _, row in df.iterrows():
            title = f"{row['ProbesNames']}"
            sequence = str(row["Seq"])
            fasta_file.write(f">{title}\n")
            fasta_file.write(f"{sequence}\n")


# Simple file selection
root = tk.Tk()
root.withdraw()

input_file = filedialog.askopenfilename(title="Select your Excel or CSV file")

if input_file:
    # Generate output filename: remove extension, add .fasta
    output_file = os.path.splitext(input_file)[0] + ".fasta"

    try:
        file_to_fasta(input_file, output_file)
        print(f"Success! FASTA file saved as: {output_file}")
    except Exception as e:
        print(f"Error: {e}")
else:
    print("No file selected.")
