import pandas as pd
import os


def xlsx_to_fasta(input_xlsx, output_fasta):
    """
    Reads an Excel file, extracts columns "#", "ProbesNames", and "Seq",
    and writes them to a FASTA file with the specified title format.

    Parameters:
        input_xlsx (str): Path to the input Excel (.xlsx) file.
        output_fasta (str): Path to the output FASTA file.
    """
    # Read the Excel file
    df = pd.read_excel(input_xlsx)

    # Ensure required columns are present
    required_columns = ["ProbesNames", "Seq"]
    if not all(col in df.columns for col in required_columns):
        raise ValueError(f"Input Excel file must contain columns: {required_columns}")

    # Open the output FASTA file for writing
    with open(output_fasta, "w") as fasta_file:
        # Iterate through each row in the DataFrame
        for _, row in df.iterrows():
            # Create the FASTA title and sequence
            title = f"{row['ProbesNames']}"  # Title format: ProbesNames-#
            sequence = str(row["Seq"])  # Ensure sequence is a string

            # Write to the FASTA file
            fasta_file.write(f">{title}\n")  # Write the title line
            fasta_file.write(f"{sequence}\n")  # Write the sequence line


def csv_to_fasta(input_csv, output_fasta):
    """
    Reads an Excel file, extracts columns "#", "ProbesNames", and "Seq",
    and writes them to a FASTA file with the specified title format.

    Parameters:
        input_xlsx (str): Path to the input Excel (.xlsx) file.
        output_fasta (str): Path to the output FASTA file.
    """
    # Read the Excel file
    df = pd.read_csv(input_csv)

    # Ensure required columns are present
    required_columns = ["ProbesNames", "Seq"]
    if not all(col in df.columns for col in required_columns):
        raise ValueError(f"Input Excel file must contain columns: {required_columns}")

    # Open the output FASTA file for writing
    with open(output_fasta, "w") as fasta_file:
        # Iterate through each row in the DataFrame
        for _, row in df.iterrows():
            # Create the FASTA title and sequence
            title = f"{row['ProbesNames']}"  # Title format: ProbesNames-#
            sequence = str(row["Seq"])  # Ensure sequence is a string

            # Write to the FASTA file
            fasta_file.write(f">{title}\n")  # Write the title line
            fasta_file.write(f"{sequence}\n")  # Write the sequence line


folder = "/Users/gmgao/Dropbox/Caltech_PostDoc_GuttmanLab/constructs_and_smiFISH/Malat1"
fname_xlsx = "CandidateProbes_humanMalat1.xlsx"
# fname_csv = "Probes_C9orf72-final.csv"
os.chdir(folder)
xlsx_to_fasta(fname_xlsx, fname_xlsx[:-5] + ".fasta")
# csv_to_fasta(fname_csv, fname_csv[:-4] + ".fasta")
