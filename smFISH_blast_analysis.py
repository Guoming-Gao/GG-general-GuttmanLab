import pandas as pd
import re
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import os


def parse_blast_results(blast_text):
    """
    Parse BLAST result text and extract required information
    Returns DataFrame with columns: ProbeName, ProbeSequence, PercentAlignment,
    NumberOfHits, UniqueHitName, Start, End
    """
    # Split by queries
    queries = re.split(r"Query #\d+: ", blast_text)[1:]  # skip first empty split

    records = []

    for query in queries:
        # Extract full probe name (everything before "Query ID:")
        probe_name_search = re.search(r"^(.+?)\s+Query ID:", query)
        probe_name = probe_name_search.group(1).strip() if probe_name_search else None

        # CORRECTED: Count alignment sections starting with ">"
        # This counts actual genomic hits, not table entries
        alignment_headers = re.findall(r"^>([^\n]+)", query, re.MULTILINE)
        num_hits = len(alignment_headers)

        # Extract unique hit name if only one hit
        unique_hit_name = None
        if num_hits == 1:
            unique_hit_name = alignment_headers[0].strip()

        # Extract start and end positions (from first alignment)
        start_end_search = re.search(r"Range 1: (\d+) to (\d+)", query)
        start = int(start_end_search.group(1)) if start_end_search else None
        end = int(start_end_search.group(2)) if start_end_search else None

        # Extract percentage identity
        perc_identity_search = re.search(r"Identities:\s*\d+/\d+\((\d+)%\)", query)
        perc_identity = (
            int(perc_identity_search.group(1)) if perc_identity_search else None
        )

        # Extract probe sequence (longest query sequence from alignment)
        query_seqs = re.findall(r"Query\s+\d+\s+([A-Z]+)\s+\d+", query)
        probe_seq = max(query_seqs, key=len) if query_seqs else None

        records.append(
            {
                "ProbeName": probe_name,
                "ProbeSequence": probe_seq,
                "PercentAlignment": perc_identity,
                "NumberOfHits": num_hits,
                "UniqueHitName": unique_hit_name,
                "Start": start,
                "End": end,
            }
        )

    return pd.DataFrame(records)


def main():
    # Create GUI window (hidden)
    root = Tk()
    root.withdraw()  # Hide main window

    print("Starting BLAST results analysis...")

    # Step 1: Select BLAST result text file
    print("Please select your BLAST result text file...")
    blast_file_path = askopenfilename(
        title="Select BLAST result text file",
        filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
    )

    if not blast_file_path:
        print("No BLAST file selected. Exiting.")
        return

    # Read and parse BLAST results
    print(f"Reading BLAST results from: {blast_file_path}")
    with open(blast_file_path, "r") as f:
        blast_text = f.read()

    print("Parsing BLAST results...")
    blast_df = parse_blast_results(blast_text)
    print(f"Parsed {len(blast_df)} probe results")

    # Step 2: Select CSV file with ProbesNames
    print("Please select your CSV file with ProbesNames column...")
    csv_file_path = askopenfilename(
        title="Select CSV file with ProbesNames column",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
    )

    if not csv_file_path:
        print("No CSV file selected. Exiting.")
        return

    # Read CSV file
    print(f"Reading CSV file from: {csv_file_path}")
    csv_df = pd.read_csv(csv_file_path)

    # Check if ProbesNames column exists
    if "ProbesNames" not in csv_df.columns:
        print("Error: 'ProbesNames' column not found in CSV file")
        print(f"Available columns: {list(csv_df.columns)}")
        return

    # Step 3: Merge dataframes
    print("Merging data...")
    merged_df = csv_df.merge(
        blast_df, left_on="ProbesNames", right_on="ProbeName", how="left"
    )

    # Step 4: Filter to keep only unique hits
    print("Filtering for unique hits only...")
    cleaned_df = merged_df[merged_df["NumberOfHits"] == 1].copy()

    # Drop duplicate ProbeName column
    if "ProbeName" in cleaned_df.columns:
        cleaned_df.drop(columns=["ProbeName"], inplace=True)

    print(f"Original rows: {len(merged_df)}")
    print(f"Rows with unique hits: {len(cleaned_df)}")
    print(f"Rows removed: {len(merged_df) - len(cleaned_df)}")

    # AUTOMATED SAVING PROCESS

    # Save BLAST results dataframe: substitute .txt with .csv
    blast_csv_path = re.sub(r"\.txt$", ".csv", blast_file_path, flags=re.IGNORECASE)
    blast_df.to_csv(blast_csv_path, index=False)
    print(f"BLAST results dataframe automatically saved to: {blast_csv_path}")

    # Save cleaned CSV file: substitute .csv with -unique.csv
    base_csv_name = os.path.splitext(csv_file_path)[0]  # Remove .csv extension
    cleaned_csv_path = base_csv_name + "-unique.csv"
    cleaned_df.to_csv(cleaned_csv_path, index=False)
    print(f"Cleaned CSV automatically saved to: {cleaned_csv_path}")

    print("Analysis complete!")

    # Display summary
    print("\n=== SUMMARY ===")
    print(f"Total probes processed: {len(blast_df)}")
    print(f"Probes with unique hits: {len(blast_df[blast_df['NumberOfHits'] == 1])}")
    print(f"Probes with multiple hits: {len(blast_df[blast_df['NumberOfHits'] > 1])}")

    # Show detailed hit count distribution
    print("\n=== HIT COUNT DISTRIBUTION ===")
    hit_counts = blast_df["NumberOfHits"].value_counts().sort_index()
    for hits, count in hit_counts.items():
        print(f"Probes with {hits} hit(s): {count}")

    root.destroy()


if __name__ == "__main__":
    main()
