import pandas as pd
import openpyxl
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import os


def annotate_plate_wellpos(csv_df, xlsx_path):
    """
    Annotate CSV dataframe with PlateName, WellPos, and IDTseq from Excel file
    """
    # Load the Excel file
    xlsx = pd.ExcelFile(xlsx_path)

    # Prepare columns for annotation
    plate_names = []
    well_positions = []
    idt_seqs = []

    # For efficient searching, load all sheets into dict of DataFrames
    sheets_dict = {}
    for sheet_name in xlsx.sheet_names:
        try:
            sheets_dict[sheet_name] = xlsx.parse(sheet_name)
        except:
            print(f"Warning: Could not load sheet '{sheet_name}', skipping...")
            continue

    print(f"Loaded {len(sheets_dict)} sheets from Excel file")

    # For each Seq in csv, find matching sheet and well position
    for i, seq in enumerate(csv_df["Seq"]):
        found_plate = None
        found_well = None
        found_idtseq = None

        # Search in each sheet
        for sheet_name, df_sheet in sheets_dict.items():
            # FIXED: Check for different possible well position column names
            well_pos_col = None
            for col in df_sheet.columns:
                if "well" in col.lower() and "position" in col.lower():
                    well_pos_col = col
                    break

            # Check if required columns exist
            if "Sequence" in df_sheet.columns and well_pos_col is not None:
                # Find rows where Sequence contains seq (case-insensitive)
                try:
                    match_rows = df_sheet[
                        df_sheet["Sequence"].str.contains(seq, na=False, case=False)
                    ]
                    if not match_rows.empty:
                        # Take the first match
                        found_plate = sheet_name
                        found_well = match_rows.iloc[0][well_pos_col]
                        found_idtseq = match_rows.iloc[0]["Sequence"]
                        break
                except:
                    continue

        plate_names.append(found_plate)
        well_positions.append(found_well)
        idt_seqs.append(found_idtseq)

        # Progress indicator
        if (i + 1) % 10 == 0:
            print(f"Processed {i + 1}/{len(csv_df)} sequences...")

    # Create a copy of the dataframe to avoid modifying the original
    result_df = csv_df.copy()

    # FIXED: Insert all new columns at the beginning
    result_df.insert(0, "IDTseq", idt_seqs)
    result_df.insert(0, "WellPos", well_positions)
    result_df.insert(0, "PlateName", plate_names)

    return result_df


def add_odd_even_column(csv_df):
    """
    Add OddEven column based on theStartPos sorting
    """
    # Sort by theStartPos first
    csv_df_sorted = csv_df.sort_values(by="theStartPos").reset_index(drop=True)

    # Label Odd or Even based on position in sorted order (1-based indexing)
    # 1st position = Odd, 2nd position = Even, 3rd position = Odd, etc.
    csv_df_sorted["OddEven"] = [
        "Odd" if (i % 2) == 0 else "Even" for i in range(len(csv_df_sorted))
    ]

    # FIXED: Move OddEven column to the beginning (after PlateName, WellPos, IDTseq)
    oddeven_col = csv_df_sorted.pop("OddEven")
    csv_df_sorted.insert(3, "OddEven", oddeven_col)

    return csv_df_sorted


def main():
    # Create GUI window (hidden)
    root = Tk()
    root.withdraw()

    print("Starting plate and well position annotation...")

    # Step 1: Select CSV file
    print("Please select your CSV file to annotate...")
    csv_path = askopenfilename(
        title="Select CSV file to annotate",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
    )

    if not csv_path:
        print("No CSV file selected. Exiting.")
        return

    # Step 2: Select Excel file
    print("Please select your Excel file with plate data...")
    xlsx_path = askopenfilename(
        title="Select Excel file with plate data",
        filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")],
    )

    if not xlsx_path:
        print("No Excel file selected. Exiting.")
        return

    # Step 3: Read CSV file
    print(f"Reading CSV file: {csv_path}")
    csv_df = pd.read_csv(csv_path)
    print(f"Loaded {len(csv_df)} rows from CSV file")

    # Step 4: Annotate with PlateName, WellPos, and IDTseq
    print("Annotating with plate and well position information...")
    csv_df = annotate_plate_wellpos(csv_df, xlsx_path)

    # Step 5: Add OddEven column
    print("Adding OddEven column based on theStartPos...")
    csv_df = add_odd_even_column(csv_df)

    # Step 6: Sort by OddEven (Odd first), then PlateName, then WellPos
    print("Sorting table by OddEven (Odd first), PlateName, and WellPos...")
    # Custom sort order for OddEven - Odd first, then Even
    csv_df["OddEven"] = pd.Categorical(
        csv_df["OddEven"], categories=["Odd", "Even"], ordered=True
    )

    csv_df = csv_df.sort_values(
        by=["OddEven", "PlateName", "WellPos"], na_position="last"
    ).reset_index(drop=True)

    # Step 7: Save as Excel file with updated postfix
    base_name = os.path.splitext(csv_path)[0]
    save_path = base_name + "-withPlateWellPos.xlsx"

    print(f"Saving annotated file as: {save_path}")
    csv_df.to_excel(save_path, index=False)

    # Display summary
    print("\n=== ANNOTATION SUMMARY ===")
    print(f"Total rows processed: {len(csv_df)}")
    print(f"Rows with plate matches: {len(csv_df.dropna(subset=['PlateName']))}")
    print(f"Rows without matches: {len(csv_df[csv_df['PlateName'].isna()])}")

    # Show plate distribution
    if not csv_df["PlateName"].isna().all():
        print("\n=== PLATE DISTRIBUTION ===")
        plate_counts = csv_df["PlateName"].value_counts()
        for plate, count in plate_counts.items():
            print(f"{plate}: {count} probes")

    # Show OddEven distribution
    print("\n=== ODD/EVEN DISTRIBUTION ===")
    oddeven_counts = csv_df["OddEven"].value_counts()
    for category, count in oddeven_counts.items():
        print(f"{category}: {count} probes")

    # Show specific check for C9orf72-intron1-123
    probe_123 = csv_df[csv_df["ProbesNames"] == "C9orf72-intron1-123"]
    if not probe_123.empty:
        print(f"\n=== C9orf72-intron1-123 ANNOTATION ===")
        print(f"Plate: {probe_123.iloc[0]['PlateName']}")
        print(f"Well: {probe_123.iloc[0]['WellPos']}")
        print(f"IDT Sequence: {probe_123.iloc[0]['IDTseq'][:60]}...")

    print(f"\nAnnotation complete! File saved as: {save_path}")

    root.destroy()


if __name__ == "__main__":
    main()
