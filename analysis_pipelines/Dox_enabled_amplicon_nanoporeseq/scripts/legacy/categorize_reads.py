import os
import subprocess
import argparse
import json
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def run_cutadapt_categorize(fastq_in, fwd_primer, rev_primer_rc, output_dir, sample_name):
    """
    Uses cutadapt to count Head (Fwd), Tail (RC of 2PBC), and Both.
    Logic: We run cutadapt in several passes or use its internal logic to flag reads.
    For categorization, we'll use cutadapt's --info-file or simply run specific filters.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Paths for temporary/intermediate files
    both_fastq = os.path.join(output_dir, f"{sample_name}.both.fastq.gz")
    info_file = os.path.join(output_dir, f"{sample_name}.cutadapt_info.txt")

    # Linked adapter PASS (enforces order: FWD ... REV_RC)
    # Using --revcomp to handle mixed orientation
    # We use --action=none because the user wants ALL reads for global alignment later,
    # but we need the labels for categorization.

    linked_adapter = f"{fwd_primer}...{rev_primer_rc}"

    cmd = [
        "cutadapt",
        "-g", linked_adapter,
        "--revcomp",
        "-e", "0.2",
        "--overlap", "15",
        "--info-file", info_file,
        "-o", "/dev/null", # We don't need the trimmed output here, just the info
        fastq_in
    ]

    print(f"[{sample_name}] Categorizing reads with Cutadapt...")
    subprocess.run(["conda", "run", "-n", "bioinfo"] + cmd, check=True, capture_output=True)

    # Parse info file to categorize
    # Columns in info-file:
    # 1: read name
    # 2: number of adapters found
    # 3: start of first match...

    stats = {"Both": 0, "Head-only": 0, "Tail-only": 0, "Neither": 0, "Total": 0}
    read_to_cat = {}

    with open(info_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2: continue
            read_id = parts[0]
            num_adapters = int(parts[1])
            stats["Total"] += 1

            if num_adapters == 0:
                cat = "Neither"
            elif num_adapters == 1:
                # Need to check which one matched if possible from info file,
                # but linked adapter -g A...B usually only counts as 1 "match" if both found?
                # Actually cutadapt -g A...B reports 1 if the linked pair is found.
                # To get partials, we would need separate passes.
                cat = "Both"
            elif num_adapters == -1: # Info file placeholder
                cat = "Neither"
            else:
                cat = "Both"

            stats[cat] += 1
            read_to_cat[read_id] = cat

    # Since the linked pass only tells us "Both", let's do independent passes for Head and Tail
    # to find true "Head-only" and "Tail-only".

    # Pass 1: Head search
    cmd_head = ["cutadapt", "-g", fwd_primer, "--revcomp", "-e", "0.2", "--overlap", "15", "--info-file", f"{info_file}.head", "-o", "/dev/null", fastq_in]
    subprocess.run(["conda", "run", "-n", "bioinfo"] + cmd_head, check=True, capture_output=True)

    # Pass 2: Tail search
    cmd_tail = ["cutadapt", "-a", rev_primer_rc, "--revcomp", "-e", "0.2", "--overlap", "15", "--info-file", f"{info_file}.tail", "-o", "/dev/null", fastq_in]
    subprocess.run(["conda", "run", "-n", "bioinfo"] + cmd_tail, check=True, capture_output=True)

    head_reads = set()
    with open(f"{info_file}.head", "r") as f:
        for line in f:
            p = line.split("\t")
            if len(p) > 1 and int(p[1]) > 0: head_reads.add(p[0])

    tail_reads = set()
    with open(f"{info_file}.tail", "r") as f:
        for line in f:
            p = line.split("\t")
            if len(p) > 1 and int(p[1]) > 0: tail_reads.add(p[0])

    # Final Categorization
    final_read_cat = []
    stats = {"Both": 0, "Head-only": 0, "Tail-only": 0, "Neither": 0, "Total": 0}

    # Get all read IDs from one of the info files (they should be the same)
    all_reads = []
    with open(f"{info_file}.head", "r") as f:
        for line in f:
            all_reads.append(line.split("\t")[0])

    for rid in all_reads:
        has_h = rid in head_reads
        has_t = rid in tail_reads
        stats["Total"] += 1
        if has_h and has_t:
            cat = "Both"
        elif has_h:
            cat = "Head-only"
        elif has_t:
            cat = "Tail-only"
        else:
            cat = "Neither"

        stats[cat] += 1
        final_read_cat.append({"ReadID": rid, "Category": cat})

    # Save mapping
    df_cat = pd.DataFrame(final_read_cat)
    cat_mapping_path = os.path.join(output_dir, f"{sample_name}.categorization.csv")
    df_cat.to_csv(cat_mapping_path, index=False)

    return stats

def main():
    parser = argparse.ArgumentParser(description="Categorize reads into Head/Tail/Both/None buckets.")
    parser.add_argument("--data_dir", required=True)
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--f_primer", required=True)
    parser.add_argument("--r_primer", required=True)
    args = parser.parse_args()

    cat_dir = os.path.join(args.results_dir, "categorized")
    os.makedirs(cat_dir, exist_ok=True)

    fastq_files = glob.glob(os.path.join(args.data_dir, "*.fastq")) + glob.glob(os.path.join(args.data_dir, "*.fastq.gz"))

    fwd = args.f_primer
    # The tail we look for is the RC of the 2PBC sequence in the R_primer
    rev_rc = str(Seq(args.r_primer).reverse_complement())

    all_stats = []
    for f in fastq_files:
        sample = os.path.basename(f).split(".")[0]
        stats = run_cutadapt_categorize(f, fwd, rev_rc, cat_dir, sample)
        stats["Sample"] = sample
        all_stats.append(stats)

    df_summary = pd.DataFrame(all_stats)
    summary_path = os.path.join(cat_dir, "categorization_summary.csv")
    df_summary.to_csv(summary_path, index=False)

    print("\nCategorization Summary:")
    print(df_summary.to_string())

if __name__ == "__main__":
    main()
