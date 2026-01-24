import os
import glob
import subprocess
import argparse

def run_command(cmd):
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
    return result

def main():
    parser = argparse.ArgumentParser(description="Align Nanopore reads to amplicon reference.")
    parser.add_argument("--data_dir", default="data", help="Directory containing input FASTQ files (default: data)")
    parser.add_argument("--results_dir", default=None, help="Results directory (default: ./results)")
    parser.add_argument("--ref_file", default=None, help="Path to reference FASTA (default: results/ref_seq/target_amplicon.fa)")
    args = parser.parse_args()

    results_dir = args.results_dir if args.results_dir else os.path.join(".", "results")
    align_dir = os.path.join(results_dir, "aligned")
    ref_file = args.ref_file if args.ref_file else os.path.join(results_dir, "ref_seq", "target_amplicon.fa")

    fastq_files = glob.glob(os.path.join(args.data_dir, "*.fastq"))

    if not fastq_files:
        print(f"No FASTQ files found in {args.data_dir}")
        return

    if not os.path.exists(align_dir):
        os.makedirs(align_dir, exist_ok=True)

    for fastq in fastq_files:
        base_name = os.path.basename(fastq).replace(".fastq", "")
        bam_output = os.path.join(align_dir, f"{base_name}.bam")
        sorted_bam = os.path.join(align_dir, f"{base_name}.sorted.bam")

        # Alignment
        # -ax map-ont for Nanopore
        # -t 4 for threads
        minimap_cmd = [
            "minimap2", "-ax", "map-ont", "-t", "4",
            ref_file, fastq
        ]

        sam_output = os.path.join(align_dir, f"{base_name}.sam")
        with open(sam_output, "w") as f:
            subprocess.run(minimap_cmd, stdout=f, check=True)

        # Convert to BAM, Sort, and Index
        run_command(["samtools", "view", "-bS", sam_output, "-o", bam_output])
        run_command(["samtools", "sort", bam_output, "-o", sorted_bam])
        run_command(["samtools", "index", sorted_bam])

        # Clean up intermediate files
        os.remove(sam_output)
        os.remove(bam_output)

        print(f"Finished {base_name}")

if __name__ == "__main__":
    main()
