
import os
import glob
import subprocess

DATA_DIR = "data"
REF_FILE = "ref_seq/xist_amplicon_ref.fa"
OUTPUT_DIR = "aligned"
FASTQ_FILES = glob.glob(os.path.join(DATA_DIR, "*.fastq"))

def run_command(cmd):
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
    return result

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    for fastq in FASTQ_FILES:
        base_name = os.path.basename(fastq).replace(".fastq", "")
        bam_output = os.path.join(OUTPUT_DIR, f"{base_name}.bam")
        sorted_bam = os.path.join(OUTPUT_DIR, f"{base_name}.sorted.bam")

        # Alignment
        # -ax map-ont for Nanopore
        # -t 4 for threads
        minimap_cmd = [
            "minimap2", "-ax", "map-ont", "-t", "4",
            REF_FILE, fastq
        ]

        # We'll pipe to samtools view to get BAM directly if possible,
        # but let's just do it in steps for robustness
        sam_output = os.path.join(OUTPUT_DIR, f"{base_name}.sam")
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
