import os
import subprocess
import argparse
import glob

def align_minimap2(fastq, genome_fasta, output_bam, threads=4):
    """Aligns FASTQ to genome using minimap2 -ax map-ont."""
    print(f"Aligning {os.path.basename(fastq)} to {os.path.basename(genome_fasta)}...")

    # Run minimap2 and pipe to samtools for sorting and indexing
    cmd = [
        "minimap2", "-ax", "map-ont", "-t", str(threads),
        genome_fasta, fastq
    ]

    # We pipe minimap2 output directly to samtools sort
    with open(output_bam.replace(".sorted.bam", ".sam"), "w") as f_sam:
        subprocess.run(cmd, stdout=f_sam, check=True)

    sam_file = output_bam.replace(".sorted.bam", ".sam")

    # Convert to BAM, sort, and index
    sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", output_bam, sam_file]
    subprocess.run(sort_cmd, check=True)

    index_cmd = ["samtools", "index", output_bam]
    subprocess.run(index_cmd, check=True)

    # Cleanup SAM
    if os.path.exists(sam_file):
        os.remove(sam_file)

def main():
    parser = argparse.ArgumentParser(description="Global alignment of reads to mm10.")
    parser.add_argument("--data_dir", required=True)
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--genome_fasta", default="/Volumes/guttman/genomes/mm10/fasta/mm10.fa")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    align_dir = os.path.join(args.results_dir, "aligned")
    os.makedirs(align_dir, exist_ok=True)

    fastq_files = glob.glob(os.path.join(args.data_dir, "*.fastq")) + glob.glob(os.path.join(args.data_dir, "*.fastq.gz"))

    for f in fastq_files:
        sample = os.path.basename(f).split(".")[0]
        output_bam = os.path.join(align_dir, f"{sample}.genomic.sorted.bam")

        if os.path.exists(output_bam):
            print(f"Skipping {sample}, output already exists.")
            continue

        align_minimap2(f, args.genome_fasta, output_bam, args.threads)

    print("\nGlobal Alignment Complete.")

if __name__ == "__main__":
    main()
