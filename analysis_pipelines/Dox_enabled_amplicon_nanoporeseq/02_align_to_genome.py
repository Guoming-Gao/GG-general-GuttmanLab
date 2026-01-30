
import os
import subprocess
import glob
import argparse

def run_alignment(fastq_dir, results_dir, genome_fasta, threads=8):
    align_dir = os.path.join(results_dir, "aligned")
    os.makedirs(align_dir, exist_ok=True)

    fastq_files = glob.glob(os.path.join(fastq_dir, "*.fastq"))

    for f in fastq_files:
        sample = os.path.basename(f).split(".")[0]
        out_bam = os.path.join(align_dir, f"{sample}.sorted.bam")

        if os.path.exists(out_bam):
            print(f"Skipping {sample}, output already exists.")
            continue

        print(f"Aligning {sample}...")
        temp_sam = out_bam + ".sam"

        # minimap2 -ax map-ont
        cmd_align = [
            "/opt/miniconda3/envs/bioinfo/bin/minimap2", "-ax", "map-ont", "-t", str(threads),
            genome_fasta, f
        ]

        with open(temp_sam, "w") as sam_file:
            subprocess.run(cmd_align, stdout=sam_file, check=True)

        # samtools sort
        print(f"Sorting {sample}...")
        subprocess.run(["/opt/miniconda3/envs/bioinfo/bin/samtools", "sort", "-@", str(threads), "-o", out_bam, temp_sam], check=True)

        # samtools index
        print(f"Indexing {sample}...")
        subprocess.run(["/opt/miniconda3/envs/bioinfo/bin/samtools", "index", out_bam], check=True)

        # cleanup
        os.remove(temp_sam)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", required=True)
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    run_alignment(args.data_dir, args.results_dir, args.genome, args.threads)
