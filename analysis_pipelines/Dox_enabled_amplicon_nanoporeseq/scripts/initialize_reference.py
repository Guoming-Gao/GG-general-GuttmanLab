import os
import subprocess
import pysam
import argparse
import json
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path

# Lab standard paths (consistent with probe-designer config)
GENOME_FASTA = "/Volumes/guttman/genomes/mm10/fasta/mm10.fa"
VCF_FILE = "/Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
B6_SAMPLE = "C57BL_6NJ"
CAST_SAMPLE = "CAST_EiJ"
PRIMER_FILE = "ValidatedPrimers.fa"

def find_primer_genomic(primer_seq):
    """Finds the first genomic occurrence of a primer sequence in chrX mm10."""
    # We search the known Xist locus +/- 10kb
    start, end = 103460373 - 10000, 103483233 + 10000
    region = f"chrX:{start}-{end}"

    cmd = ["samtools", "faidx", GENOME_FASTA, region]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    seq = "".join(result.stdout.strip().split("\n")[1:]).upper()

    pos = seq.find(primer_seq)
    if pos != -1:
        return ("+", start + pos)

    pos_rc = seq.find(str(Seq(primer_seq).reverse_complement()))
    if pos_rc != -1:
        return ("-", start + pos_rc)

    return (None, None)

def get_vcf_snps(chrom, start, end):
    """Identifies all B6 vs Cast SNPs in a genomic region."""
    vcf = pysam.VariantFile(VCF_FILE)
    target_chrom = chrom.replace("chr", "") if "chr" in chrom else chrom
    if target_chrom not in vcf.header.contigs:
        target_chrom = "chr" + target_chrom

    snps = []
    for record in vcf.fetch(target_chrom, start - 1, end):
        b6_gt = record.samples[B6_SAMPLE].alleles
        cast_gt = record.samples[CAST_SAMPLE].alleles

        if b6_gt[0] is not None and cast_gt[0] is not None:
            b6_allele = "/".join(b6_gt)
            cast_allele = "/".join(cast_gt)
            if b6_allele != cast_allele:
                snps.append({
                    "genomic_pos": record.pos,
                    "ref": record.ref,
                    "alt": list(record.alts),
                    "b6": b6_allele,
                    "cast": cast_allele
                })
    return snps

def main():
    parser = argparse.ArgumentParser(description="Initialize reference sequence and SNPs from primers.")
    parser.add_argument("--output_dir", default=".", help="Project output directory (default: current directory)")
    parser.add_argument("--f_primer", default="AC_XistExAmp_5SNPs-F", help="Forward primer ID in ValidatedPrimers.fa")
    parser.add_argument("--r_primer", default="AC_XistExAmp_5SNPs-R", help="Reverse primer ID in ValidatedPrimers.fa")
    args = parser.parse_args()

    # 1. Load Primers
    if not os.path.exists(PRIMER_FILE):
        print(f"Error: {PRIMER_FILE} not found.")
        return

    primers = list(SeqIO.parse(PRIMER_FILE, "fasta"))
    try:
        f_seq = str(next(p for p in primers if p.id == args.f_primer).seq).upper()
        r_seq = str(next(p for p in primers if p.id == args.r_primer).seq).upper()
    except StopIteration:
        print("Error: Primers not found in file.")
        return

    print(f"Targeting: {args.f_primer} / {args.r_primer}")

    # 2. Map primers
    strand_f, pos_f = find_primer_genomic(f_seq)
    strand_r, pos_r = find_primer_genomic(r_seq)

    if not pos_f or not pos_r:
        print("Error: Could not map primers to Xist locus.")
        return

    # Amplicon range
    start_genomic = min(pos_f, pos_r)
    if pos_f > pos_r:
        end_genomic = pos_f + len(f_seq) - 1
    else:
        end_genomic = pos_r + len(r_seq) - 1

    print(f"Amplicon defined at chrX:{start_genomic}-{end_genomic}")

    # 3. Extract Reference
    cmd = ["samtools", "faidx", GENOME_FASTA, f"chrX:{start_genomic}-{end_genomic}"]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    ref_seq = "".join(result.stdout.strip().split("\n")[1:]).upper()

    # 4. Find SNPs
    snps = get_vcf_snps("chrX", start_genomic, end_genomic)
    print(f"Found {len(snps)} B6/Cast SNPs in amplicon.")

    # 5. Output
    results_ref_dir = os.path.join(args.output_dir, "results", "ref_seq")
    os.makedirs(results_ref_dir, exist_ok=True)

    fasta_path = os.path.join(results_ref_dir, "target_amplicon.fa")
    with open(fasta_path, "w") as f:
        f.write(f">Xist_Amplicon chrX:{start_genomic}-{end_genomic}\n")
        f.write(ref_seq + "\n")

    snp_path = os.path.join(results_ref_dir, "snps.json")
    with open(snp_path, "w") as f:
        for s in snps:
            s["local_pos"] = s["genomic_pos"] - start_genomic + 1
        json.dump(snps, f, indent=4)

    # 6. Align back to genome for validation (SAM)
    print("Aligning amplicon back to genome for validation...")
    sam_path = os.path.join(results_ref_dir, "amplicon_to_genome.sam")
    minimap_cmd = [
        "minimap2", "-ax", "sr",
        GENOME_FASTA, fasta_path
    ]
    with open(sam_path, "w") as f:
        subprocess.run(minimap_cmd, stdout=f, check=True)

    print(f"Reference, SNPs, and SAM validation created in {results_ref_dir}/")

if __name__ == "__main__":
    main()
