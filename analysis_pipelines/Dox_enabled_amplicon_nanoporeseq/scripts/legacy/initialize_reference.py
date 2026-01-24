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

def find_primer_genomic(primer_seq, genome_fasta, region="chrX"):
    """Finds the first genomic occurrence of a primer sequence in a region."""
    # Search Xist locus +/- 20kb
    start, end = 103460373 - 20000, 103483233 + 20000
    region_str = f"{region}:{start}-{end}"

    cmd = ["samtools", "faidx", genome_fasta, region_str]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    seq = "".join(result.stdout.strip().split("\n")[1:]).upper()

    pos = seq.find(primer_seq.upper())
    if pos != -1:
        return ("+", start + pos)

    pos_rc = seq.find(str(Seq(primer_seq).reverse_complement()).upper())
    if pos_rc != -1:
        return ("-", start + pos_rc)

    return (None, None)

def get_vcf_snps(vcf_file, b6_sample, cast_sample, chrom, start, end):
    """Identifies all B6 vs Cast SNPs in a genomic region."""
    vcf = pysam.VariantFile(vcf_file)
    target_chrom = chrom.replace("chr", "") if "chr" in chrom else chrom
    if target_chrom not in vcf.header.contigs:
        target_chrom = "chr" + target_chrom

    snps = []
    for record in vcf.fetch(target_chrom, start - 1, end):
        try:
            b6_gt = record.samples[b6_sample].alleles
            cast_gt = record.samples[cast_sample].alleles
        except KeyError:
            continue

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
    parser = argparse.ArgumentParser(description="Initialize genomic reference info and SNPs.")
    parser.add_argument("--results_dir", required=True, help="Results directory")
    parser.add_argument("--f_primer", required=True, help="Forward primer sequence")
    parser.add_argument("--r_primer", required=True, help="Reverse primer sequence (containing 2PBC)")
    parser.add_argument("--genome_fasta", default="/Volumes/guttman/genomes/mm10/fasta/mm10.fa")
    parser.add_argument("--vcf_file", default="/Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz")
    parser.add_argument("--b6_sample", default="C57BL_6NJ")
    parser.add_argument("--cast_sample", default="CAST_EiJ")
    args = parser.parse_args()

    results_ref_dir = os.path.join(args.results_dir, "ref_seq")
    os.makedirs(results_ref_dir, exist_ok=True)

    print(f"Mapping primers to mm10...")
    strand_f, pos_f = find_primer_genomic(args.f_primer, args.genome_fasta)
    strand_r, pos_r = find_primer_genomic(args.r_primer, args.genome_fasta)

    # Note: R_primer might be a random primer, let's try just the constant part if it's the RT primer
    if not pos_r and "N" in args.r_primer.upper():
        constant_part = args.r_primer.upper().split("N")[0]
        if len(constant_part) > 15:
            print(f"R_primer mapping failed with Ns, trying constant part: {constant_part}")
            strand_r, pos_r = find_primer_genomic(constant_part, args.genome_fasta)

    if not pos_f or not pos_r:
        print(f"Warning: Could not map one or both primers exactly (F:{pos_f}, R:{pos_r}).")
        # For Xist, we know the approximate locus, we can proceed if we have at least one anchor or just use the locus
        start_genomic, end_genomic = 103460373, 103483233
    else:
        start_genomic = min(pos_f, pos_r)
        end_genomic = max(pos_f + len(args.f_primer), pos_r + len(args.r_primer))
        print(f"Amplicon defined at chrX:{start_genomic}-{end_genomic}")

    # 3. Fetch SNPs
    print(f"Fetching SNPs for chrX:{start_genomic}-{end_genomic}...")
    snps = get_vcf_snps(args.vcf_file, args.b6_sample, args.cast_sample, "chrX", start_genomic, end_genomic)
    print(f"Found {len(snps)} B6/Cast SNPs.")

    # 4. Output
    snp_path = os.path.join(results_ref_dir, "snps.json")
    with open(snp_path, "w") as f:
        json.dump(snps, f, indent=4)

    # Save a small reference fasta for the Xist locus for quick visualization
    fasta_path = os.path.join(results_ref_dir, "xist_locus.fa")
    cmd = ["samtools", "faidx", args.genome_fasta, f"chrX:{start_genomic}-{end_genomic}"]
    with open(fasta_path, "w") as f:
        subprocess.run(cmd, stdout=f, check=True)

    print(f"Step 1 Complete. SNP data saved to {snp_path}")

if __name__ == "__main__":
    main()
