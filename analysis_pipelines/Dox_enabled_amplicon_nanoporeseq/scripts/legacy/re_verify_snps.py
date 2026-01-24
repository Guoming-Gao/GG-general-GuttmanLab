
import pysam
import os
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess

# Paths
GENOM_FASTA = "/Volumes/guttman/genomes/mm10/fasta/mm10.fa"
VCF_FILE = "/Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
GB_FILE = "/Users/gmgao/GGscripts/GG-general-GuttmanLab/analysis_pipelines/Dox_enabled_amplicon_nanoporeseq/ref_seq/xist-ref-seq-with-snps-annot-1712-2295.gb"
B6_SAMPLE = "C57BL_6NJ"
CAST_SAMPLE = "CAST_EiJ"

def get_snps_in_region(chrom, start, end, vcf_path, b6_sample, cast_sample):
    vcf = pysam.VariantFile(vcf_path)
    # Check contig naming
    contigs = list(vcf.header.contigs)
    if chrom not in contigs:
        alt_chrom = chrom.replace("chr", "")
        if alt_chrom in contigs:
            chrom = alt_chrom

    snps = []
    try:
        for record in vcf.fetch(chrom, start, end):
            # record.pos is 1-based
            b6_gt = record.samples[b6_sample].alleles
            cast_gt = record.samples[cast_sample].alleles

            # Skip if missing or same
            if b6_gt[0] is None or cast_gt[0] is None:
                continue

            b6_str = "/".join(b6_gt)
            cast_str = "/".join(cast_gt)

            if b6_str != cast_str:
                snps.append({
                    "pos": record.pos,
                    "ref": record.ref,
                    "alt": list(record.alts),
                    "b6": b6_str,
                    "cast": cast_str
                })
    except Exception as e:
        print(f"Error fetching VCF: {e}")
    return snps

def main():
    # Genomic match for the amplicon from previous run: chrX:103466603-103467186
    chrom = "chrX"
    start = 103466603
    end = 103467186

    print(f"Scanning range {chrom}:{start}-{end} for ANY SNPs where B6 != Cast...")
    vcf_snps = get_snps_in_region(chrom, start, end, VCF_FILE, B6_SAMPLE, CAST_SAMPLE)

    print(f"\nFound {len(vcf_snps)} SNPs in VCF where B6 != Cast:")
    print("-" * 80)
    print(f"{'Genomic Pos':<12} {'Ref':<5} {'Alt':<10} {'B6 Allele':<10} {'Cast Allele':<10}")
    for s in vcf_snps:
        print(f"{s['pos']:<12} {s['ref']:<5} {str(s['alt']):<10} {s['b6']:<10} {s['cast']:<10}")

    # Load GenBank to compare
    gb_record = SeqIO.read(GB_FILE, "genbank")

    # Map GenBank pos to Genomic pos (since it was identified as PLUS strand mapping)
    print("\nMapping Back to GenBank Positions:")
    print("-" * 80)
    print(f"{'Genomic Pos':<12} {'GB Pos':<10} {'B6':<5} {'Cast':<5}")
    for s in vcf_snps:
        gb_pos = s['pos'] - start + 1
        print(f"{s['pos']:<12} {gb_pos:<10} {s['b6']:<5} {s['cast']:<5}")

if __name__ == "__main__":
    main()
