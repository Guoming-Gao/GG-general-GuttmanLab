
import os
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess

# Paths from config
GENOM_FASTA = "/Volumes/guttman/genomes/mm10/fasta/mm10.fa"
VCF_FILE = "/Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
GB_FILE = "ref_seq/xist-ref-seq-with-snps-annot-1712-2295.gb"

# Allele samples
B6_SAMPLE = "C57BL_6NJ"
CAST_SAMPLE = "CAST_EiJ"

def get_vcf_alleles(chrom, pos, ref, vcf_path, b6_sample, cast_sample):
    vcf = pysam.VariantFile(vcf_path)

    # Try both 'chrX' and 'X'
    contigs = list(vcf.header.contigs)
    if chrom not in contigs:
        alt_chrom = chrom.replace("chr", "")
        if alt_chrom in contigs:
            chrom = alt_chrom
        else:
            print(f"Warning: Contig {chrom} not found in VCF header.")
            return None

    region = f"{chrom}:{pos}-{pos}"
    alleles = {"B6": None, "Cast": None, "Ref": None, "Alt": None}

    try:
        # Use fetch with explicit coordinates
        for record in vcf.fetch(chrom, pos-1, pos):
            if record.pos == pos:
                alleles["Ref"] = record.ref
                alleles["Alt"] = list(record.alts)

                b6_gt = record.samples[b6_sample].alleles
                cast_gt = record.samples[cast_sample].alleles

                alleles["B6"] = "/".join(b6_gt) if b6_gt[0] is not None else "N/A"
                alleles["Cast"] = "/".join(cast_gt) if cast_gt[0] is not None else "N/A"
                return alleles
    except Exception as e:
        print(f"Error fetching VCF for {region}: {e}")
    return None

def main():
    # 1. Load GenBank sequence
    gb_record = SeqIO.read(GB_FILE, "genbank")
    gb_seq = str(gb_record.seq).upper()
    print(f"Loaded GenBank sequence: {len(gb_seq)} bp")

    snps_in_gb = []
    for feature in gb_record.features:
        if feature.type == "misc_feature":
            label = feature.qualifiers.get("label", ["Unknown"])[0]
            pos = int(feature.location.start) + 1
            snps_in_gb.append({"gb_pos": pos, "label": label})

    print(f"Found {len(snps_in_gb)} SNPs annotated in GenBank.")

    xist_chrom = "chrX"
    xist_start = 103460373
    xist_end = 103483233

    cmd = ["samtools", "faidx", GENOM_FASTA, f"{xist_chrom}:{xist_start}-{xist_end}"]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    genomic_xist = "".join(result.stdout.strip().split("\n")[1:]).upper()

    rna_sense_xist = str(Seq(genomic_xist).reverse_complement())

    match_index = rna_sense_xist.find(gb_seq)
    found_strand = None
    if match_index != -1:
        print(f"Found match in Xist RNA sense sequence at relative position {match_index}")
        genomic_5p = xist_end - match_index
        found_strand = "-"
    else:
        match_index_plus = genomic_xist.find(gb_seq)
        if match_index_plus != -1:
            print(f"Found match in Xist genomic (+) sequence at relative position {match_index_plus}")
            genomic_5p = xist_start + match_index_plus
            found_strand = "+"
        else:
            print("ERROR: Could not find GenBank sequence in Xist locus!")
            return

    print(f"Genomic range (1-indexed): {xist_chrom}:{genomic_5p} ({found_strand} strand)")

    print("\nSNP Verification Report:")
    print("-" * 80)
    print(f"{'GB Pos':<8} {'Genomic':<12} {'Label':<10} {'B6':<10} {'Cast':<10} {'Ref':<5} {'Alt':<10}")
    for snp in snps_in_gb:
        if found_strand == "-":
            genomic_snp_pos = genomic_5p - (snp["gb_pos"] - 1)
        else:
            genomic_snp_pos = genomic_5p + (snp["gb_pos"] - 1)

        vcf_info = get_vcf_alleles(xist_chrom, genomic_snp_pos, "", VCF_FILE, B6_SAMPLE, CAST_SAMPLE)

        if vcf_info:
            print(f"{snp['gb_pos']:<8} {genomic_snp_pos:<12} {snp['label']:<10} {vcf_info['B6']:<10} {vcf_info['Cast']:<10} {vcf_info['Ref']:<5} {str(vcf_info['Alt']):<10}")
        else:
            print(f"{snp['gb_pos']:<8} {genomic_snp_pos:<12} {snp['label']:<10} {'NOT FOUND':<10}")

if __name__ == "__main__":
    main()
