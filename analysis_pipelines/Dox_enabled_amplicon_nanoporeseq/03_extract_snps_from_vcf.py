
import pysam
import os
import argparse

def extract_snps(vcf_path, region, b6_sample, cast_sample, output_path):
    vcf = pysam.VariantFile(vcf_path)
    samples = [b6_sample, cast_sample]
    for s in samples:
        if s not in vcf.header.samples:
            print(f"Error: Sample {s} not in VCF header.")
            return

    vcf.subset_samples(samples)

    chrom = region.split(":")[0]
    start = int(region.split(":")[1].split("-")[0])
    end = int(region.split("-")[1])

    try:
        records = list(vcf.fetch(chrom, start, end))
    except ValueError:
        chrom = chrom.replace("chr", "")
        records = list(vcf.fetch(chrom, start, end))

    with open(output_path, "w") as out:
        out.write("POS\tREF\tALT\tB6_GT\tCAST_GT\tB6_BASE\tCAST_BASE\n")
        for record in records:
            b6_gt = record.samples[b6_sample]['GT']
            cast_gt = record.samples[cast_sample]['GT']

            if b6_gt != cast_gt:
                b6_base = record.ref if b6_gt == (0, 0) else record.alleles[b6_gt[0]] if b6_gt[0] is not None else "N"
                cast_base = record.ref if cast_gt == (0, 0) else record.alleles[cast_gt[0]] if cast_gt[0] is not None else "N"

                if b6_base != cast_base:
                    alts = ",".join(record.alts) if record.alts else "."
                    out.write(f"{record.pos}\t{record.ref}\t{alts}\t{record.samples[b6_sample]['GT']}\t{record.samples[cast_sample]['GT']}\t{b6_base}\t{cast_base}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--b6", required=True)
    parser.add_argument("--cast", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    extract_snps(args.vcf, args.region, args.b6, args.cast, args.output)
