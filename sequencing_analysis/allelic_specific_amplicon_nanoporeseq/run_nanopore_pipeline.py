#!/usr/bin/env python3
"""Run the Nanopore amplicon pipeline with optional inline barcode demux."""

import argparse
import subprocess
from pathlib import Path


def run(cmd):
    print(" ".join(str(part) for part in cmd))
    subprocess.run([str(part) for part in cmd], check=True)


def has_blank_sample_names(sample_map):
    with open(sample_map) as handle:
        header = handle.readline().rstrip("\n").split(",")
        try:
            sample_idx = header.index("SampleName")
        except ValueError as exc:
            raise ValueError("sample map must include a SampleName column") from exc
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split(",")
            if sample_idx >= len(fields) or not fields[sample_idx].strip():
                return True
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, help="Dataset directory")
    parser.add_argument("--fastq", help="Single multiplexed FASTQ/FASTQ.GZ. Required with --sample-map.")
    parser.add_argument("--sample-map", help="Sample map for inline barcode demultiplexing")
    parser.add_argument("--fastq-dir", help="Already-demultiplexed FASTQ directory. Defaults to data-dir or demux output.")
    parser.add_argument("--results-dir", help="Results directory. Default: data-dir/results/mode_1")
    parser.add_argument("--genome", required=True, help="Genome FASTA for minimap2")
    parser.add_argument("--vcf", required=True, help="Indexed VCF for SNP extraction")
    parser.add_argument("--region", default="chrX:103460373-103483233", help="Target genomic region")
    parser.add_argument("--b6", default="C57BL_6NJ", help="B6 sample name in VCF")
    parser.add_argument("--cast", default="CAST_EiJ", help="Cast sample name in VCF")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--python", default="python3", help="Python executable for step scripts")
    parser.add_argument("--max-mismatches", type=int, default=1, help="Barcode mismatch tolerance for demux")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    data_dir = Path(args.data_dir)
    results_dir = Path(args.results_dir) if args.results_dir else data_dir / "results" / "mode_1"
    fastq_dir = Path(args.fastq_dir) if args.fastq_dir else data_dir

    if args.sample_map:
        if not args.fastq:
            raise ValueError("--fastq is required when --sample-map is provided")
        if has_blank_sample_names(args.sample_map):
            raise ValueError("Sample map has blank SampleName values")
        fastq_dir = data_dir / "fastq_by_sample"
        audit_dir = data_dir / "demux_audit"
        run([
            args.python,
            script_dir / "00_demultiplex_inline_barcodes.py",
            "--fastq", args.fastq,
            "--sample-map", args.sample_map,
            "--output-dir", fastq_dir,
            "--summary", data_dir / "demux_summary.csv",
            "--unassigned-fastq", audit_dir / "unassigned.fastq",
            "--max-mismatches", args.max_mismatches,
        ])

    run([args.python, script_dir / "01_fastq_quality_metrics.py", "--data_dir", fastq_dir, "--results_dir", results_dir])
    run([
        args.python,
        script_dir / "02_align_to_genome.py",
        "--data_dir", fastq_dir,
        "--results_dir", results_dir,
        "--genome", args.genome,
        "--threads", args.threads,
    ])
    run([
        args.python,
        script_dir / "03_extract_snps_from_vcf.py",
        "--vcf", args.vcf,
        "--region", args.region,
        "--b6", args.b6,
        "--cast", args.cast,
        "--output", results_dir / "xist_snps.txt",
    ])
    run([
        args.python,
        script_dir / "04_quantify_alleles.py",
        "--bam_dir", results_dir / "aligned",
        "--snp_file", results_dir / "xist_snps.txt",
        "--output_dir", results_dir / "quantification",
        "--region", args.region,
    ])
    run([
        args.python,
        script_dir / "05_plot_allele_summary.py",
        "--results-dir", results_dir,
        "--sample-map", args.sample_map or "",
    ])


if __name__ == "__main__":
    main()
