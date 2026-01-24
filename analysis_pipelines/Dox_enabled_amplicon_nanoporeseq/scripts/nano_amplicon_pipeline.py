
import os
import sys
import subprocess
import argparse
import json
import glob
import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pysam
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

# --- Configuration & Constants ---
XIST_CHROM = "chrX"
XIST_START = 103460373
XIST_END = 103483233
DEFAULT_GENOME = "/Volumes/guttman/genomes/mm10/fasta/mm10.fa"
DEFAULT_GTF = "/Volumes/guttman/genomes/mm10/annotation/mm10.refGene.gtf.gz"
BIN_SIZE = 100000  # 100kb bins for coarse Gene mapping

class NanoAmpliconPipeline:
    def __init__(self, data_dir, results_dir, f_primer, r_primer, genome_fasta=DEFAULT_GENOME, gtf_path=DEFAULT_GTF):
        self.data_dir = data_dir
        self.results_dir = results_dir
        self.f_primer = f_primer
        self.r_primer = r_primer
        self.genome_fasta = genome_fasta
        self.gtf_path = gtf_path

        self.cat_dir = os.path.join(results_dir, "categorized")
        self.align_dir = os.path.join(results_dir, "aligned")
        self.spec_dir = os.path.join(results_dir, "specificity")
        self.graph_dir = os.path.join(results_dir, "graphs")
        self.report_dir = os.path.join(results_dir, "reports")

        for d in [self.cat_dir, self.align_dir, self.spec_dir, self.graph_dir, self.report_dir]:
            os.makedirs(d, exist_ok=True)

    # --- Step 1: Categorization ---
    def categorize_reads(self, threads=4):
        print("\n>>> Phase 1: Read Categorization (Cutadapt)")
        fastq_files = glob.glob(os.path.join(self.data_dir, "*.fastq*"))
        fwd_primer = self.f_primer.upper()
        rev_rc = str(Seq(self.r_primer).reverse_complement()).upper()

        summary_stats = []

        for f in fastq_files:
            sample = os.path.basename(f).split(".")[0]
            print(f" Processing {sample}...")

            # Linked search
            info_linked = os.path.join(self.cat_dir, f"{sample}.linked.info.txt")
            cmd_linked = ["cutadapt", "-g", f"{fwd_primer}...{rev_rc}", "--revcomp", "-e", "0.25", "--overlap", "10", "--info-file", info_linked, "-o", "/dev/null", f]
            subprocess.run(cmd_linked, check=True, capture_output=True)

            # Independent Head
            info_head = os.path.join(self.cat_dir, f"{sample}.head.info.txt")
            cmd_head = ["cutadapt", "-b", fwd_primer, "--revcomp", "-e", "0.25", "--overlap", "10", "--info-file", info_head, "-o", "/dev/null", f]
            subprocess.run(cmd_head, check=True, capture_output=True)

            # Independent Tail
            info_tail = os.path.join(self.cat_dir, f"{sample}.tail.info.txt")
            cmd_tail = ["cutadapt", "-b", rev_rc, "--revcomp", "-e", "0.25", "--overlap", "10", "--info-file", info_tail, "-o", "/dev/null", f]
            subprocess.run(cmd_tail, check=True, capture_output=True)

            # Parse IDs
            head_ids = self._parse_info_file(info_head)
            tail_ids = self._parse_info_file(info_tail)

            # Merge
            all_reads = []
            stats = {"Both": 0, "Head-only": 0, "Tail-only": 0, "Neither": 0, "Total": 0}
            final_mapping = []

            with open(info_head, "r") as info:
                for line in info:
                    rid = line.split("\t")[0]
                    has_h = rid in head_ids
                    has_t = rid in tail_ids
                    stats["Total"] += 1
                    if has_h and has_t: cat = "Both"
                    elif has_h: cat = "Head-only"
                    elif has_t: cat = "Tail-only"
                    else: cat = "Neither"

                    stats[cat] += 1
                    final_mapping.append({"ReadID": rid, "Category": cat})

            pd.DataFrame(final_mapping).to_csv(os.path.join(self.cat_dir, f"{sample}.categorization.csv"), index=False)
            stats["Sample"] = sample
            summary_stats.append(stats)

        df_summary = pd.DataFrame(summary_stats)
        df_summary.to_csv(os.path.join(self.cat_dir, "categorization_summary.csv"), index=False)
        return df_summary

    def _parse_info_file(self, path):
        ids = set()
        with open(path, "r") as f:
            for line in f:
                p = line.split("\t")
                if len(p) > 1 and int(p[1]) > 0: ids.add(p[0])
        return ids

    # --- Step 2: Global Alignment ---
    def align_global(self, threads=8):
        print("\n>>> Phase 2: Global Alignment (Minimap2)")
        fastq_files = glob.glob(os.path.join(self.data_dir, "*.fastq*"))
        for f in fastq_files:
            sample = os.path.basename(f).split(".")[0]
            out_bam = os.path.join(self.align_dir, f"{sample}.genomic.sorted.bam")
            if os.path.exists(out_bam): continue

            print(f" Aligning {sample}...")
            temp_sam = out_bam + ".sam"
            cmd = ["minimap2", "-ax", "map-ont", "-t", str(threads), self.genome_fasta, f]
            with open(temp_sam, "w") as sam:
                subprocess.run(cmd, stdout=sam, check=True)

            subprocess.run(["samtools", "sort", "-@", str(threads), "-o", out_bam, temp_sam], check=True)
            subprocess.run(["samtools", "index", out_bam], check=True)
            os.remove(temp_sam)

    # --- Step 3: Specificity & Gene Mapping ---
    def profile_specificity(self):
        print("\n>>> Phase 3: Specificity Profiling & Gene Mapping")
        gene_map = self._load_gtf_bins()

        bam_files = glob.glob(os.path.join(self.align_dir, "*.bam"))
        all_metrics = []
        global_stats = {cat: Counter() for cat in ["Both", "Head-only", "Tail-only", "Neither"]}

        for bam in bam_files:
            sample = os.path.basename(bam).replace(".genomic.sorted.bam", "")
            cat_csv = os.path.join(self.cat_dir, f"{sample}.categorization.csv")
            df_cat = pd.read_csv(cat_csv)
            read_to_cat = dict(zip(df_cat["ReadID"], df_cat["Category"]))
            cat_totals = Counter(df_cat["Category"])

            sample_stats = {cat: Counter() for cat in ["Both", "Head-only", "Tail-only", "Neither"]}
            on_target_counts = Counter()

            with pysam.AlignmentFile(bam, "rb") as sam:
                for read in sam.fetch(until_eof=True):
                    if read.is_unmapped or read.query_name not in read_to_cat: continue
                    if read.is_secondary or read.is_supplementary: continue

                    cat = read_to_cat[read.query_name]
                    if read.reference_name == XIST_CHROM and read.reference_start < XIST_END and read.reference_end > XIST_START:
                        label = "Xist"
                        on_target_counts[cat] += 1
                    else:
                        bin_idx = read.reference_start // BIN_SIZE
                        gene_hits = gene_map.get((read.reference_name, bin_idx), [])
                        label = gene_hits[0] if gene_hits else f"{read.reference_name}:{bin_idx*BIN_SIZE/1e6:.1f}Mb"

                    sample_stats[cat][label] += 1
                    global_stats[cat][label] += 1

            # Summarize sample
            for cat in ["Both", "Head-only", "Tail-only", "Neither"]:
                top_hits = sample_stats[cat].most_common(10)
                genes_to_report = [g for g, c in top_hits]
                if "Xist" not in genes_to_report: genes_to_report.append("Xist")

                for gene in genes_to_report:
                    count = sample_stats[cat][gene]
                    all_metrics.append({
                        "Sample": sample, "Category": cat, "Gene": gene, "Count": int(count),
                        "OnTarget": (gene == "Xist"), "CatTotal": int(cat_totals[cat])
                    })

        df_results = pd.DataFrame(all_metrics)
        df_results.to_csv(os.path.join(self.spec_dir, "gene_specificity_metrics.csv"), index=False)
        self._generate_plots(df_results)
        self._generate_step3_report(df_results, global_stats)

    def _load_gtf_bins(self):
        print(f" Loading GTF for gene mapping: {os.path.basename(self.gtf_path)}...")
        gene_map = defaultdict(list)
        with gzip.open(self.gtf_path, "rt") as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.split("\t")
                if parts[2] != "transcript": continue

                chrom, start, end = parts[0], int(parts[3]), int(parts[4])
                attrs = parts[8]
                gene_name = "Unknown"
                if 'gene_name "' in attrs:
                    gene_name = attrs.split('gene_name "')[1].split('"')[0]

                # Binning
                for b in range(start // BIN_SIZE, (end // BIN_SIZE) + 1):
                    if gene_name not in gene_map[(chrom, b)]:
                        gene_map[(chrom, b)].append(gene_name)
        return gene_map

    def _generate_plots(self, df):
        # 1. Histogram of top accumulation for a representative sample
        # We pick the sample with the most 'Tail-only' hits as representative
        tail_counts = df[df['Category'] == 'Tail-only'].groupby('Sample')['Count'].sum()
        if tail_counts.empty: return
        best_sample = tail_counts.idxmax()

        plt.figure(figsize=(12, 7))
        # Top 10 genes for "Neither" category (the majority of noise)
        plot_data = df[(df['Sample'] == best_sample) & (df['Category'] == 'Neither')].sort_values('Count', ascending=False).head(10)

        if not plot_data.empty:
            bars = plt.bar(plot_data['Gene'], plot_data['Count'], color='#e74c3c', edgecolor='black', alpha=0.8)
            plt.title(f"Off-Target Accumulation Histogram: {best_sample}\n(Top 10 Loci for 'Neither' Category)", fontsize=14, fontweight='bold')
            plt.ylabel("Read Count", fontsize=12)
            plt.xlabel("Target Gene / Genomic Locus", fontsize=12)
            plt.xticks(rotation=45, ha='right')
            plt.grid(axis='y', linestyle='--', alpha=0.7)

            # Add counts on top of bars
            for bar in bars:
                height = bar.get_height()
                plt.text(bar.get_x() + bar.get_width()/2., height + 5, f'{int(height)}', ha='center', va='bottom', fontsize=10)

            plt.tight_layout()
            plt.savefig(os.path.join(self.graph_dir, "top_offtarget_histogram.png"), dpi=300)

        # 2. Xist Recovery Plot (On-Target)
        plt.figure(figsize=(14, 7))
        # Aggregate counts per sample and category
        xist_df = df[df['OnTarget']].pivot_table(index='Sample', columns='Category', values='Count', aggfunc='sum', fill_value=0)

        # Ensure all categories are present for consistent plotting
        for cat in ["Both", "Head-only", "Tail-only", "Neither"]:
            if cat not in xist_df.columns:
                xist_df[cat] = 0

        xist_df = xist_df[["Both", "Head-only", "Tail-only", "Neither"]]
        xist_df.plot(kind='bar', stacked=True, figsize=(14, 7), colormap='viridis', edgecolor='black', alpha=0.9)

        plt.title("Xist (On-Target) Recovery Across 8 Conditions", fontsize=16, fontweight='bold')
        plt.ylabel("Read Count", fontsize=13)
        plt.xlabel("Experimental Condition", fontsize=13)
        plt.legend(title="Read Category", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', linestyle='--', alpha=0.5)
        plt.tight_layout()
        plt.savefig(os.path.join(self.graph_dir, "xist_recovery_comparison.png"), dpi=300)

    def _generate_step3_report(self, df, global_stats):
        # Determine aggregate top hits
        report = f"""# Step 3 Validation: Global Alignment & Specificity Report

This report evaluates the genomic destination of reads aggregated across all 8 experimental conditions.

## 1. Global Off-Target Accumulation Profile (Aggregated across 8 conditions)

The following histogram identifying the primary genomic hotspots where non-specific reads accumulate across the entire library.

![Top Off-Target Histogram](../graphs/top_offtarget_histogram.png)

### Top Loci per Category (Aggregated Metrics)

| Category | Gene / Genomic Locus | Total Read Count | **On-Target (Xist)?** |
| :--- | :--- | :---: | :---: |
"""
        for cat in ["Both", "Head-only", "Tail-only", "Neither"]:
            top_hits = global_stats[cat].most_common(5)
            # Ensure Xist is shown if it has ANY hits but isn't in top 5
            genes_to_show = [g for g, c in top_hits]
            if global_stats[cat]["Xist"] > 0 and "Xist" not in genes_to_show:
                genes_to_show.append("Xist")

            if not genes_to_show:
                report += f"| {cat} | - | 0 | - |\n"
            else:
                # Re-sort to show highest on top
                sorted_hits = sorted([(g, global_stats[cat][g]) for g in genes_to_show], key=lambda x: x[1], reverse=True)
                for gene, count in sorted_hits:
                    report += f"| {cat} | {gene} | {count} | {'**YES**' if gene == 'Xist' else 'No'} |\n"

        report += f"\n## 2. On-Target (*Xist*) Recovery Summary Across Conditions\n\n![Xist Recovery](../graphs/xist_recovery_comparison.png)\n\n"
        report += "### Per-Sample Xist Hits\n\n"
        report += "| Sample | Complete (Both) | Tail-only | Head-only | Neither | **Total Xist Hits** |\n"
        report += "| :--- | :---: | :---: | :---: | :---: | :---: |\n"

        xist_pivot = df[df['OnTarget']].pivot_table(index='Sample', columns='Category', values='Count', aggfunc='sum', fill_value=0)
        # Ensure all columns exist
        for cat in ["Both", "Tail-only", "Head-only", "Neither"]:
            if cat not in xist_pivot.columns: xist_pivot[cat] = 0

        for sample, row in xist_pivot.iterrows():
            total = row.sum()
            report += f"| {sample} | {int(row['Both'])} | {int(row['Tail-only'])} | {int(row['Head-only'])} | {int(row['Neither'])} | **{int(total)}** |\n"

        report += "\n> [!NOTE]\n> The 'Neither' category represents reads that map to the *Xist* locus but do not contain identifiable primer/barcode sequences within current error tolerances."

        with open(os.path.join(self.report_dir, "Step3_Alignment_Report.md"), "w") as f:
            f.write(report)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["all", "init", "categorize", "align", "profile"])
    parser.add_argument("--data_dir", required=True)
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--f_primer", required=True)
    parser.add_argument("--r_primer", required=True)
    args = parser.parse_args()

    pipeline = NanoAmpliconPipeline(args.data_dir, args.results_dir, args.f_primer, args.r_primer)

    if args.command in ["all", "categorize"]:
        pipeline.categorize_reads()
    if args.command in ["all", "align"]:
        pipeline.align_global()
    if args.command in ["all", "profile"]:
        pipeline.profile_specificity()

if __name__ == "__main__":
    main()
