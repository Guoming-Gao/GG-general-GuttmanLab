#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam


ROOT = Path("/Volumes/guttman/users/gmgao/Data_seq/20260629-AmpSeq-KOveri on polyclonal-cross-junction SNPs-aware on WT vs Dox")
ORIGIN = ROOT / "origin_analysis"
REFS = ORIGIN / "references"
LOGS = ORIGIN / "logs"
FIGS = ORIGIN / "figures"
NCRNA_ALIGN = ORIGIN / "ncRNA_alignment"
RNA_STAR = ORIGIN / "RNA_mm10_star"
RNA_STAR_ALIGNED = RNA_STAR / "aligned"
RNA_MM10 = ORIGIN / "RNA_mm10_bowtie2"
RNA_MM10_ALIGNED = RNA_MM10 / "aligned"
GDNA = ROOT / "gDNA_reanalysis"
MM10_GDNA = GDNA / "mm10_alignment" / "aligned"

MM10_NCRNA_FA = Path("/Volumes/guttman/genomes/mm10/fasta/mm10_ncRNA.fa")
MM10_NCRNA_INDEX = REFS / "mm10_ncRNA"
MM10_BOWTIE2_INDEX = "/Volumes/guttman/genomes/mm10/bowtie2_index/mm10"
STAR_INDEX = Path("/private/tmp/origin_analysis_STAR_mm10_2.7.11b")
MM10_FASTA = Path("/Volumes/guttman/genomes/mm10/fasta/mm10.fa")
GTF = Path("/Volumes/guttman/genomes/mm10/annotation/gencode.vM25.annotation.gtf")

SAMPLES = ["gDNA_WT", "gDNA_polyclonal_dRNF12", "gDNA_polyclonal_dTsix", "RNA_WT_diff", "RNA_Dox_24h"]
GDNA_SAMPLES = ["gDNA_WT", "gDNA_polyclonal_dRNF12", "gDNA_polyclonal_dTsix"]
RNA_SAMPLES = ["RNA_WT_diff", "RNA_Dox_24h"]
INTENDED_RNA_GENES = {"Xist", "Kdm5c", "Kdm6a", "Tsix", "Rlim", "Rbmx"}
BIN = 100_000

LOCI = {
    "RNF12_Koveri_P2": ("chrX", 103965939, 103966938),
    "Tsix_Koveri_region": ("chrX", 103447173, 103448802),
}


def ensure_dirs() -> None:
    for path in [ORIGIN, REFS, LOGS, FIGS, NCRNA_ALIGN, RNA_STAR, RNA_STAR_ALIGNED, RNA_MM10, RNA_MM10_ALIGNED]:
        path.mkdir(parents=True, exist_ok=True)


def pct(num: int | float, den: int | float) -> str:
    return f"{(100 * num / den):.2f}" if den else "0.00"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run_command(args: list[str], log_path: Path, stdout_path: Path | None = None, cwd: Path | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        if stdout_path:
            with stdout_path.open("wb") as out:
                proc = subprocess.run(args, stdout=out, stderr=log, cwd=cwd)
        else:
            proc = subprocess.run(args, stdout=log, stderr=log, cwd=cwd)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(args)}; see {log_path}")


def load_sample_manifest() -> dict[str, dict[str, str]]:
    return {row["sample"]: row for row in read_csv(ROOT / "metadata" / "sample_manifest.csv")}


def load_total_read_pairs() -> dict[str, int]:
    totals: dict[str, int] = {}
    for row in read_csv(ROOT / "qc" / "fastp" / "fastp_summary.csv"):
        totals[row["sample"]] = int(int(row["total_reads"]) / 2)
    return totals


def fastq_ids(r1_path: str) -> set[str]:
    ids: set[str] = set()
    opener = gzip.open if r1_path.endswith(".gz") else open
    with opener(r1_path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            handle.readline()
            handle.readline()
            handle.readline()
            name = header.strip().split()[0].lstrip("@")
            if name.endswith("/1") or name.endswith("/2"):
                name = name[:-2]
            ids.add(name)
    return ids


def copy_driver_script() -> None:
    try:
        source = Path(__file__).resolve()
        target = ORIGIN / "run_origin_analysis.py"
        if source != target:
            shutil.copy2(source, target)
    except Exception:
        pass


def ensure_ncRNA_index() -> None:
    if (Path(str(MM10_NCRNA_INDEX) + ".1.bt2")).exists():
        return
    run_command(
        ["bowtie2-build", str(MM10_NCRNA_FA), str(MM10_NCRNA_INDEX)],
        LOGS / "mm10_ncRNA.bowtie2_build.log",
    )


def ensure_star_index() -> None:
    if (STAR_INDEX / "Genome").exists() and (STAR_INDEX / "SA").exists():
        return
    if STAR_INDEX.exists():
        shutil.rmtree(STAR_INDEX)
    STAR_INDEX.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            "STAR",
            "--runMode",
            "genomeGenerate",
            "--runThreadN",
            "4",
            "--genomeDir",
            str(STAR_INDEX),
            "--genomeFastaFiles",
            str(MM10_FASTA),
            "--sjdbGTFfile",
            str(GTF),
            "--sjdbOverhang",
            "129",
        ],
        LOGS / "STAR_mm10_2.7.11b_genomeGenerate.log",
    )


def align_ncRNA(samples: dict[str, dict[str, str]]) -> dict[str, Path]:
    bam_paths: dict[str, Path] = {}
    for sample in SAMPLES:
        sorted_bam = NCRNA_ALIGN / f"{sample}.mm10_ncRNA.sorted.bam"
        if sorted_bam.exists() and Path(str(sorted_bam) + ".bai").exists():
            bam_paths[sample] = sorted_bam
            continue
        sam = NCRNA_ALIGN / f"{sample}.mm10_ncRNA.sam"
        bam = NCRNA_ALIGN / f"{sample}.mm10_ncRNA.bam"
        run_command(
            [
                "bowtie2",
                "--very-sensitive-local",
                "--no-unal",
                "-x",
                str(MM10_NCRNA_INDEX),
                "-1",
                samples[sample]["r1"],
                "-2",
                samples[sample]["r2"],
                "-S",
                str(sam),
            ],
            LOGS / f"{sample}.mm10_ncRNA.bowtie2.log",
        )
        run_command(["samtools", "view", "-bS", str(sam), "-o", str(bam)], LOGS / f"{sample}.mm10_ncRNA.samtools_view.log")
        run_command(["samtools", "sort", "-o", str(sorted_bam), str(bam)], LOGS / f"{sample}.mm10_ncRNA.samtools_sort.log")
        run_command(["samtools", "index", str(sorted_bam)], LOGS / f"{sample}.mm10_ncRNA.samtools_index.log")
        sam.unlink(missing_ok=True)
        bam.unlink(missing_ok=True)
        bam_paths[sample] = sorted_bam
    return bam_paths


def align_rna_mm10_bowtie2(samples: dict[str, dict[str, str]]) -> dict[str, Path]:
    bam_paths: dict[str, Path] = {}
    for sample in RNA_SAMPLES:
        sorted_bam = RNA_MM10_ALIGNED / f"{sample}.mm10.local.sorted.bam"
        if sorted_bam.exists() and Path(str(sorted_bam) + ".bai").exists():
            bam_paths[sample] = sorted_bam
            continue
        sam = RNA_MM10_ALIGNED / f"{sample}.mm10.local.sam"
        bam = RNA_MM10_ALIGNED / f"{sample}.mm10.local.bam"
        run_command(
            [
                "bowtie2",
                "--very-sensitive-local",
                "-x",
                MM10_BOWTIE2_INDEX,
                "-1",
                samples[sample]["r1"],
                "-2",
                samples[sample]["r2"],
                "-S",
                str(sam),
            ],
            LOGS / f"{sample}.mm10_local.bowtie2.log",
        )
        run_command(["samtools", "view", "-bS", str(sam), "-o", str(bam)], LOGS / f"{sample}.mm10_local.samtools_view.log")
        run_command(["samtools", "sort", "-o", str(sorted_bam), str(bam)], LOGS / f"{sample}.mm10_local.samtools_sort.log")
        run_command(["samtools", "index", str(sorted_bam)], LOGS / f"{sample}.mm10_local.samtools_index.log")
        sam.unlink(missing_ok=True)
        bam.unlink(missing_ok=True)
        bam_paths[sample] = sorted_bam
    return bam_paths


def group_primary_reads(bam_path: Path) -> dict[str, list[pysam.AlignedSegment]]:
    grouped: dict[str, list[pysam.AlignedSegment]] = defaultdict(list)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            grouped[read.query_name].append(read)
    return grouped


def overlaps(read: pysam.AlignedSegment, chrom: str, start1: int, end1: int) -> bool:
    if read.is_unmapped or read.reference_name != chrom or read.reference_end is None:
        return False
    return read.reference_start < end1 and read.reference_end > start1 - 1


def gdna_mm10_origin(reads: list[pysam.AlignedSegment]) -> str:
    mapped = [r for r in reads if not r.is_unmapped]
    if not mapped:
        return "unmapped_to_mm10"
    if max(r.mapping_quality for r in mapped) < 10:
        return "multi_mapping_or_low_MAPQ"
    if any(overlaps(r, *LOCI["RNF12_Koveri_P2"]) for r in mapped):
        return "intended_RNF12_locus"
    if any(overlaps(r, *LOCI["Tsix_Koveri_region"]) for r in mapped):
        return "intended_Tsix_locus"
    if any(r.reference_name == "chrX" for r in mapped):
        return "other_chrX"
    return "autosomal_or_other_contig"


def ncRNA_best_hit(reads: list[pysam.AlignedSegment]) -> tuple[str, int]:
    mapped = [r for r in reads if not r.is_unmapped]
    if not mapped:
        return "", -1
    best = max(mapped, key=lambda r: (r.mapping_quality, r.query_alignment_length or 0))
    return best.reference_name or "", int(best.mapping_quality)


def bam_mapped_query_names(bam_path: Path) -> set[str]:
    names: set[str] = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            names.add(read.query_name)
    return names


def ensure_star_rna_alignments(samples: dict[str, dict[str, str]]) -> dict[str, Path]:
    if not shutil.which("STAR"):
        raise RuntimeError("STAR is not available on PATH inside this environment.")
    bam_paths: dict[str, Path] = {}
    for sample in RNA_SAMPLES:
        final_bam = RNA_STAR_ALIGNED / f"{sample}.mm10.STAR.sorted.bam"
        if final_bam.exists() and Path(str(final_bam) + ".bai").exists():
            bam_paths[sample] = final_bam
            continue
        out_dir = RNA_STAR / sample
        out_dir.mkdir(parents=True, exist_ok=True)
        prefix = str(out_dir / f"{sample}.")
        out_tmp = Path("/private/tmp") / f"origin_analysis_STARtmp_{sample}"
        if out_tmp.exists():
            shutil.rmtree(out_tmp)
        tmp_r1 = Path("/private/tmp") / f"origin_analysis_{sample}_R1.fastq"
        tmp_r2 = Path("/private/tmp") / f"origin_analysis_{sample}_R2.fastq"
        for source, target in [(samples[sample]["r1"], tmp_r1), (samples[sample]["r2"], tmp_r2)]:
            with gzip.open(source, "rb") as src, target.open("wb") as dst:
                shutil.copyfileobj(src, dst)
        try:
            run_command(
                [
                    "STAR",
                    "--runThreadN",
                    "4",
                    "--genomeDir",
                    str(STAR_INDEX),
                    "--readFilesIn",
                    str(tmp_r1),
                    str(tmp_r2),
                    "--outFileNamePrefix",
                    prefix,
                    "--outTmpDir",
                    str(out_tmp),
                    "--outSAMtype",
                    "BAM",
                    "SortedByCoordinate",
                    "--outSAMunmapped",
                    "Within",
                    "--outSAMattributes",
                    "NH",
                    "HI",
                    "AS",
                    "nM",
                    "XS",
                ],
                LOGS / f"{sample}.STAR.log",
                cwd=out_dir,
            )
        finally:
            tmp_r1.unlink(missing_ok=True)
            tmp_r2.unlink(missing_ok=True)
        star_bam = out_dir / f"{sample}.Aligned.sortedByCoord.out.bam"
        if not star_bam.exists():
            raise RuntimeError(f"STAR did not create expected BAM: {star_bam}")
        star_bam.replace(final_bam)
        shutil.rmtree(out_tmp, ignore_errors=True)
        run_command(["samtools", "index", str(final_bam)], LOGS / f"{sample}.STAR.samtools_index.log")
        bam_paths[sample] = final_bam
    return bam_paths


ATTR_RE = re.compile(r'(\S+) "([^"]+)"')


def parse_gtf_gene_index() -> dict[str, dict[int, list[tuple[int, int, str]]]]:
    index: dict[str, dict[int, list[tuple[int, int, str]]]] = defaultdict(lambda: defaultdict(list))
    with GTF.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom, start_s, end_s, attrs_s = fields[0], fields[3], fields[4], fields[8]
            attrs = dict(ATTR_RE.findall(attrs_s))
            gene = attrs.get("gene_name") or attrs.get("gene_id")
            if not gene:
                continue
            start, end = int(start_s), int(end_s)
            for bin_id in range((start - 1) // BIN, (end - 1) // BIN + 1):
                index[chrom][bin_id].append((start, end, gene))
    return index


def gene_hits_for_reads(reads: list[pysam.AlignedSegment], gene_index: dict[str, dict[int, list[tuple[int, int, str]]]]) -> set[str]:
    genes: set[str] = set()
    for read in reads:
        if read.is_unmapped or read.mapping_quality < 10:
            continue
        chrom = read.reference_name
        for start0, end0 in read.get_blocks():
            start, end = start0 + 1, end0
            for bin_id in range((start - 1) // BIN, (end - 1) // BIN + 1):
                for gene_start, gene_end, gene in gene_index.get(chrom, {}).get(bin_id, []):
                    if gene_start <= end and gene_end >= start:
                        genes.add(gene)
    return genes


def rna_mm10_origin(reads: list[pysam.AlignedSegment], gene_hits: set[str]) -> str:
    mapped = [r for r in reads if not r.is_unmapped]
    if not mapped:
        return "unmapped_to_mm10"
    if max(r.mapping_quality for r in mapped) < 10:
        return "multi_mapping_or_low_MAPQ"
    if gene_hits & INTENDED_RNA_GENES:
        return "intended_RNA_target_gene_locus"
    if any(gene_hits) and any(r.reference_name == "chrX" for r in mapped if not r.is_unmapped):
        return "other_chrX_gene"
    if gene_hits:
        return "autosomal_or_other_gene"
    return "intergenic_no_gene_overlap"


def summarize_gdna(
    sample_ids: dict[str, set[str]],
    totals: dict[str, int],
    ncRNA_groups: dict[str, dict[str, list[pysam.AlignedSegment]]],
) -> tuple[list[dict], list[dict], Counter]:
    assigned: dict[str, set[str]] = defaultdict(set)
    for row in read_csv(GDNA / "classification" / "gDNA_fragment_KO_classification.csv"):
        assigned[row["sample"]].add(row["read_pair_id"])

    summary_rows: list[dict] = []
    assignment_rows: list[dict] = []
    top_ncRNA: Counter = Counter()

    for sample in GDNA_SAMPLES:
        mm10_grouped = group_primary_reads(MM10_GDNA / f"{sample}.mm10.sorted.bam")
        not_assigned = sample_ids[sample] - assigned[sample]
        counts: Counter = Counter()
        for qname in sorted(sample_ids[sample]):
            if qname in assigned[sample]:
                final_origin = "gDNA_expected_amplicon_assigned"
                mm10_origin = ""
                ncRNA_origin = ""
                contig = ""
                mapq = ""
            else:
                mm10_origin = gdna_mm10_origin(mm10_grouped.get(qname, []))
                contig, nc_mapq = ncRNA_best_hit(ncRNA_groups[sample].get(qname, []))
                if mm10_origin == "unmapped_to_mm10" and contig:
                    final_origin = "mapped_to_mm10_ncRNA_reference"
                    ncRNA_origin = "mapped_to_mm10_ncRNA_reference"
                    top_ncRNA[(sample, "gDNA", contig)] += 1
                    mapq = nc_mapq
                elif mm10_origin == "unmapped_to_mm10":
                    final_origin = "still_unmapped_after_mm10_and_ncRNA"
                    ncRNA_origin = "not_mapped_to_mm10_ncRNA_reference"
                    mapq = ""
                else:
                    final_origin = mm10_origin
                    ncRNA_origin = "not_tested_mm10_resolved"
                    mapq = ""
                counts[final_origin] += 1
            assignment_rows.append(
                {
                    "sample": sample,
                    "group": "gDNA",
                    "read_pair_id": qname,
                    "assignment_scope": "gDNA_synthetic_expected_amplicon",
                    "is_primary_assigned": qname in assigned[sample],
                    "mm10_origin": mm10_origin,
                    "ncRNA_origin": ncRNA_origin,
                    "final_origin": final_origin,
                    "best_ncRNA_contig": contig,
                    "best_ncRNA_mapq": mapq,
                }
            )
        for origin in [
            "intended_RNF12_locus",
            "intended_Tsix_locus",
            "other_chrX",
            "autosomal_or_other_contig",
            "multi_mapping_or_low_MAPQ",
            "mapped_to_mm10_ncRNA_reference",
            "still_unmapped_after_mm10_and_ncRNA",
        ]:
            count = counts[origin]
            summary_rows.append(
                {
                    "sample": sample,
                    "origin": origin,
                    "read_pairs": count,
                    "synthetic_not_assigned_read_pairs": len(not_assigned),
                    "total_sample_read_pairs": totals[sample],
                    "percent_of_synthetic_not_assigned": pct(count, len(not_assigned)),
                    "percent_of_all_read_pairs": pct(count, totals[sample]),
                }
            )
    return summary_rows, assignment_rows, top_ncRNA


def summarize_rna(
    sample_ids: dict[str, set[str]],
    totals: dict[str, int],
    ncRNA_groups: dict[str, dict[str, list[pysam.AlignedSegment]]],
    star_bams: dict[str, Path],
    gene_index: dict[str, dict[int, list[tuple[int, int, str]]]],
) -> tuple[list[dict], list[dict], list[dict], Counter]:
    summary_rows: list[dict] = []
    assignment_rows: list[dict] = []
    top_genes_by_sample: dict[str, Counter] = defaultdict(Counter)
    top_ncRNA: Counter = Counter()

    for sample in RNA_SAMPLES:
        amplicon_assigned = bam_mapped_query_names(ROOT / "aligned" / f"{sample}.bam")
        star_grouped = group_primary_reads(star_bams[sample])
        not_amplicon = sample_ids[sample] - amplicon_assigned
        counts: Counter = Counter()
        for qname in sorted(sample_ids[sample]):
            if qname in amplicon_assigned:
                final_origin = "RNA_amplicon_assigned"
                mm10_origin = ""
                ncRNA_origin = ""
                contig = ""
                mapq = ""
            else:
                reads = star_grouped.get(qname, [])
                genes = gene_hits_for_reads(reads, gene_index)
                for gene in genes:
                    top_genes_by_sample[sample][gene] += 1
                mm10_origin = rna_mm10_origin(reads, genes)
                contig, nc_mapq = ncRNA_best_hit(ncRNA_groups[sample].get(qname, []))
                if mm10_origin == "unmapped_to_mm10" and contig:
                    final_origin = "mapped_to_mm10_ncRNA_reference"
                    ncRNA_origin = "mapped_to_mm10_ncRNA_reference"
                    top_ncRNA[(sample, "RNA", contig)] += 1
                    mapq = nc_mapq
                elif mm10_origin == "unmapped_to_mm10":
                    final_origin = "still_unmapped_after_mm10_and_ncRNA"
                    ncRNA_origin = "not_mapped_to_mm10_ncRNA_reference"
                    mapq = ""
                else:
                    final_origin = mm10_origin
                    ncRNA_origin = "not_tested_mm10_resolved"
                    mapq = ""
                counts[final_origin] += 1
            assignment_rows.append(
                {
                    "sample": sample,
                    "group": "RNA",
                    "read_pair_id": qname,
                    "assignment_scope": "RNA_expected_amplicon",
                    "is_primary_assigned": qname in amplicon_assigned,
                    "mm10_origin": mm10_origin,
                    "ncRNA_origin": ncRNA_origin,
                    "final_origin": final_origin,
                    "best_ncRNA_contig": contig,
                    "best_ncRNA_mapq": mapq,
                }
            )
        for origin in [
            "intended_RNA_target_gene_locus",
            "other_chrX_gene",
            "autosomal_or_other_gene",
            "intergenic_no_gene_overlap",
            "multi_mapping_or_low_MAPQ",
            "mapped_to_mm10_ncRNA_reference",
            "still_unmapped_after_mm10_and_ncRNA",
        ]:
            count = counts[origin]
            summary_rows.append(
                {
                    "sample": sample,
                    "origin": origin,
                    "read_pairs": count,
                    "not_RNA_amplicon_assigned_read_pairs": len(not_amplicon),
                    "total_sample_read_pairs": totals[sample],
                    "percent_of_not_RNA_amplicon_assigned": pct(count, len(not_amplicon)),
                    "percent_of_all_read_pairs": pct(count, totals[sample]),
                }
            )
    gene_rows: list[dict] = []
    for sample, counts in top_genes_by_sample.items():
        not_amp_total = next(int(r["not_RNA_amplicon_assigned_read_pairs"]) for r in summary_rows if r["sample"] == sample)
        for gene, count in counts.most_common(30):
            gene_rows.append(
                {
                    "sample": sample,
                    "gene_name": gene,
                    "read_pairs": count,
                    "percent_of_not_RNA_amplicon_assigned": pct(count, not_amp_total),
                }
            )
    return summary_rows, assignment_rows, gene_rows, top_ncRNA


def write_top_ncRNA(top_ncRNA: Counter, summary_lookup: dict[tuple[str, str], int]) -> None:
    rows: list[dict] = []
    for (sample, group, contig), count in top_ncRNA.most_common():
        total = summary_lookup.get((sample, group), 0)
        rows.append(
            {
                "sample": sample,
                "group": group,
                "ncRNA_contig": contig,
                "read_pairs": count,
                "percent_of_unresolved_input": pct(count, total),
            }
        )
    write_csv(
        ORIGIN / "top_ncRNA_origin_hits.csv",
        rows,
        ["sample", "group", "ncRNA_contig", "read_pairs", "percent_of_unresolved_input"],
    )


def plot_stacked(path: Path, rows: list[dict], samples: list[str], denominator_col: str, title: str) -> None:
    origins = list(dict.fromkeys(row["origin"] for row in rows))
    colors = plt.cm.tab20.colors
    fig, ax = plt.subplots(figsize=(10, 3.2))
    left = [0.0] * len(samples)
    for idx, origin in enumerate(origins):
        vals: list[float] = []
        labels: list[str] = []
        for sample in samples:
            row = next((r for r in rows if r["sample"] == sample and r["origin"] == origin), None)
            vals.append(float(row[f"percent_of_{denominator_col}"]) if row else 0.0)
            labels.append(row["read_pairs"] if row else "0")
        ax.barh(samples, vals, left=left, label=origin.replace("_", " "), color=colors[idx % len(colors)])
        for i, (val, label) in enumerate(zip(vals, labels)):
            if val >= 7:
                ax.text(left[i] + val / 2, i, label, ha="center", va="center", fontsize=6, color="white")
        left = [l + v for l, v in zip(left, vals)]
    ax.set_xlim(0, 100)
    ax.set_xlabel(f"% of {denominator_col.replace('_', ' ')}", fontsize=8)
    ax.set_title(title, fontsize=9, loc="left")
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_top_ncRNA() -> None:
    rows = read_csv(ORIGIN / "top_ncRNA_origin_hits.csv")
    top = rows[:20]
    if not top:
        return
    labels = [f'{r["sample"]}: {r["ncRNA_contig"]}' for r in top][::-1]
    vals = [int(r["read_pairs"]) for r in top][::-1]
    colors = ["#0072B2" if r["group"] == "gDNA" else "#D55E00" for r in top][::-1]
    fig, ax = plt.subplots(figsize=(9, max(3, len(top) * 0.24)))
    ax.barh(labels, vals, color=colors)
    ax.set_xlabel("read pairs", fontsize=8)
    ax.set_title("Top ncRNA-reference hits among mm10-unresolved reads", fontsize=9, loc="left")
    ax.tick_params(labelsize=6)
    fig.tight_layout()
    fig.savefig(FIGS / "top_ncRNA_origin_hits.png", dpi=220)
    plt.close(fig)


def append_report_section(path: Path, title: str, body: str) -> None:
    start = "<!-- origin-analysis-start -->"
    end = "<!-- origin-analysis-end -->"
    section = f"\n{start}\n\n## {title}\n\n{body.strip()}\n\n{end}\n"
    text = path.read_text()
    if start in text and end in text:
        pre = text.split(start)[0].rstrip()
        post = text.split(end, 1)[1].lstrip()
        path.write_text(pre + section + ("\n" + post if post else ""))
    else:
        path.write_text(text.rstrip() + "\n" + section)


def update_reports() -> None:
    gdna = read_csv(ORIGIN / "gDNA_origin_summary_with_ncRNA.csv")
    rna = read_csv(ORIGIN / "RNA_origin_summary_with_ncRNA.csv")

    def mini_table(rows: list[dict], sample_col: str) -> str:
        lines = ["|sample|mapped to ncRNA ref|still unresolved|", "|---|---:|---:|"]
        for sample in sorted(set(r["sample"] for r in rows)):
            nc = next((r for r in rows if r["sample"] == sample and r["origin"] == "mapped_to_mm10_ncRNA_reference"), None)
            st = next((r for r in rows if r["sample"] == sample and r["origin"] == "still_unmapped_after_mm10_and_ncRNA"), None)
            nc_text = f'{nc["read_pairs"]} ({nc[f"percent_of_{sample_col}"]}%)' if nc else "0 (0.00%)"
            st_text = f'{st["read_pairs"]} ({st[f"percent_of_{sample_col}"]}%)' if st else "0 (0.00%)"
            lines.append(f"|{sample}|{nc_text}|{st_text}|")
        return "\n".join(lines)

    caveat = (
        "`mm10_ncRNA.fa` is used here as a local ncRNA/repeat-reference test. "
        "Reads in the ncRNA category match that reference; reads that remain unresolved are not automatically proven to be low-complexity."
    )
    gdna_body = (
        "The gDNA origin audit now splits the previous `unmapped_to_mm10` category by remapping those unresolved reads to the local "
        "`mm10_ncRNA.fa` reference.\n\n"
        + mini_table(gdna, "synthetic_not_assigned")
        + "\n\n"
        + caveat
        + "\n\n![gDNA origin with ncRNA split](origin_analysis/figures/gDNA_origin_with_ncRNA_split.png)\n"
        + "\n![Top ncRNA hits](origin_analysis/figures/top_ncRNA_origin_hits.png)"
    )
    rna_body = (
        "RNA reads not assigned to the small RNA amplicon reference were mapped to mm10 with `bowtie2 --very-sensitive-local`, annotated by overlap with GENCODE vM25 genes, "
        "and then mm10-unresolved reads were tested against `mm10_ncRNA.fa`. This is an origin audit, not a replacement for transcript-aware allele quantification.\n\n"
        + mini_table(rna, "not_RNA_amplicon_assigned")
        + "\n\n"
        + caveat
        + "\n\n![RNA origin with ncRNA split](origin_analysis/figures/RNA_origin_with_ncRNA_split.png)\n"
        + "\n![Top ncRNA hits](origin_analysis/figures/top_ncRNA_origin_hits.png)"
    )
    append_report_section(ROOT / "gDNA_ampseq_KO_and_allele_report.md", "Origin Audit With ncRNA Reference", gdna_body)
    append_report_section(ROOT / "RNA_ampseq_processing_and_results_report.md", "Origin Audit With ncRNA Reference", rna_body)


def main() -> None:
    ensure_dirs()
    copy_driver_script()
    samples = load_sample_manifest()
    totals = load_total_read_pairs()
    sample_ids = {sample: fastq_ids(samples[sample]["r1"]) for sample in SAMPLES}
    ensure_ncRNA_index()
    ncRNA_bams = align_ncRNA(samples)
    ncRNA_groups = {sample: group_primary_reads(path) for sample, path in ncRNA_bams.items()}

    gdna_summary, gdna_assignments, gdna_top_ncRNA = summarize_gdna(sample_ids, totals, ncRNA_groups)
    write_csv(
        ORIGIN / "gDNA_origin_summary_with_ncRNA.csv",
        gdna_summary,
        [
            "sample",
            "origin",
            "read_pairs",
            "synthetic_not_assigned_read_pairs",
            "total_sample_read_pairs",
            "percent_of_synthetic_not_assigned",
            "percent_of_all_read_pairs",
        ],
    )

    rna_mm10_bams = align_rna_mm10_bowtie2(samples)
    gene_index = parse_gtf_gene_index()
    rna_summary, rna_assignments, gene_rows, rna_top_ncRNA = summarize_rna(sample_ids, totals, ncRNA_groups, rna_mm10_bams, gene_index)
    write_csv(
        ORIGIN / "RNA_origin_summary_with_ncRNA.csv",
        rna_summary,
        [
            "sample",
            "origin",
            "read_pairs",
            "not_RNA_amplicon_assigned_read_pairs",
            "total_sample_read_pairs",
            "percent_of_not_RNA_amplicon_assigned",
            "percent_of_all_read_pairs",
        ],
    )
    write_csv(
        ORIGIN / "RNA_top_gene_origin_hits.csv",
        gene_rows,
        ["sample", "gene_name", "read_pairs", "percent_of_not_RNA_amplicon_assigned"],
    )
    with (ORIGIN / "read_pair_origin_assignments.csv").open("w", newline="") as handle:
        fields = [
            "sample",
            "group",
            "read_pair_id",
            "assignment_scope",
            "is_primary_assigned",
            "mm10_origin",
            "ncRNA_origin",
            "final_origin",
            "best_ncRNA_contig",
            "best_ncRNA_mapq",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(gdna_assignments)
        writer.writerows(rna_assignments)

    unresolved_totals: dict[tuple[str, str], int] = {}
    for row in gdna_summary:
        if row["origin"] in {"mapped_to_mm10_ncRNA_reference", "still_unmapped_after_mm10_and_ncRNA"}:
            unresolved_totals[(row["sample"], "gDNA")] = unresolved_totals.get((row["sample"], "gDNA"), 0) + int(row["read_pairs"])
    for row in rna_summary:
        if row["origin"] in {"mapped_to_mm10_ncRNA_reference", "still_unmapped_after_mm10_and_ncRNA"}:
            unresolved_totals[(row["sample"], "RNA")] = unresolved_totals.get((row["sample"], "RNA"), 0) + int(row["read_pairs"])
    write_top_ncRNA(gdna_top_ncRNA + rna_top_ncRNA, unresolved_totals)

    plot_stacked(
        FIGS / "gDNA_origin_with_ncRNA_split.png",
        gdna_summary,
        GDNA_SAMPLES,
        "synthetic_not_assigned",
        "gDNA synthetic-not-assigned read origins with ncRNA split",
    )
    plot_stacked(
        FIGS / "RNA_origin_with_ncRNA_split.png",
        rna_summary,
        RNA_SAMPLES,
        "not_RNA_amplicon_assigned",
        "RNA not-amplicon read origins with ncRNA split",
    )
    plot_top_ncRNA()
    update_reports()


if __name__ == "__main__":
    main()
