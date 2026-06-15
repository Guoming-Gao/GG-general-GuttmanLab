"""High-level KOveri primer design orchestration."""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

import pandas as pd

from .guide_locator import locate_guides_exact, select_coherent_target
from .io_utils import read_guide_fasta, write_dataframe, write_fasta
from .models import SnpRecord, TargetRegion
from .post_filters import _annotate_row, annotate_post_filters
from .primer3_utils import run_primer3_pair_design
from .report import write_html_report
from .sequence_utils import extract_reference_sequence
from .snp_utils import SnpDatabase
from .specificity import check_primer_specificity
from .tools import require_file, resolve_tool


def design_koveri_primers(input_fasta: str | Path, output_dir: str | Path, config: dict) -> dict:
    """Run the full KOveri primer design pipeline."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    genome_fasta = require_file(config["genome_fasta"], "Genome FASTA")
    blast_db = require_file(config["blast_database"] + ".nsq", "BLAST database").removesuffix(".nsq")
    snp_vcf = require_file(config["snp_vcf"], "SNP VCF")
    samtools = resolve_tool(config, "samtools", "samtools")
    blastn = resolve_tool(config, "blastn", "blastn")
    primer3_core = resolve_tool(config, "primer3_core", "primer3_core")
    ntthal = resolve_tool(config, "ntthal", "ntthal")

    guides = read_guide_fasta(input_fasta)
    guide_hits, partial_counts = locate_guides_exact(guides, blast_db, blastn)
    target = select_coherent_target(guides, guide_hits, config["max_cluster_span"])

    snp_db = SnpDatabase(snp_vcf, config["vcf_b6_sample"], config["vcf_cast_sample"])
    try:
        all_rows = _design_candidates(
            target=target,
            snp_db=snp_db,
            genome_fasta=genome_fasta,
            samtools=samtools,
            blast_db=blast_db,
            blastn=blastn,
            primer3_core=primer3_core,
            config=config,
            ntthal=ntthal,
        )
    finally:
        snp_db.close()

    candidate_rows = _dedupe_candidates(all_rows)
    candidate_rows.sort(key=_rank_key)
    if config.get("blast_specificity", True) and candidate_rows:
        _annotate_specificity(candidate_rows, blast_db, blastn, config)
    filter_info = annotate_post_filters(candidate_rows, config, ntthal)
    candidate_rows.sort(key=_rank_key)
    top_rows = _select_top_rows(candidate_rows, filter_info["selected_tier"], int(config["top_n_output"]))
    if top_rows:
        candidate_rows.sort(key=_rank_key)

    guide_rows = [
        {
            "Guide_ID": hit.guide_id,
            "Guide_Seq": hit.sequence,
            "Chromosome": hit.chrom,
            "Guide_Start": hit.start,
            "Guide_End": hit.end,
            "Guide_Strand": hit.strand,
            "Exact_Hit_Count_For_Guide": sum(h.guide_id == hit.guide_id for h in guide_hits),
            "Partial_Hit_Count_For_Guide": partial_counts.get(hit.guide_id, 0),
            "Selected_For_Target": hit in target.guide_hits,
        }
        for hit in guide_hits
    ]
    write_dataframe(guide_rows, output_dir / "KOveri_guide_hits.csv")
    write_dataframe(candidate_rows, output_dir / "KOveri_primer_candidates.csv")
    _write_top_candidates(top_rows, candidate_rows, output_dir / "KOveri_primer_top_candidates.csv")
    _write_final_fasta(top_rows, output_dir / "KOveri_primers.fasta")
    _write_bed(top_rows, target, output_dir / "KOveri_amplicons.bed")
    write_html_report(
        output_dir / "KOveri_report.html",
        target=target,
        guide_hits=list(target.guide_hits),
        top_rows=top_rows,
        candidate_count=len(candidate_rows),
        config=config,
        filter_info=filter_info,
    )
    _write_summary(
        output_dir / "KOveri_summary.txt",
        input_fasta=input_fasta,
        output_dir=output_dir,
        target=target,
        guides=guides,
        guide_hits=guide_hits,
        rows=candidate_rows,
        top_rows=top_rows,
        config=config,
        filter_info=filter_info,
    )

    return {
        "target": target,
        "guide_hits": guide_hits,
        "candidate_rows": candidate_rows,
        "top_rows": top_rows,
        "filter_info": filter_info,
        "output_dir": output_dir,
    }


def _design_candidates(
    target: TargetRegion,
    snp_db: SnpDatabase,
    genome_fasta: str,
    samtools: str,
    blast_db: str,
    blastn: str,
    primer3_core: str,
    config: dict,
    ntthal: str,
) -> list[dict]:
    rows: list[dict] = []
    min_snps = int(config["min_informative_snps"])
    goals = list(range(min_snps, 0, -1))
    stop_after_goal_success = False
    thermo_cache = {}

    for snp_goal in goals:
        for flank in config["flank_steps"]:
            region_start = max(1, target.start - int(flank))
            region_end = target.end + int(flank)
            snps = snp_db.query_informative(target.chrom, region_start, region_end)
            spans = _target_spans(target, snps, snp_goal, config["max_target_spans_per_flank"])
            if not spans:
                continue

            template = extract_reference_sequence(target.chrom, region_start, region_end, genome_fasta, samtools)
            for span_index, (span_start, span_end, span_snps) in enumerate(spans):
                if span_end - span_start + 1 > int(config["max_amplicon_size"]):
                    continue
                sequence_id = f"{target.chrom}_{region_start}_{region_end}_flank{flank}_span{span_index}"
                try:
                    pairs = run_primer3_pair_design(
                        sequence_id=sequence_id,
                        template=template,
                        target_start_abs=span_start,
                        target_end_abs=span_end,
                        template_start_abs=region_start,
                        excluded_snps=snps,
                        config=config,
                        primer3_core=primer3_core,
                    )
                except Exception as exc:
                    rows.append(_failure_row(target, flank, snp_goal, span_start, span_end, str(exc)))
                    continue

                for pair in pairs:
                    row = _score_pair(
                        pair=pair,
                        target=target,
                        region_start=region_start,
                        region_end=region_end,
                        flank=flank,
                        snp_goal=snp_goal,
                        all_snps=snps,
                        target_span_snps=span_snps,
                        ideal_amplicon_max=int(config["ideal_amplicon_max"]),
                    )
                    if row["Primer_Overlaps_SNP"]:
                        row["Failure_Reason"] = "primer_overlaps_snp"
                    else:
                        row["Failure_Reason"] = ""

                    # Early hard-filter annotation lets the search keep expanding
                    # until it finds orderable GC-clamped candidates.
                    _annotate_row(row, config, ntthal, thermo_cache)
                    rows.append(row)

            if any(r.get("Hard_Filter_Pass") is True and r.get("SNP_Goal_For_Attempt") == snp_goal for r in rows):
                stop_after_goal_success = True
                break
        if stop_after_goal_success:
            break

    return rows


def _target_spans(
    target: TargetRegion,
    snps: list[SnpRecord],
    snp_goal: int,
    max_spans: int,
) -> list[tuple[int, int, list[SnpRecord]]]:
    if snp_goal <= 0:
        return [(target.start, target.end, [])]
    sorted_snps = sorted(snps, key=lambda s: s.pos)
    spans: list[tuple[int, int, list[SnpRecord]]] = []
    for i in range(0, len(sorted_snps) - snp_goal + 1):
        subset = sorted_snps[i : i + snp_goal]
        span_start = min(target.start, subset[0].pos)
        span_end = max(target.end, subset[-1].pos)
        spans.append((span_start, span_end, subset))

    # Also try denser windows with extra SNPs if they do not greatly expand the span.
    for i in range(len(sorted_snps)):
        for j in range(i + snp_goal, min(len(sorted_snps), i + snp_goal + 4) + 1):
            subset = sorted_snps[i:j]
            if len(subset) < snp_goal:
                continue
            span_start = min(target.start, subset[0].pos)
            span_end = max(target.end, subset[-1].pos)
            spans.append((span_start, span_end, subset))

    unique = OrderedDict()
    for span_start, span_end, subset in spans:
        unique[(span_start, span_end)] = (span_start, span_end, subset)
    ranked = sorted(
        unique.values(),
        key=lambda x: (x[1] - x[0] + 1, -len(x[2]), abs((x[0] + x[1]) / 2 - (target.start + target.end) / 2)),
    )
    return ranked[:max_spans]


def _score_pair(
    pair: dict,
    target: TargetRegion,
    region_start: int,
    region_end: int,
    flank: int,
    snp_goal: int,
    all_snps: list[SnpRecord],
    target_span_snps: list[SnpRecord],
    ideal_amplicon_max: int,
) -> dict:
    amp_start = pair["amplicon_start"]
    amp_end = pair["amplicon_end"]
    snps_in_amplicon = [s for s in all_snps if amp_start <= s.pos <= amp_end]
    left_overlaps = [
        s for s in all_snps if pair["left_genomic_start"] <= s.pos <= pair["left_genomic_end"]
    ]
    right_overlaps = [
        s for s in all_snps if pair["right_genomic_start"] <= s.pos <= pair["right_genomic_end"]
    ]
    covers_guides = amp_start <= target.start and amp_end >= target.end
    tm_delta = abs(pair["left_tm"] - pair["right_tm"])
    ideal_penalty = max(0, pair["amplicon_size"] - ideal_amplicon_max)
    return {
        "Chromosome": target.chrom,
        "Protected_Edit_Start": target.start,
        "Protected_Edit_End": target.end,
        "Protected_Edit_Length": target.length,
        "Design_Window_Start": region_start,
        "Design_Window_End": region_end,
        "Flank_Used": flank,
        "SNP_Goal_For_Attempt": snp_goal,
        "Target_SNPs_For_Attempt": _format_snps(target_span_snps),
        "Amplicon_Start": amp_start,
        "Amplicon_End": amp_end,
        "Amplicon_Size": pair["amplicon_size"],
        "Covers_All_Guides": covers_guides,
        "SNP_Count_In_Amplicon": len(snps_in_amplicon),
        "SNP_Positions_In_Amplicon": _format_snps(snps_in_amplicon),
        "Left_Primer_Seq": pair["left_sequence"],
        "Left_Primer_Start": pair["left_genomic_start"],
        "Left_Primer_End": pair["left_genomic_end"],
        "Left_Primer_Tm": pair["left_tm"],
        "Left_Primer_GC": pair["left_gc"],
        "Right_Primer_Seq": pair["right_sequence"],
        "Right_Primer_Start": pair["right_genomic_start"],
        "Right_Primer_End": pair["right_genomic_end"],
        "Right_Primer_Tm": pair["right_tm"],
        "Right_Primer_GC": pair["right_gc"],
        "Primer_Tm_Delta": tm_delta,
        "Primer3_Pair_Penalty": pair["pair_penalty"],
        "Primer_Overlaps_SNP": bool(left_overlaps or right_overlaps),
        "Left_Primer_SNP_Overlaps": _format_snps(left_overlaps),
        "Right_Primer_SNP_Overlaps": _format_snps(right_overlaps),
        "Ideal_Size_Overage": ideal_penalty,
        "Left_Primer_PerfectOrNear_Hits": "",
        "Left_Primer_Noise_Hits": "",
        "Left_Primer_Specificity_Pass": "",
        "Right_Primer_PerfectOrNear_Hits": "",
        "Right_Primer_Noise_Hits": "",
        "Right_Primer_Specificity_Pass": "",
        "Primer_Pair_Specificity_Pass": "",
    }


def _failure_row(
    target: TargetRegion,
    flank: int,
    snp_goal: int,
    span_start: int,
    span_end: int,
    reason: str,
) -> dict:
    return {
        "Chromosome": target.chrom,
        "Protected_Edit_Start": target.start,
        "Protected_Edit_End": target.end,
        "Flank_Used": flank,
        "SNP_Goal_For_Attempt": snp_goal,
        "Amplicon_Start": "",
        "Amplicon_End": "",
        "Amplicon_Size": "",
        "Covers_All_Guides": False,
        "SNP_Count_In_Amplicon": 0,
        "Left_Primer_Seq": "",
        "Right_Primer_Seq": "",
        "Primer_Overlaps_SNP": "",
        "Failure_Reason": f"primer3_failed_for_target_span_{span_start}_{span_end}: {reason}",
    }


def _format_snps(snps: list[SnpRecord]) -> str:
    return ";".join(f"{s.chrom}:{s.pos}:{s.ref}>{s.alt}:{s.genotype}" for s in snps)


def _dedupe_candidates(rows: list[dict]) -> list[dict]:
    seen = set()
    deduped = []
    for row in rows:
        key = (
            row.get("Left_Primer_Seq"),
            row.get("Right_Primer_Seq"),
            row.get("Amplicon_Start"),
            row.get("Amplicon_End"),
        )
        if not key[0] or key in seen:
            continue
        seen.add(key)
        deduped.append(row)
    return deduped


def _rank_key(row: dict) -> tuple:
    specificity = row.get("Primer_Pair_Specificity_Pass")
    specificity_pass = specificity is True or specificity == "True"
    tier_rank = {
        "strict_all": 0,
        "relaxed_amplicon": 1,
        "relaxed_specificity": 2,
        "relaxed_gc_tm": 3,
        "failed": 4,
    }.get(row.get("Filter_Tier"), 5)
    return (
        tier_rank,
        not row.get("Covers_All_Guides", False),
        row.get("Primer_Overlaps_SNP", True),
        -int(row.get("SNP_Count_In_Amplicon") or 0),
        not specificity_pass,
        int(row.get("Ideal_Size_Overage") or 0),
        int(row.get("Amplicon_Size") or 999999),
        float(row.get("Primer_Tm_Delta") or 999),
        float(row.get("Primer3_Pair_Penalty") or 999),
    )


def _select_top_rows(rows: list[dict], selected_tier: str, top_n: int) -> list[dict]:
    if selected_tier == "failed":
        return []
    return [row for row in rows if row.get("Filter_Tier") == selected_tier][:top_n]


def _write_top_candidates(top_rows: list[dict], candidate_rows: list[dict], path: Path) -> None:
    if top_rows:
        pd.DataFrame(top_rows).to_csv(path, index=False)
        return
    columns = list(candidate_rows[0].keys()) if candidate_rows else []
    pd.DataFrame(columns=columns).to_csv(path, index=False)


def _annotate_specificity(rows: list[dict], blast_db: str, blastn: str, config: dict) -> None:
    cache = {}
    for row in rows:
        left_seq = row.get("Left_Primer_Seq")
        right_seq = row.get("Right_Primer_Seq")
        if not left_seq or not right_seq:
            continue
        if left_seq not in cache:
            cache[left_seq] = check_primer_specificity(left_seq, blast_db, blastn, config)
        if right_seq not in cache:
            cache[right_seq] = check_primer_specificity(right_seq, blast_db, blastn, config)
        left_spec = cache[left_seq]
        right_spec = cache[right_seq]
        row.update(
            {
                "Left_Primer_PerfectOrNear_Hits": left_spec.perfect_or_near_hits,
                "Left_Primer_Noise_Hits": left_spec.noise_hits,
                "Left_Primer_Specificity_Pass": left_spec.pass_specificity,
                "Right_Primer_PerfectOrNear_Hits": right_spec.perfect_or_near_hits,
                "Right_Primer_Noise_Hits": right_spec.noise_hits,
                "Right_Primer_Specificity_Pass": right_spec.pass_specificity,
                "Primer_Pair_Specificity_Pass": (
                    left_spec.pass_specificity and right_spec.pass_specificity
                ),
            }
        )


def _write_final_fasta(rows: list[dict], path: Path) -> None:
    records: list[tuple[str, str]] = []
    for idx, row in enumerate(rows, start=1):
        records.append((f"KOveri_pair{idx}_LEFT_{row['Chromosome']}:{row['Left_Primer_Start']}-{row['Left_Primer_End']}", row["Left_Primer_Seq"]))
        records.append((f"KOveri_pair{idx}_RIGHT_{row['Chromosome']}:{row['Right_Primer_Start']}-{row['Right_Primer_End']}", row["Right_Primer_Seq"]))
    write_fasta(records, path)


def _write_bed(rows: list[dict], target: TargetRegion, path: Path) -> None:
    with path.open("w") as handle:
        handle.write("#chrom\tstart0\tend1\tname\tscore\tstrand\tfeature\n")
        handle.write(
            f"{target.chrom}\t{target.start - 1}\t{target.end}\tprotected_edit_interval\t0\t.\tguide_interval\n"
        )
        for idx, row in enumerate(rows, start=1):
            handle.write(
                f"{row['Chromosome']}\t{int(row['Amplicon_Start']) - 1}\t{row['Amplicon_End']}\t"
                f"KOveri_pair{idx}_amplicon\t{row['SNP_Count_In_Amplicon']}\t.\tamplicon\n"
            )
            handle.write(
                f"{row['Chromosome']}\t{int(row['Left_Primer_Start']) - 1}\t{row['Left_Primer_End']}\t"
                f"KOveri_pair{idx}_left_primer\t0\t+\tprimer\n"
            )
            handle.write(
                f"{row['Chromosome']}\t{int(row['Right_Primer_Start']) - 1}\t{row['Right_Primer_End']}\t"
                f"KOveri_pair{idx}_right_primer\t0\t-\tprimer\n"
            )


def _write_summary(
    path: Path,
    input_fasta: str | Path,
    output_dir: Path,
    target: TargetRegion,
    guides: list,
    guide_hits: list,
    rows: list[dict],
    top_rows: list[dict],
    config: dict,
    filter_info: dict,
) -> None:
    best = top_rows[0] if top_rows else None
    with path.open("w") as handle:
        handle.write("KOveri primer design summary\n")
        handle.write("============================\n\n")
        handle.write(f"Input FASTA: {input_fasta}\n")
        handle.write(f"Output directory: {output_dir}\n")
        handle.write(f"Guides: {len(guides)}\n")
        handle.write(f"Exact guide hits: {len(guide_hits)}\n")
        handle.write(f"Selected protected edit interval: {target.chrom}:{target.start}-{target.end}\n")
        handle.write(f"Minimum informative SNP goal: {config['min_informative_snps']}\n")
        handle.write(f"Candidate primer pairs: {len(rows)}\n")
        handle.write(f"Top primer pairs reported: {len(top_rows)}\n\n")
        handle.write("Post-Primer3 filters\n")
        handle.write("--------------------\n")
        handle.write(f"Selected filter tier: {filter_info['selected_tier']}\n")
        handle.write(f"Tier counts: {filter_info['tier_counts']}\n")
        handle.write(f"Hard-filter passing candidates: {filter_info['hard_pass_count']}\n")
        handle.write(f"Hard-filter failing candidates: {filter_info['hard_fail_count']}\n")
        handle.write("Tier explanations:\n")
        for tier, explanation in filter_info["tier_explanations"].items():
            handle.write(f"  - {tier}: {explanation}\n")
        handle.write("\n")
        if best:
            handle.write("Best ranked primer pair\n")
            handle.write("-----------------------\n")
            handle.write(f"Left primer:  {best['Left_Primer_Seq']} ({best['Chromosome']}:{best['Left_Primer_Start']}-{best['Left_Primer_End']})\n")
            handle.write(f"Right primer: {best['Right_Primer_Seq']} ({best['Chromosome']}:{best['Right_Primer_Start']}-{best['Right_Primer_End']})\n")
            handle.write(f"Amplicon: {best['Chromosome']}:{best['Amplicon_Start']}-{best['Amplicon_End']} ({best['Amplicon_Size']} bp)\n")
            handle.write(f"Informative SNPs in amplicon: {best['SNP_Count_In_Amplicon']}\n")
            handle.write(f"Primer overlaps informative SNP: {best['Primer_Overlaps_SNP']}\n")
            handle.write(f"Primer pair specificity pass: {best.get('Primer_Pair_Specificity_Pass')}\n")
            handle.write(f"Filter tier: {best.get('Filter_Tier')}\n")
            handle.write(f"Relaxed criteria used: {best.get('Relaxed_Criteria_Used', '')}\n")
            handle.write(f"Worst thermodynamic dG kcal/mol: {best.get('Worst_Thermo_dG_kcal')}\n")
            handle.write(f"Hard-filter failure reasons: {best.get('Hard_Filter_Failure_Reasons', '')}\n")
            handle.write(f"Soft relaxed reasons: {best.get('Soft_Filter_Relaxed_Reasons', '')}\n")
            handle.write(f"SNPs: {best['SNP_Positions_In_Amplicon']}\n")
        else:
            handle.write("No orderable primer pair passed the hard selectable filters. Inspect KOveri_primer_candidates.csv for failure reasons.\n")
