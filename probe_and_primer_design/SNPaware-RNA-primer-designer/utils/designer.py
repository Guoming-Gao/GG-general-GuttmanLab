"""High-level SNP-aware RNA primer design orchestration."""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

import pandas as pd
from rich.progress import track

from .io_utils import write_dataframe, write_fasta
from .models import CommonJunction, SnpRecord, TranscriptModel, TranscriptTemplate
from .post_filters import _annotate_row, annotate_post_filters
from .primer3_utils import run_primer3_cdna_pair_design
from .report import write_html_report
from .snp_utils import SnpDatabase
from .specificity import check_primer_specificity
from .tools import require_file, resolve_tool
from .transcripts import (
    build_transcript_template,
    cdna_interval_to_genomic_segments,
    common_junctions_for_template,
    format_segments,
    load_gene_transcripts,
    map_genomic_pos_to_cdna,
    select_template_transcript,
)


def design_snpaware_rna_primers(
    gene_name: str,
    output_dir: str | Path,
    config: dict,
    progress_callback=None,
) -> dict:
    """Run the full SNP-aware mature-RNA primer design pipeline."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gtf_path = require_file(config["gtf_path"], "GTF annotation")
    genome_fasta = require_file(config["genome_fasta"], "Genome FASTA")
    blast_db = require_file(config["blast_database"] + ".nsq", "BLAST database").removesuffix(".nsq")
    snp_vcf = require_file(config["snp_vcf"], "SNP VCF")
    samtools = resolve_tool(config, "samtools", "samtools")
    blastn = resolve_tool(config, "blastn", "blastn")
    primer3_core = resolve_tool(config, "primer3_core", "primer3_core")
    ntthal = resolve_tool(config, "ntthal", "ntthal")

    all_transcripts, curated_transcripts = load_gene_transcripts(
        gtf_path,
        gene_name,
        tuple(config["curated_transcript_prefixes"]),
    )
    if not all_transcripts:
        raise ValueError(f"No transcripts found for gene: {gene_name}")
    selected_transcript = select_template_transcript(curated_transcripts)
    template = build_transcript_template(selected_transcript, genome_fasta, samtools)
    common_junctions = common_junctions_for_template(curated_transcripts, selected_transcript, template)
    if not common_junctions:
        raise ValueError(f"No exon-exon junctions are common to curated multi-exon transcripts for {gene_name}.")

    snp_db = SnpDatabase(snp_vcf, config["vcf_b6_sample"], config["vcf_cast_sample"])
    try:
        snps, snp_cdna = _query_template_snps(template, snp_db)
    finally:
        snp_db.close()

    crossing_rows = _design_mode_candidates(
        mode="junction_crossing",
        gene_name=gene_name,
        selected_transcript=selected_transcript,
        template=template,
        common_junctions=common_junctions,
        snps=snps,
        snp_cdna=snp_cdna,
        primer3_core=primer3_core,
        config=config,
        ntthal=ntthal,
        progress_callback=progress_callback,
    )
    candidate_rows = _dedupe_candidates(crossing_rows)
    candidate_rows.sort(key=_pre_specificity_rank_key)
    if config.get("blast_specificity", True) and candidate_rows:
        _annotate_specificity(candidate_rows, blast_db, blastn, config)
    filter_info = annotate_post_filters(candidate_rows, config, ntthal)

    fallback_used = False
    if filter_info["selected_tier"] == "failed":
        fallback_used = True
        fallback_rows = _design_mode_candidates(
            mode="junction_spanning_fallback",
            gene_name=gene_name,
            selected_transcript=selected_transcript,
            template=template,
            common_junctions=common_junctions,
            snps=snps,
            snp_cdna=snp_cdna,
            primer3_core=primer3_core,
            config=config,
            ntthal=ntthal,
            progress_callback=progress_callback,
        )
        candidate_rows = _dedupe_candidates(candidate_rows + fallback_rows)
        candidate_rows.sort(key=_pre_specificity_rank_key)
        if config.get("blast_specificity", True) and candidate_rows:
            _annotate_specificity(candidate_rows, blast_db, blastn, config)
        filter_info = annotate_post_filters(candidate_rows, config, ntthal)

    candidate_rows.sort(key=_rank_key)
    top_rows = _select_top_rows(candidate_rows, filter_info["selected_tier"], int(config["top_n_output"]))

    transcript_rows = _transcript_rows(all_transcripts, selected_transcript)
    junction_rows = _junction_rows(common_junctions)
    write_dataframe(transcript_rows, output_dir / "SNPaware_RNA_transcripts.csv")
    write_dataframe(junction_rows, output_dir / "SNPaware_RNA_common_junctions.csv")
    write_dataframe(candidate_rows, output_dir / "SNPaware_RNA_primer_candidates.csv")
    _write_top_candidates(top_rows, candidate_rows, output_dir / "SNPaware_RNA_primer_top_candidates.csv")
    _write_final_fasta(top_rows, output_dir / "SNPaware_RNA_primers.fasta")
    _write_bed(top_rows, common_junctions, output_dir / "SNPaware_RNA_amplicons.bed")
    write_html_report(
        output_dir / "SNPaware_RNA_report.html",
        gene_name=gene_name,
        selected_transcript=selected_transcript,
        template=template,
        common_junctions=common_junctions,
        top_rows=top_rows,
        candidate_count=len(candidate_rows),
        config=config,
        filter_info=filter_info,
        fallback_used=fallback_used,
    )
    _write_summary(
        output_dir / "SNPaware_RNA_summary.txt",
        gene_name=gene_name,
        output_dir=output_dir,
        selected_transcript=selected_transcript,
        curated_transcripts=curated_transcripts,
        common_junctions=common_junctions,
        rows=candidate_rows,
        top_rows=top_rows,
        config=config,
        filter_info=filter_info,
        fallback_used=fallback_used,
    )

    return {
        "gene_name": gene_name,
        "selected_transcript": selected_transcript,
        "template": template,
        "common_junctions": common_junctions,
        "candidate_rows": candidate_rows,
        "top_rows": top_rows,
        "filter_info": filter_info,
        "fallback_used": fallback_used,
        "output_dir": output_dir,
    }


def design_snpaware_rna_primers_for_genes(
    gene_names: list[str],
    output_root: str | Path,
    config: dict,
    show_progress: bool = True,
) -> list[dict]:
    """Run the primer designer for multiple genes, writing each gene to its own folder."""
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    results: list[dict] = []
    iterator = track(gene_names, description="Designing RNA primers") if show_progress else gene_names
    for gene_name in iterator:
        gene_output = output_root / f"{gene_name}_SNPaware_RNA_primers"
        try:
            result = design_snpaware_rna_primers(gene_name, gene_output, config)
            results.append(
                {
                    "Gene": gene_name,
                    "Status": "success",
                    "Output_Dir": str(gene_output),
                    "Selected_Transcript": result["selected_transcript"].transcript_id,
                    "Common_Junctions": len(result["common_junctions"]),
                    "Candidate_Pairs": len(result["candidate_rows"]),
                    "Top_Pairs": len(result["top_rows"]),
                    "Selected_Filter_Tier": result["filter_info"]["selected_tier"],
                    "Fallback_Used": result["fallback_used"],
                    "Error": "",
                }
            )
        except Exception as exc:
            results.append(
                {
                    "Gene": gene_name,
                    "Status": "failed",
                    "Output_Dir": str(gene_output),
                    "Selected_Transcript": "",
                    "Common_Junctions": 0,
                    "Candidate_Pairs": 0,
                    "Top_Pairs": 0,
                    "Selected_Filter_Tier": "failed",
                    "Fallback_Used": "",
                    "Error": str(exc),
                }
            )
    pd.DataFrame(results).to_csv(output_root / "SNPaware_RNA_batch_summary.csv", index=False)
    return results


def _design_mode_candidates(
    mode: str,
    gene_name: str,
    selected_transcript: TranscriptModel,
    template: TranscriptTemplate,
    common_junctions: list[CommonJunction],
    snps: list[SnpRecord],
    snp_cdna: dict[tuple[str, int], int],
    primer3_core: str,
    config: dict,
    ntthal: str,
    progress_callback=None,
) -> list[dict]:
    rows: list[dict] = []
    min_snps = int(config["min_informative_snps"])
    goals = [min_snps] if mode == "junction_crossing" else list(range(min_snps, 0, -1))
    thermo_cache = {}
    attempt_total = _count_primer3_attempts(mode, goals, common_junctions, snps, snp_cdna, config)
    attempt_done = 0
    if progress_callback:
        progress_callback(
            "mode_start",
            mode=mode,
            total_attempts=attempt_total,
            message=f"{mode}: {attempt_total} Primer3 windows",
        )

    for snp_goal in goals:
        for junction in common_junctions:
            if mode == "junction_crossing":
                spans = [(junction.cdna_left, junction.cdna_right, [])]
            else:
                spans = _target_spans(junction, snps, snp_cdna, snp_goal, config["max_target_spans_per_flank"])
            if not spans:
                continue
            for flank in config["flank_steps"]:
                for span_index, (span_start, span_end, target_snps) in enumerate(spans):
                    attempt_done += 1
                    if progress_callback:
                        progress_callback(
                            "attempt_start",
                            mode=mode,
                            attempt_index=attempt_done,
                            total_attempts=attempt_total,
                            message=f"{mode} {attempt_done}/{attempt_total}: cDNA{junction.cdna_left}|{junction.cdna_right}, flank {flank}",
                        )
                    if span_end - span_start + 1 > int(config["max_amplicon_size"]):
                        if progress_callback:
                            progress_callback(
                                "attempt_done",
                                mode=mode,
                                attempt_index=attempt_done,
                                total_attempts=attempt_total,
                            )
                        continue
                    template_start = max(1, span_start - int(flank))
                    template_end = min(template.length, span_end + int(flank))
                    sequence = template.sequence[template_start - 1 : template_end]
                    sequence_id = (
                        f"{gene_name}_{selected_transcript.transcript_id}_{mode}_"
                        f"j{junction.cdna_left}_flank{flank}_span{span_index}"
                    )
                    try:
                        pairs = run_primer3_cdna_pair_design(
                            sequence_id=sequence_id,
                            template=sequence,
                            target_start_cdna=span_start,
                            target_end_cdna=span_end,
                            template_start_cdna=template_start,
                            excluded_positions_cdna=list(snp_cdna.values()),
                            config=config,
                            primer3_core=primer3_core,
                        )
                    except Exception as exc:
                        rows.append(
                            _failure_row(
                                gene_name,
                                selected_transcript,
                                mode,
                                junction,
                                flank,
                                snp_goal,
                                span_start,
                                span_end,
                                str(exc),
                            )
                        )
                        if progress_callback:
                            progress_callback(
                                "attempt_done",
                                mode=mode,
                                attempt_index=attempt_done,
                                total_attempts=attempt_total,
                            )
                        continue

                    for pair in pairs:
                        row = _score_pair(
                            pair=pair,
                            mode=mode,
                            gene_name=gene_name,
                            selected_transcript=selected_transcript,
                            template=template,
                            common_junctions=common_junctions,
                            target_junction=junction,
                            template_start=template_start,
                            template_end=template_end,
                            flank=flank,
                            snp_goal=snp_goal,
                            all_snps=snps,
                            snp_cdna=snp_cdna,
                            target_span_snps=target_snps,
                            ideal_amplicon_max=int(config["ideal_amplicon_max"]),
                            config=config,
                        )
                        if not row["Junction_Target_Pass"]:
                            continue
                        row["Failure_Reason"] = "primer_overlaps_snp" if row["Primer_Overlaps_SNP"] else ""
                        _annotate_row(row, config, ntthal, thermo_cache)
                        rows.append(row)
                    if progress_callback:
                        progress_callback(
                            "attempt_done",
                            mode=mode,
                            attempt_index=attempt_done,
                            total_attempts=attempt_total,
                        )

                if any(
                    row.get("Hard_Filter_Pass") is True
                    and row.get("SNP_Goal_For_Attempt") == snp_goal
                    and row.get("Design_Mode") == mode
                    for row in rows
                ):
                    if progress_callback:
                        progress_callback(
                            "mode_done",
                            mode=mode,
                            attempt_index=attempt_done,
                            total_attempts=attempt_total,
                        )
                    return rows
    if progress_callback:
        progress_callback(
            "mode_done",
            mode=mode,
            attempt_index=attempt_done,
            total_attempts=attempt_total,
        )
    return rows


def _count_primer3_attempts(
    mode: str,
    goals: list[int],
    common_junctions: list[CommonJunction],
    snps: list[SnpRecord],
    snp_cdna: dict[tuple[str, int], int],
    config: dict,
) -> int:
    total = 0
    for snp_goal in goals:
        for junction in common_junctions:
            if mode == "junction_crossing":
                spans = [(junction.cdna_left, junction.cdna_right, [])]
            else:
                spans = _target_spans(junction, snps, snp_cdna, snp_goal, config["max_target_spans_per_flank"])
            total += len(spans) * len(config["flank_steps"])
    return total


def _query_template_snps(template: TranscriptTemplate, snp_db: SnpDatabase) -> tuple[list[SnpRecord], dict[tuple[str, int], int]]:
    by_key: dict[tuple[str, int], SnpRecord] = {}
    cdna_by_key: dict[tuple[str, int], int] = {}
    for segment in template.segments:
        for snp in snp_db.query_informative(segment.exon.chrom, segment.exon.start, segment.exon.end):
            cdna_pos = map_genomic_pos_to_cdna(template, snp.pos)
            if cdna_pos is None:
                continue
            key = (snp.chrom, snp.pos)
            by_key[key] = snp
            cdna_by_key[key] = cdna_pos
    snps = sorted(by_key.values(), key=lambda snp: cdna_by_key[(snp.chrom, snp.pos)])
    return snps, cdna_by_key


def _target_spans(
    junction: CommonJunction,
    snps: list[SnpRecord],
    snp_cdna: dict[tuple[str, int], int],
    snp_goal: int,
    max_spans: int,
) -> list[tuple[int, int, list[SnpRecord]]]:
    if snp_goal <= 0:
        return [(junction.cdna_left, junction.cdna_right, [])]
    sorted_snps = sorted(snps, key=lambda snp: snp_cdna[(snp.chrom, snp.pos)])
    spans: list[tuple[int, int, list[SnpRecord]]] = []
    for i in range(0, len(sorted_snps) - snp_goal + 1):
        subset = sorted_snps[i : i + snp_goal]
        subset_positions = [snp_cdna[(snp.chrom, snp.pos)] for snp in subset]
        span_start = min(junction.cdna_left, *subset_positions)
        span_end = max(junction.cdna_right, *subset_positions)
        spans.append((span_start, span_end, subset))

    for i in range(len(sorted_snps)):
        for j in range(i + snp_goal, min(len(sorted_snps), i + snp_goal + 4) + 1):
            subset = sorted_snps[i:j]
            if len(subset) < snp_goal:
                continue
            subset_positions = [snp_cdna[(snp.chrom, snp.pos)] for snp in subset]
            span_start = min(junction.cdna_left, *subset_positions)
            span_end = max(junction.cdna_right, *subset_positions)
            spans.append((span_start, span_end, subset))

    unique = OrderedDict()
    for span_start, span_end, subset in spans:
        unique[(span_start, span_end)] = (span_start, span_end, subset)
    ranked = sorted(
        unique.values(),
        key=lambda item: (
            item[1] - item[0] + 1,
            -len(item[2]),
            abs((item[0] + item[1]) / 2 - (junction.cdna_left + junction.cdna_right) / 2),
        ),
    )
    return ranked[:max_spans]


def _score_pair(
    pair: dict,
    mode: str,
    gene_name: str,
    selected_transcript: TranscriptModel,
    template: TranscriptTemplate,
    common_junctions: list[CommonJunction],
    target_junction: CommonJunction,
    template_start: int,
    template_end: int,
    flank: int,
    snp_goal: int,
    all_snps: list[SnpRecord],
    snp_cdna: dict[tuple[str, int], int],
    target_span_snps: list[SnpRecord],
    ideal_amplicon_max: int,
    config: dict,
) -> dict:
    amp_start = pair["amplicon_cdna_start"]
    amp_end = pair["amplicon_cdna_end"]
    left_start = pair["left_cdna_start"]
    left_end = pair["left_cdna_end"]
    right_start = pair["right_cdna_start"]
    right_end = pair["right_cdna_end"]
    snps_in_amplicon = [s for s in all_snps if amp_start <= snp_cdna[(s.chrom, s.pos)] <= amp_end]
    left_overlaps = [s for s in all_snps if left_start <= snp_cdna[(s.chrom, s.pos)] <= left_end]
    right_overlaps = [s for s in all_snps if right_start <= snp_cdna[(s.chrom, s.pos)] <= right_end]
    left_crosses = _primer_crosses_junction(
        left_start,
        left_end,
        "left",
        common_junctions,
        config,
    )
    right_crosses = _primer_crosses_junction(
        right_start,
        right_end,
        "right",
        common_junctions,
        config,
    )
    left_target_crosses = _primer_crosses_junction(
        left_start,
        left_end,
        "left",
        [target_junction],
        config,
    )
    right_target_crosses = _primer_crosses_junction(
        right_start,
        right_end,
        "right",
        [target_junction],
        config,
    )
    spanning_junctions = _amplicon_spanning_junctions(amp_start, amp_end, common_junctions)
    target_pass = bool(left_target_crosses or right_target_crosses) if mode == "junction_crossing" else bool(spanning_junctions)
    primary_crossing = left_target_crosses[0] if left_target_crosses else right_target_crosses[0] if right_target_crosses else None
    crossing_side = "left" if left_target_crosses else "right" if right_target_crosses else ""
    left_arm_lengths = _junction_arm_lengths(left_start, left_end, "left", primary_crossing) if crossing_side == "left" else {}
    right_arm_lengths = _junction_arm_lengths(right_start, right_end, "right", primary_crossing) if crossing_side == "right" else {}
    tm_delta = abs(pair["left_tm"] - pair["right_tm"])
    ideal_penalty = max(0, pair["amplicon_size"] - ideal_amplicon_max)
    left_segments = cdna_interval_to_genomic_segments(template, left_start, left_end)
    right_segments = cdna_interval_to_genomic_segments(template, right_start, right_end)
    amp_segments = cdna_interval_to_genomic_segments(template, amp_start, amp_end)
    genomic_starts = [segment[1] for segment in amp_segments]
    genomic_ends = [segment[2] for segment in amp_segments]
    genomic_span_start = min(genomic_starts) if genomic_starts else ""
    genomic_span_end = max(genomic_ends) if genomic_ends else ""
    genomic_span = genomic_span_end - genomic_span_start + 1 if genomic_starts else ""
    return {
        "Gene": gene_name,
        "Transcript_ID": selected_transcript.transcript_id,
        "Chromosome": template.chrom,
        "Strand": template.strand,
        "Design_Mode": mode,
        "Target_Junction_cDNA_Left": target_junction.cdna_left,
        "Target_Junction_cDNA_Right": target_junction.cdna_right,
        "Target_Junction_Genomic": _format_junction(target_junction),
        "Primary_Crossing_Junction": _format_junction(primary_crossing) if primary_crossing else "",
        "Crossing_Primer_Side": crossing_side,
        "Crossing_Junctions": _format_crossings(left_crosses + right_crosses),
        "Spanned_Junctions": _format_crossings(spanning_junctions),
        "Design_Window_cDNA_Start": template_start,
        "Design_Window_cDNA_End": template_end,
        "Flank_Used": flank,
        "SNP_Goal_For_Attempt": snp_goal,
        "Target_SNPs_For_Attempt": _format_snps(target_span_snps, snp_cdna),
        "Amplicon_cDNA_Start": amp_start,
        "Amplicon_cDNA_End": amp_end,
        "Amplicon_cDNA_Size": pair["amplicon_size"],
        "Amplicon_Size": pair["amplicon_size"],
        "Amplicon_Genomic_Segments": format_segments(amp_segments),
        "Amplicon_Genomic_Span_Start": genomic_span_start,
        "Amplicon_Genomic_Span_End": genomic_span_end,
        "Amplicon_Genomic_Span": genomic_span,
        "Junction_Target_Pass": target_pass,
        "Covers_All_Guides": target_pass,
        "SNP_Count_In_Amplicon": len(snps_in_amplicon),
        "SNP_Positions_In_Amplicon": _format_snps(snps_in_amplicon, snp_cdna),
        "Left_Primer_Seq": pair["left_sequence"],
        "Left_Primer_cDNA_Start": left_start,
        "Left_Primer_cDNA_End": left_end,
        "Left_Primer_Genomic_Segments": format_segments(left_segments),
        "Left_Primer_Start": min((segment[1] for segment in left_segments), default=""),
        "Left_Primer_End": max((segment[2] for segment in left_segments), default=""),
        "Left_Primer_Crosses_Junction": bool(left_crosses),
        "Left_Primer_Crosses_Target_Junction": bool(left_target_crosses),
        "Left_Primer_Junction_5p_Arm_Bases": left_arm_lengths.get("five_prime", ""),
        "Left_Primer_Junction_3p_Arm_Bases": left_arm_lengths.get("three_prime", ""),
        "Left_Primer_Junction_Upstream_Arm_Bases": left_arm_lengths.get("upstream", ""),
        "Left_Primer_Junction_Downstream_Arm_Bases": left_arm_lengths.get("downstream", ""),
        "Left_Primer_Tm": pair["left_tm"],
        "Left_Primer_GC": pair["left_gc"],
        "Right_Primer_Seq": pair["right_sequence"],
        "Right_Primer_cDNA_Start": right_start,
        "Right_Primer_cDNA_End": right_end,
        "Right_Primer_Genomic_Segments": format_segments(right_segments),
        "Right_Primer_Start": min((segment[1] for segment in right_segments), default=""),
        "Right_Primer_End": max((segment[2] for segment in right_segments), default=""),
        "Right_Primer_Crosses_Junction": bool(right_crosses),
        "Right_Primer_Crosses_Target_Junction": bool(right_target_crosses),
        "Right_Primer_Junction_5p_Arm_Bases": right_arm_lengths.get("five_prime", ""),
        "Right_Primer_Junction_3p_Arm_Bases": right_arm_lengths.get("three_prime", ""),
        "Right_Primer_Junction_Upstream_Arm_Bases": right_arm_lengths.get("upstream", ""),
        "Right_Primer_Junction_Downstream_Arm_Bases": right_arm_lengths.get("downstream", ""),
        "Right_Primer_Tm": pair["right_tm"],
        "Right_Primer_GC": pair["right_gc"],
        "Primer_Tm_Delta": tm_delta,
        "Primer3_Pair_Penalty": pair["pair_penalty"],
        "Primer_Overlaps_SNP": bool(left_overlaps or right_overlaps),
        "Left_Primer_SNP_Overlaps": _format_snps(left_overlaps, snp_cdna),
        "Right_Primer_SNP_Overlaps": _format_snps(right_overlaps, snp_cdna),
        "Ideal_Size_Overage": ideal_penalty,
        "Left_Primer_PerfectOrNear_Hits": "",
        "Left_Primer_Noise_Hits": "",
        "Left_Primer_Specificity_Pass": "",
        "Left_Primer_Relaxed_Specificity_Pass": "",
        "Right_Primer_PerfectOrNear_Hits": "",
        "Right_Primer_Noise_Hits": "",
        "Right_Primer_Specificity_Pass": "",
        "Right_Primer_Relaxed_Specificity_Pass": "",
        "Primer_Pair_Specificity_Pass": "",
        "Primer_Pair_Relaxed_Specificity_Pass": "",
    }


def _primer_crosses_junction(
    start: int,
    end: int,
    primer_side: str,
    junctions: list[CommonJunction],
    config: dict,
) -> list[CommonJunction]:
    min_each_side = int(config["min_junction_overlap_bases"])
    min_3p = int(config["min_junction_3p_anchor_bases"])
    crossed = []
    for junction in junctions:
        if not (start <= junction.cdna_left and end >= junction.cdna_right):
            continue
        left_bases = junction.cdna_left - start + 1
        right_bases = end - junction.cdna_right + 1
        if left_bases < min_each_side or right_bases < min_each_side:
            continue
        three_prime_bases = right_bases if primer_side == "left" else left_bases
        if three_prime_bases < min_3p:
            continue
        crossed.append(junction)
    return crossed


def _junction_arm_lengths(
    start: int,
    end: int,
    primer_side: str,
    junction: CommonJunction | None,
) -> dict[str, int]:
    if junction is None:
        return {}
    upstream = junction.cdna_left - start + 1
    downstream = end - junction.cdna_right + 1
    if primer_side == "left":
        return {
            "upstream": upstream,
            "downstream": downstream,
            "five_prime": upstream,
            "three_prime": downstream,
        }
    return {
        "upstream": upstream,
        "downstream": downstream,
        "five_prime": downstream,
        "three_prime": upstream,
    }


def _amplicon_spanning_junctions(start: int, end: int, junctions: list[CommonJunction]) -> list[CommonJunction]:
    return [junction for junction in junctions if start <= junction.cdna_left and end >= junction.cdna_right]


def _failure_row(
    gene_name: str,
    selected_transcript: TranscriptModel,
    mode: str,
    junction: CommonJunction,
    flank: int,
    snp_goal: int,
    span_start: int,
    span_end: int,
    reason: str,
) -> dict:
    return {
        "Gene": gene_name,
        "Transcript_ID": selected_transcript.transcript_id,
        "Chromosome": selected_transcript.chrom,
        "Strand": selected_transcript.strand,
        "Design_Mode": mode,
        "Target_Junction_cDNA_Left": junction.cdna_left,
        "Target_Junction_cDNA_Right": junction.cdna_right,
        "Target_Junction_Genomic": _format_junction(junction),
        "Flank_Used": flank,
        "SNP_Goal_For_Attempt": snp_goal,
        "Amplicon_Size": "",
        "Junction_Target_Pass": False,
        "SNP_Count_In_Amplicon": 0,
        "Left_Primer_Seq": "",
        "Right_Primer_Seq": "",
        "Primer_Overlaps_SNP": "",
        "Failure_Reason": f"primer3_failed_for_target_span_{span_start}_{span_end}: {reason}",
    }


def _format_snps(snps: list[SnpRecord], snp_cdna: dict[tuple[str, int], int]) -> str:
    return ";".join(
        f"{s.chrom}:{s.pos}:cDNA{snp_cdna[(s.chrom, s.pos)]}:{s.ref}>{s.alt}:{s.genotype}"
        for s in snps
    )


def _format_junction(junction: CommonJunction) -> str:
    return (
        f"{junction.chrom}:{junction.upstream_exon_end_genomic}|"
        f"{junction.downstream_exon_start_genomic}:{junction.strand}:"
        f"cDNA{junction.cdna_left}|{junction.cdna_right}"
    )


def _format_crossings(junctions: list[CommonJunction]) -> str:
    return ";".join(_format_junction(junction) for junction in junctions)


def _dedupe_candidates(rows: list[dict]) -> list[dict]:
    seen = set()
    deduped = []
    for row in rows:
        key = (
            row.get("Left_Primer_Seq"),
            row.get("Right_Primer_Seq"),
            row.get("Amplicon_cDNA_Start"),
            row.get("Amplicon_cDNA_End"),
            row.get("Design_Mode"),
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
    mode_rank = 0 if row.get("Design_Mode") == "junction_crossing" else 1
    return (
        tier_rank,
        mode_rank,
        not row.get("Junction_Target_Pass", False),
        row.get("Primer_Overlaps_SNP", True),
        -int(row.get("SNP_Count_In_Amplicon") or 0),
        not specificity_pass,
        int(row.get("Ideal_Size_Overage") or 0),
        int(row.get("Amplicon_Size") or 999999),
        float(row.get("Primer_Tm_Delta") or 999),
        float(row.get("Primer3_Pair_Penalty") or 999),
    )


def _pre_specificity_rank_key(row: dict) -> tuple:
    mode_rank = 0 if row.get("Design_Mode") == "junction_crossing" else 1
    hard_pass = row.get("Hard_Filter_Pass") is True
    return (
        not hard_pass,
        mode_rank,
        not row.get("Junction_Target_Pass", False),
        row.get("Primer_Overlaps_SNP", True),
        -int(row.get("SNP_Count_In_Amplicon") or 0),
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
    for row in rows[: int(config.get("specificity_candidates", len(rows)))]:
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
        left_pass, left_relaxed = _specificity_passes(
            left_spec.perfect_or_near_hits,
            left_spec.noise_hits,
            bool(row.get("Left_Primer_Crosses_Junction")),
        )
        right_pass, right_relaxed = _specificity_passes(
            right_spec.perfect_or_near_hits,
            right_spec.noise_hits,
            bool(row.get("Right_Primer_Crosses_Junction")),
        )
        row.update(
            {
                "Left_Primer_PerfectOrNear_Hits": left_spec.perfect_or_near_hits,
                "Left_Primer_Noise_Hits": left_spec.noise_hits,
                "Left_Primer_Specificity_Pass": left_pass,
                "Left_Primer_Relaxed_Specificity_Pass": left_relaxed,
                "Right_Primer_PerfectOrNear_Hits": right_spec.perfect_or_near_hits,
                "Right_Primer_Noise_Hits": right_spec.noise_hits,
                "Right_Primer_Specificity_Pass": right_pass,
                "Right_Primer_Relaxed_Specificity_Pass": right_relaxed,
                "Primer_Pair_Specificity_Pass": left_pass and right_pass,
                "Primer_Pair_Relaxed_Specificity_Pass": left_relaxed and right_relaxed,
            }
        )


def _specificity_passes(perfect_or_near_hits: int, noise_hits: int, crosses_junction: bool) -> tuple[bool, bool]:
    if crosses_junction:
        return perfect_or_near_hits == 0 and noise_hits == 0, perfect_or_near_hits == 0
    return perfect_or_near_hits == 1 and noise_hits == 0, perfect_or_near_hits >= 1


def _write_final_fasta(rows: list[dict], path: Path) -> None:
    records: list[tuple[str, str]] = []
    for idx, row in enumerate(rows, start=1):
        records.append((
            f"SNPaware_RNA_pair{idx}_LEFT_{row['Gene']}_{row['Transcript_ID']}_"
            f"cDNA{row['Left_Primer_cDNA_Start']}-{row['Left_Primer_cDNA_End']}_"
            f"genomic={row['Left_Primer_Genomic_Segments']}",
            row["Left_Primer_Seq"],
        ))
        records.append((
            f"SNPaware_RNA_pair{idx}_RIGHT_{row['Gene']}_{row['Transcript_ID']}_"
            f"cDNA{row['Right_Primer_cDNA_Start']}-{row['Right_Primer_cDNA_End']}_"
            f"genomic={row['Right_Primer_Genomic_Segments']}",
            row["Right_Primer_Seq"],
        ))
    write_fasta(records, path)


def _write_bed(rows: list[dict], common_junctions: list[CommonJunction], path: Path) -> None:
    with path.open("w") as handle:
        handle.write("#chrom\tstart0\tend1\tname\tscore\tstrand\tfeature\n")
        for idx, junction in enumerate(common_junctions, start=1):
            left = min(junction.upstream_exon_end_genomic, junction.downstream_exon_start_genomic)
            right = max(junction.upstream_exon_end_genomic, junction.downstream_exon_start_genomic)
            handle.write(f"{junction.chrom}\t{left - 1}\t{right}\tcommon_junction_{idx}\t0\t{junction.strand}\tcommon_junction\n")
        for idx, row in enumerate(rows, start=1):
            _write_segments_bed(handle, row["Amplicon_Genomic_Segments"], f"SNPaware_RNA_pair{idx}_amplicon", row["SNP_Count_In_Amplicon"], "amplicon")
            _write_segments_bed(handle, row["Left_Primer_Genomic_Segments"], f"SNPaware_RNA_pair{idx}_left_primer", 0, "primer")
            _write_segments_bed(handle, row["Right_Primer_Genomic_Segments"], f"SNPaware_RNA_pair{idx}_right_primer", 0, "primer")


def _write_segments_bed(handle, segments_text: str, name: str, score: int, feature: str) -> None:
    for segment in str(segments_text).split(";"):
        if not segment:
            continue
        coord, strand = segment.rsplit(":", 1)
        chrom, span = coord.split(":", 1)
        start, end = span.split("-", 1)
        handle.write(f"{chrom}\t{int(start) - 1}\t{end}\t{name}\t{score}\t{strand}\t{feature}\n")


def _write_summary(
    path: Path,
    gene_name: str,
    output_dir: Path,
    selected_transcript: TranscriptModel,
    curated_transcripts: list[TranscriptModel],
    common_junctions: list[CommonJunction],
    rows: list[dict],
    top_rows: list[dict],
    config: dict,
    filter_info: dict,
    fallback_used: bool,
) -> None:
    best = top_rows[0] if top_rows else None
    with path.open("w") as handle:
        handle.write("SNP-aware RNA primer design summary\n")
        handle.write("===================================\n\n")
        handle.write(f"Gene: {gene_name}\n")
        handle.write(f"Output directory: {output_dir}\n")
        handle.write(f"Selected transcript: {selected_transcript.transcript_id}\n")
        handle.write(f"Selected transcript strand: {selected_transcript.strand}\n")
        handle.write(f"Curated multi-exon transcripts: {len(curated_transcripts)}\n")
        handle.write(f"Common junctions: {len(common_junctions)}\n")
        handle.write(f"Minimum informative SNP goal: {config['min_informative_snps']}\n")
        handle.write(f"Fallback junction-spanning search used: {fallback_used}\n")
        handle.write(f"Candidate primer pairs: {len(rows)}\n")
        handle.write(f"Top primer pairs reported: {len(top_rows)}\n\n")
        handle.write("Post-Primer3 filters\n")
        handle.write("--------------------\n")
        handle.write(f"Selected filter tier: {filter_info['selected_tier']}\n")
        handle.write(f"Tier counts: {filter_info['tier_counts']}\n")
        handle.write(f"Hard-filter passing candidates: {filter_info['hard_pass_count']}\n")
        handle.write(f"Hard-filter failing candidates: {filter_info['hard_fail_count']}\n\n")
        if best:
            handle.write("Best ranked primer pair\n")
            handle.write("-----------------------\n")
            handle.write(f"Design mode: {best['Design_Mode']}\n")
            handle.write(f"Left primer:  {best['Left_Primer_Seq']} ({best['Left_Primer_Genomic_Segments']})\n")
            handle.write(f"Right primer: {best['Right_Primer_Seq']} ({best['Right_Primer_Genomic_Segments']})\n")
            handle.write(f"cDNA amplicon: {best['Amplicon_cDNA_Start']}-{best['Amplicon_cDNA_End']} ({best.get('Amplicon_cDNA_Size', best['Amplicon_Size'])} bp)\n")
            handle.write(f"Genomic amplicon segments: {best['Amplicon_Genomic_Segments']}\n")
            handle.write(f"Genomic span, not PCR product length: {best.get('Amplicon_Genomic_Span', '')} bp\n")
            handle.write(f"Target junction: {best.get('Target_Junction_Genomic', '')}\n")
            handle.write(f"Primary crossing junction: {best.get('Primary_Crossing_Junction', '')}\n")
            handle.write(f"Crossing primer side: {best.get('Crossing_Primer_Side', '')}\n")
            handle.write(
                "Crossing primer arm lengths: "
                f"left 5p/3p={best.get('Left_Primer_Junction_5p_Arm_Bases', '')}/{best.get('Left_Primer_Junction_3p_Arm_Bases', '')}; "
                f"right 5p/3p={best.get('Right_Primer_Junction_5p_Arm_Bases', '')}/{best.get('Right_Primer_Junction_3p_Arm_Bases', '')}\n"
            )
            handle.write(f"Crossing junctions: {best.get('Crossing_Junctions', '')}\n")
            handle.write(f"Spanned junctions: {best.get('Spanned_Junctions', '')}\n")
            handle.write(f"Informative SNPs in amplicon: {best['SNP_Count_In_Amplicon']}\n")
            handle.write(f"Primer overlaps informative SNP: {best['Primer_Overlaps_SNP']}\n")
            handle.write(f"Primer pair specificity pass: {best.get('Primer_Pair_Specificity_Pass')}\n")
            handle.write(f"Filter tier: {best.get('Filter_Tier')}\n")
            handle.write(f"Relaxed criteria used: {best.get('Relaxed_Criteria_Used', '')}\n")
            handle.write(f"Worst thermodynamic dG kcal/mol: {best.get('Worst_Thermo_dG_kcal')}\n")
            handle.write(f"Hard-filter failure reasons: {best.get('Hard_Filter_Failure_Reasons', '')}\n")
            handle.write(f"SNPs: {best['SNP_Positions_In_Amplicon']}\n")
        else:
            handle.write("No orderable primer pair passed the hard selectable filters. Inspect SNPaware_RNA_primer_candidates.csv.\n")


def _transcript_rows(all_transcripts: list[TranscriptModel], selected: TranscriptModel) -> list[dict]:
    return [
        {
            "Gene": transcript.gene_name,
            "Transcript_ID": transcript.transcript_id,
            "Chromosome": transcript.chrom,
            "Strand": transcript.strand,
            "Transcript_Start": transcript.start,
            "Transcript_End": transcript.end,
            "Exon_Count": transcript.exon_count,
            "Spliced_Length": transcript.spliced_length,
            "Curated_MultiExon_Used_For_Common_Set": transcript.is_curated,
            "Selected_Template_Transcript": transcript.transcript_id == selected.transcript_id,
            "Exclusion_Reason": transcript.exclusion_reason,
        }
        for transcript in all_transcripts
    ]


def _junction_rows(junctions: list[CommonJunction]) -> list[dict]:
    return [
        {
            "Chromosome": junction.chrom,
            "Strand": junction.strand,
            "Upstream_Exon_End_Genomic": junction.upstream_exon_end_genomic,
            "Downstream_Exon_Start_Genomic": junction.downstream_exon_start_genomic,
            "cDNA_Left": junction.cdna_left,
            "cDNA_Right": junction.cdna_right,
            "Upstream_Exon_Number": junction.upstream_exon_number,
            "Downstream_Exon_Number": junction.downstream_exon_number,
            "Curated_Transcript_Count": junction.transcript_count,
        }
        for junction in junctions
    ]
