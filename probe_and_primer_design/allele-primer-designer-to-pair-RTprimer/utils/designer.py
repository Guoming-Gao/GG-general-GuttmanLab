"""Design allele-aware gene-specific primers to pair with an RT/poly-dT handle."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn, TimeElapsedColumn
from rich.console import Console

from .coordinate_rules import primer_to_snp_distance_toward_polydt, snp_distance_to_rna_3p
from .io_utils import write_dataframe, write_fasta
from .models import SnpRecord, TranscriptModel, TranscriptTemplate
from .primer3_single import run_primer3_left_primer_design
from .report import write_consolidated_html_report, write_html_report
from .snp_utils import SnpDatabase
from .specificity import check_primer_specificity_batch
from .tools import require_file, resolve_tool
from .transcripts import (
    build_transcript_template,
    cdna_interval_to_genomic_segments,
    format_segments,
    load_gene_transcripts,
    map_genomic_pos_to_cdna,
    select_template_transcript,
)


HANDLE_SEQUENCES = {
    "2PBC": "CAGACGTGTGCTCTTCCGATCT",
    "2PUNI": "CCTACACGACGCTCTTCCGATCT",
}
console = Console()


def design_allele_primers(
    gene_name: str,
    output_dir: str | Path,
    config: dict,
) -> dict:
    """Run the full single-insert allele primer design workflow for one gene."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gtf_path = require_file(config["gtf_path"], "GTF annotation")
    genome_fasta = require_file(config["genome_fasta"], "Genome FASTA")
    snp_vcf = require_file(config["snp_vcf"], "SNP VCF")
    blast_db = ""
    if config.get("blast_specificity", True):
        blast_db = require_file(config["blast_database"] + ".nsq", "BLAST database").removesuffix(".nsq")
    samtools = resolve_tool(config, "samtools", "samtools")
    blastn = resolve_tool(config, "blastn", "blastn")
    primer3_core = resolve_tool(config, "primer3_core", "primer3_core")

    all_transcripts, curated_transcripts = load_gene_transcripts(
        gtf_path,
        gene_name,
        tuple(config["curated_transcript_prefixes"]),
    )
    if not all_transcripts:
        raise ValueError(f"No transcripts found for gene: {gene_name}")
    selected_transcript = select_template_transcript(curated_transcripts)
    template = build_transcript_template(selected_transcript, genome_fasta, samtools)
    snps, snp_cdna = _query_template_snps(template, config)

    snp_rows = _snp_rows(template, snps, snp_cdna)
    candidate_rows = _design_candidates_for_snps(
        gene_name=gene_name,
        selected_transcript=selected_transcript,
        template=template,
        snps=snps,
        snp_cdna=snp_cdna,
        primer3_core=primer3_core,
        blast_db=blast_db,
        blastn=blastn,
        config=config,
    )
    top_insert_rows = _select_top_insert_rows(candidate_rows, int(config["top_snps"]))
    top_oligo_rows = _expand_handle(top_insert_rows, config.get("handle", "2PUNI"))

    transcript_rows = _transcript_rows(all_transcripts, selected_transcript)
    write_dataframe(transcript_rows, output_dir / "allele_primer_transcripts.csv")
    write_dataframe(snp_rows, output_dir / "allele_primer_informative_snps.csv")
    write_dataframe(candidate_rows, output_dir / "allele_primer_candidate_inserts.csv")
    write_dataframe(top_oligo_rows, output_dir / "allele_primer_top_oligos.csv")
    _write_final_fasta(top_oligo_rows, output_dir / "allele_primer_oligos.fasta")
    write_html_report(
        output_dir / "allele_primer_report.html",
        gene_name=gene_name,
        selected_transcript=selected_transcript,
        template=template,
        snp_rows=snp_rows,
        candidate_rows=candidate_rows,
        top_insert_rows=top_insert_rows,
        top_oligo_rows=top_oligo_rows,
        config=config,
    )
    _write_summary(
        output_dir / "allele_primer_summary.txt",
        gene_name=gene_name,
        selected_transcript=selected_transcript,
        template=template,
        snp_rows=snp_rows,
        candidate_rows=candidate_rows,
        top_insert_rows=top_insert_rows,
        top_oligo_rows=top_oligo_rows,
        config=config,
    )

    return {
        "gene_name": gene_name,
        "selected_transcript": selected_transcript,
        "template": template,
        "snps": snps,
        "snp_cdna": snp_cdna,
        "candidate_rows": candidate_rows,
        "top_insert_rows": top_insert_rows,
        "top_oligo_rows": top_oligo_rows,
        "output_dir": output_dir,
    }


def design_allele_primers_for_genes(
    gene_names: list[str],
    output_root: str | Path,
    config: dict,
    show_progress: bool = True,
) -> list[dict]:
    """Run the workflow for multiple genes, writing one subfolder per gene."""
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    results = []
    consolidated_rows = []
    iterator = gene_names
    for gene_index, gene_name in enumerate(iterator, start=1):
        if show_progress:
            console.print(f"[cyan][Gene {gene_index}/{len(gene_names)}][/cyan] Starting {gene_name}")
        gene_output = output_root / f"{gene_name}_allele_primers"
        try:
            result = design_allele_primers(gene_name, gene_output, config)
            consolidated_rows.extend(result["top_oligo_rows"])
            if show_progress:
                console.print(f"[green][Gene {gene_index}/{len(gene_names)}][/green] Completed {gene_name}: {len(result['top_oligo_rows'])} top oligos")
            results.append(
                {
                    "Gene": gene_name,
                    "Status": "success",
                    "Output_Dir": str(gene_output),
                    "Selected_Transcript": result["selected_transcript"].transcript_id,
                    "Strand": result["template"].strand,
                    "Informative_SNPs": len(result["snps"]),
                    "Candidate_Inserts": len(result["candidate_rows"]),
                    "Top_Inserts": len(result["top_insert_rows"]),
                    "Top_Oligos": len(result["top_oligo_rows"]),
                    "Error": "",
                }
            )
        except Exception as exc:
            if show_progress:
                console.print(f"[red][Gene {gene_index}/{len(gene_names)}][/red] Failed {gene_name}: {exc}")
            results.append(
                {
                    "Gene": gene_name,
                    "Status": "failed",
                    "Output_Dir": str(gene_output),
                    "Selected_Transcript": "",
                    "Strand": "",
                    "Informative_SNPs": 0,
                    "Candidate_Inserts": 0,
                    "Top_Inserts": 0,
                    "Top_Oligos": 0,
                    "Error": str(exc),
                }
            )
    pd.DataFrame(results).to_csv(output_root / "allele_primer_batch_summary.csv", index=False)
    write_dataframe(consolidated_rows, output_root / "allele_primer_consolidated_top_oligos.csv")
    write_consolidated_html_report(
        output_root / "allele_primer_consolidated_report.html",
        batch_rows=results,
        top_oligo_rows=consolidated_rows,
        config=config,
    )
    return results


def _query_template_snps(template: TranscriptTemplate, config: dict) -> tuple[list[SnpRecord], dict[tuple[str, int], int]]:
    snp_db = SnpDatabase(config["snp_vcf"], config["vcf_b6_sample"], config["vcf_cast_sample"])
    by_key: dict[tuple[str, int], SnpRecord] = {}
    cdna_by_key: dict[tuple[str, int], int] = {}
    try:
        for segment in template.segments:
            for snp in snp_db.query_informative(segment.exon.chrom, segment.exon.start, segment.exon.end):
                cdna_pos = map_genomic_pos_to_cdna(template, snp.pos)
                if cdna_pos is None:
                    continue
                key = (snp.chrom, snp.pos)
                by_key[key] = snp
                cdna_by_key[key] = cdna_pos
    finally:
        snp_db.close()
    snps = sorted(by_key.values(), key=lambda snp: cdna_by_key[(snp.chrom, snp.pos)])
    return snps, cdna_by_key


def _design_candidates_for_snps(
    gene_name: str,
    selected_transcript: TranscriptModel,
    template: TranscriptTemplate,
    snps: list[SnpRecord],
    snp_cdna: dict[tuple[str, int], int],
    primer3_core: str,
    blast_db: str,
    blastn: str,
    config: dict,
) -> list[dict]:
    rows = []
    window = int(config["snp_window_toward_polydt"])
    informative_positions = list(snp_cdna.values())
    ranked_snps = sorted(snps, key=lambda snp: snp_distance_to_rna_3p(template.length, snp_cdna[(snp.chrom, snp.pos)]))
    accepted_snps = 0
    progress = Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        transient=False,
        disable=not config.get("show_progress", True),
    )
    task_id = progress.add_task(f"{gene_name}: SNP targets (0/{config['top_snps']} accepted)", total=len(ranked_snps))

    with progress:
      for snp_rank, snp in enumerate(ranked_snps, start=1):
        target_cdna = snp_cdna[(snp.chrom, snp.pos)]
        design_start = max(1, target_cdna - window)
        design_end = target_cdna - 1
        if design_end - design_start + 1 < int(config["primer_min_size"]):
            progress.advance(task_id)
            continue
        design_seq = template.sequence[design_start - 1 : design_end]
        try:
            primers = run_primer3_left_primer_design(
                sequence_id=f"{gene_name}_{selected_transcript.transcript_id}_snp{target_cdna}",
                template=design_seq,
                template_start_cdna=design_start,
                excluded_positions_cdna=informative_positions,
                config=config,
                primer3_core=primer3_core,
            )
        except Exception as exc:
            rows.append(_failure_row(gene_name, selected_transcript, template, snp, target_cdna, snp_rank, str(exc), config))
            progress.advance(task_id)
            continue

        snp_accepted = False
        scored_rows = [
            _score_primer(
                gene_name=gene_name,
                selected_transcript=selected_transcript,
                template=template,
                snp=snp,
                target_cdna=target_cdna,
                snp_rank=snp_rank,
                primer=primer,
                snps=snps,
                snp_cdna=snp_cdna,
                config=config,
            )
            for primer in sorted(primers, key=lambda item: (item["penalty"], abs(item["tm"] - 60.0)))
        ]
        viable_rows = [row for row in scored_rows if row["Filter_Pass"]]
        blast_limit = int(config.get("blast_candidates_per_snp", 5))
        evaluated_rows = viable_rows[:blast_limit]
        specificity_results = []
        if config.get("blast_specificity", True):
            specificity_results = check_primer_specificity_batch(
                [row["Primer_Insert_Seq"] for row in evaluated_rows], blast_db, blastn, config
            )
        result_by_sequence = {result.sequence: result for result in specificity_results}

        for row in scored_rows:
            if config.get("blast_specificity", True) and row["Filter_Pass"]:
                result = result_by_sequence.get(row["Primer_Insert_Seq"])
                if result is None:
                    row["Filter_Pass"] = False
                    row["Failure_Reasons"] = _append_reason(row["Failure_Reasons"], "specificity_not_evaluated_candidate_limit")
                    rows.append(row)
                    continue
                row.update(
                    {
                        "Primer_PerfectOrNear_Hits": result.perfect_or_near_hits,
                        "Primer_Noise_Hits": result.noise_hits,
                        "Primer_Specificity_Pass": result.pass_specificity,
                    }
                )
                row["Filter_Pass"] = row["Filter_Pass"] and result.pass_specificity
                if not result.pass_specificity:
                    row["Failure_Reasons"] = _append_reason(row["Failure_Reasons"], "specificity_failed")
            rows.append(row)
            if row["Filter_Pass"]:
                snp_accepted = True
                break
        if snp_accepted:
            accepted_snps += 1
        progress.update(
            task_id,
            advance=1,
            description=f"{gene_name}: SNP targets ({accepted_snps}/{config['top_snps']} accepted)",
        )
        if accepted_snps >= int(config["top_snps"]):
            break
    rows.sort(key=_rank_key)
    return rows


def _score_primer(
    gene_name: str,
    selected_transcript: TranscriptModel,
    template: TranscriptTemplate,
    snp: SnpRecord,
    target_cdna: int,
    snp_rank: int,
    primer: dict,
    snps: list[SnpRecord],
    snp_cdna: dict[tuple[str, int], int],
    config: dict,
) -> dict:
    primer_start = primer["cdna_start"]
    primer_end = primer["cdna_end"]
    window = int(config["snp_window_toward_polydt"])
    primer_to_snp = primer_to_snp_distance_toward_polydt(primer_end, target_cdna, window)
    overlaps = [s for s in snps if primer_start <= snp_cdna[(s.chrom, s.pos)] <= primer_end]
    downstream_snps = [
        s
        for s in snps
        if primer_to_snp_distance_toward_polydt(primer_end, snp_cdna[(s.chrom, s.pos)], window) is not None
    ]
    segments = cdna_interval_to_genomic_segments(template, primer_start, primer_end)
    reasons = []
    if primer_to_snp is None:
        reasons.append("target_snp_not_toward_polydt_within_window")
    if overlaps:
        reasons.append("primer_overlaps_informative_snp")
    if not _in_range(primer["length"], config["valid_primer_length_min"], config["valid_primer_length_max"]):
        reasons.append("length_outside_range")
    if not _in_range(primer["gc"], config["post_filter_gc_min"], config["post_filter_gc_max"]):
        reasons.append("gc_outside_range")
    if not _in_range(primer["tm"], config["post_filter_tm_min"], config["post_filter_tm_max"]):
        reasons.append("tm_outside_range")
    if not _has_3p_gc_clamp(primer["sequence"]):
        reasons.append("missing_3p_gc_clamp")
    if _has_homopolymer(primer["sequence"], int(config["max_homopolymer_run"]) + 1):
        reasons.append("homopolymer_run_ge_4")
    if "G" * (int(config["max_quad_g_run"]) + 1) in primer["sequence"]:
        reasons.append("quad_g_ge_4")

    filter_pass = not reasons
    return {
        "Gene": gene_name,
        "Transcript_ID": selected_transcript.transcript_id,
        "Chromosome": template.chrom,
        "Strand": template.strand,
        "Template_Length_cDNA": template.length,
        "SNP_Rank_From_RNA_3p": snp_rank,
        "SNP_cDNA": target_cdna,
        "SNP_Genomic": snp.pos,
        "SNP_B6_GT": snp.b6_gt,
        "SNP_CAST_GT": snp.cast_gt,
        "SNP_Distance_To_RNA_3p": snp_distance_to_rna_3p(template.length, target_cdna),
        "Primer_To_SNP_Distance_Toward_PolyDT": primer_to_snp if primer_to_snp is not None else "",
        "SNP_Window_Toward_PolyDT": window,
        "Primer_Insert_Seq": primer["sequence"],
        "Primer_cDNA_Start": primer_start,
        "Primer_cDNA_End": primer_end,
        "Primer_Genomic_Segments": format_segments(segments),
        "Primer_Length": primer["length"],
        "Primer_Tm": round(primer["tm"], 3) if not math.isnan(primer["tm"]) else "",
        "Primer_GC": round(primer["gc"], 3) if not math.isnan(primer["gc"]) else "",
        "Primer3_Penalty": primer["penalty"],
        "Primer3_SelfAny": primer["self_any"],
        "Primer3_SelfEnd": primer["self_end"],
        "Primer3_Hairpin": primer["hairpin"],
        "Primer_Overlaps_SNP": bool(overlaps),
        "Primer_SNP_Overlaps": _format_snps(overlaps, snp_cdna),
        "Downstream_SNPs_Within_Window": _format_snps(downstream_snps, snp_cdna),
        "Downstream_SNP_Count_Within_Window": len(downstream_snps),
        "Primer_PerfectOrNear_Hits": "",
        "Primer_Noise_Hits": "",
        "Primer_Specificity_Pass": "",
        "Filter_Pass": filter_pass,
        "Failure_Reasons": ";".join(reasons),
    }


def _select_top_insert_rows(candidate_rows: list[dict], top_snps: int) -> list[dict]:
    selected = []
    seen_snps = set()
    for row in sorted(candidate_rows, key=_rank_key):
        if row.get("Filter_Pass") is not True:
            continue
        snp_key = (row["Chromosome"], row["SNP_Genomic"])
        if snp_key in seen_snps:
            continue
        selected.append(row)
        seen_snps.add(snp_key)
        if len(selected) >= top_snps:
            break
    return selected


def _expand_handle(top_insert_rows: list[dict], handle_name: str) -> list[dict]:
    if handle_name not in HANDLE_SEQUENCES:
        raise ValueError(f"Unknown handle: {handle_name}")
    handle_seq = HANDLE_SEQUENCES[handle_name]
    rows = []
    for insert_rank, row in enumerate(top_insert_rows, start=1):
        expanded = dict(row)
        expanded.update(
            {
                "Insert_Rank": insert_rank,
                "Handle_Name": handle_name,
                "Handle_Sequence": handle_seq,
                "Final_Order_Sequence": handle_seq + row["Primer_Insert_Seq"],
            }
        )
        rows.append(expanded)
    return rows


def _snp_rows(template: TranscriptTemplate, snps: list[SnpRecord], snp_cdna: dict[tuple[str, int], int]) -> list[dict]:
    rows = []
    for snp in sorted(snps, key=lambda s: snp_distance_to_rna_3p(template.length, snp_cdna[(s.chrom, s.pos)])):
        cdna = snp_cdna[(snp.chrom, snp.pos)]
        rows.append(
            {
                "Gene": template.gene_name,
                "Transcript_ID": template.transcript_id,
                "Chromosome": snp.chrom,
                "Strand": template.strand,
                "SNP_cDNA": cdna,
                "SNP_Genomic": snp.pos,
                "SNP_Ref": snp.ref,
                "SNP_Alt": snp.alt,
                "SNP_B6_GT": snp.b6_gt,
                "SNP_CAST_GT": snp.cast_gt,
                "SNP_Distance_To_RNA_3p": snp_distance_to_rna_3p(template.length, cdna),
            }
        )
    return rows


def _transcript_rows(all_transcripts: list[TranscriptModel], selected: TranscriptModel) -> list[dict]:
    rows = []
    for tx in all_transcripts:
        rows.append(
            {
                "Gene": tx.gene_name,
                "Transcript_ID": tx.transcript_id,
                "Chromosome": tx.chrom,
                "Strand": tx.strand,
                "Start": tx.start,
                "End": tx.end,
                "Exon_Count": tx.exon_count,
                "Spliced_Length": tx.spliced_length,
                "Is_Curated_Designable": tx.is_curated,
                "Selected_Template": tx.transcript_id == selected.transcript_id,
                "Exclusion_Reason": tx.exclusion_reason,
            }
        )
    return rows


def _failure_row(
    gene_name: str,
    selected_transcript: TranscriptModel,
    template: TranscriptTemplate,
    snp: SnpRecord,
    target_cdna: int,
    snp_rank: int,
    reason: str,
    config: dict,
) -> dict:
    return {
        "Gene": gene_name,
        "Transcript_ID": selected_transcript.transcript_id,
        "Chromosome": template.chrom,
        "Strand": template.strand,
        "Template_Length_cDNA": template.length,
        "SNP_Rank_From_RNA_3p": snp_rank,
        "SNP_cDNA": target_cdna,
        "SNP_Genomic": snp.pos,
        "SNP_B6_GT": snp.b6_gt,
        "SNP_CAST_GT": snp.cast_gt,
        "SNP_Distance_To_RNA_3p": snp_distance_to_rna_3p(template.length, target_cdna),
        "Primer_To_SNP_Distance_Toward_PolyDT": "",
        "SNP_Window_Toward_PolyDT": config["snp_window_toward_polydt"],
        "Primer_Insert_Seq": "",
        "Primer_cDNA_Start": "",
        "Primer_cDNA_End": "",
        "Primer_Genomic_Segments": "",
        "Primer_Length": "",
        "Primer_Tm": "",
        "Primer_GC": "",
        "Primer3_Penalty": "",
        "Primer3_SelfAny": "",
        "Primer3_SelfEnd": "",
        "Primer3_Hairpin": "",
        "Primer_Overlaps_SNP": "",
        "Primer_SNP_Overlaps": "",
        "Downstream_SNPs_Within_Window": "",
        "Downstream_SNP_Count_Within_Window": 0,
        "Primer_PerfectOrNear_Hits": "",
        "Primer_Noise_Hits": "",
        "Primer_Specificity_Pass": "",
        "Filter_Pass": False,
        "Failure_Reasons": f"primer3_error:{reason}",
    }


def _write_final_fasta(rows: list[dict], path: Path) -> None:
    records = []
    for row in rows:
        name = (
            f"{row['Gene']}_{row['Transcript_ID']}_snp{row['SNP_cDNA']}_"
            f"{row['Handle_Name']}_rank{row['Insert_Rank']}"
        )
        records.append((name, row["Final_Order_Sequence"]))
    write_fasta(records, path)


def _write_summary(
    path: Path,
    gene_name: str,
    selected_transcript: TranscriptModel,
    template: TranscriptTemplate,
    snp_rows: list[dict],
    candidate_rows: list[dict],
    top_insert_rows: list[dict],
    top_oligo_rows: list[dict],
    config: dict,
) -> None:
    lines = [
        f"Allele primer design summary for {gene_name}",
        f"Selected transcript: {selected_transcript.transcript_id}",
        f"Coordinates: {template.chrom}:{selected_transcript.start}-{selected_transcript.end} ({template.strand})",
        f"Spliced RNA/cDNA length: {template.length}",
        f"Informative exonic B6/Cast SNPs: {len(snp_rows)}",
        f"Candidate inserts: {len(candidate_rows)}",
        f"Passing candidate inserts: {sum(row.get('Filter_Pass') is True for row in candidate_rows)}",
        f"Top inserts reported: {len(top_insert_rows)}",
        f"Top oligos reported: {len(top_oligo_rows)}",
        f"SNP window toward poly-dT: {config['snp_window_toward_polydt']} nt",
        f"Selected handle: {config.get('handle', '2PUNI')}",
        "",
    ]
    if top_insert_rows:
        lines.append("Top inserts:")
        for idx, row in enumerate(top_insert_rows, start=1):
            lines.append(
                f"{idx}. SNP cDNA {row['SNP_cDNA']} genomic {row['SNP_Genomic']} "
                f"distance_to_RNA_3p {row['SNP_Distance_To_RNA_3p']} nt; "
                f"primer {row['Primer_Insert_Seq']} cDNA {row['Primer_cDNA_Start']}-{row['Primer_cDNA_End']} "
                f"primer_to_snp {row['Primer_To_SNP_Distance_Toward_PolyDT']} nt"
            )
    else:
        lines.append("No passing inserts were found.")
    path.write_text("\n".join(lines) + "\n")


def _format_snps(snps: list[SnpRecord], snp_cdna: dict[tuple[str, int], int]) -> str:
    return ";".join(
        f"{snp.chrom}:{snp.pos}:cDNA{snp_cdna[(snp.chrom, snp.pos)]}:{snp.b6_gt}>{snp.cast_gt}"
        for snp in snps
    )


def _rank_key(row: dict) -> tuple:
    return (
        row.get("Filter_Pass") is not True,
        _int(row.get("SNP_Distance_To_RNA_3p")),
        _int(row.get("Primer_To_SNP_Distance_Toward_PolyDT")),
        _float(row.get("Primer3_Penalty")),
        abs(_float(row.get("Primer_Tm")) - 60.0),
    )


def _append_reason(existing: str, reason: str) -> str:
    return reason if not existing else f"{existing};{reason}"


def _has_3p_gc_clamp(seq: str) -> bool:
    return bool(seq) and seq[-1] in {"G", "C"}


def _has_homopolymer(seq: str, length: int) -> bool:
    for base in "ACGT":
        if base * length in seq:
            return True
    return False


def _in_range(value, lower, upper) -> bool:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return False
    if math.isnan(parsed):
        return False
    return float(lower) <= parsed <= float(upper)


def _float(value) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return float("inf")
    if math.isnan(parsed):
        return float("inf")
    return parsed


def _int(value) -> int:
    try:
        if value == "":
            return 10**12
        return int(float(value))
    except (TypeError, ValueError):
        return 10**12
