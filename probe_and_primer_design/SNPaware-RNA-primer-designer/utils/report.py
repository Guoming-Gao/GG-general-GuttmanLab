"""HTML report generation for SNP-aware RNA primer outputs."""

from __future__ import annotations

import re
from html import escape
from pathlib import Path

from .models import CommonJunction, TranscriptModel, TranscriptTemplate


def write_html_report(
    path: str | Path,
    gene_name: str,
    selected_transcript: TranscriptModel,
    template: TranscriptTemplate,
    common_junctions: list[CommonJunction],
    top_rows: list[dict],
    candidate_count: int,
    config: dict,
    filter_info: dict,
    fallback_used: bool,
) -> None:
    path = Path(path)
    best = top_rows[0] if top_rows else None
    figure = _build_cdna_figure(template, common_junctions, top_rows[:5])
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>SNP-aware RNA Primer Design Report</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 28px; color: #202124; }}
    h1, h2 {{ margin: 0 0 12px; }}
    h1 {{ font-size: 26px; }}
    h2 {{ font-size: 18px; margin-top: 28px; }}
    .meta {{ color: #5f6368; margin-bottom: 20px; }}
    .summary {{ display: grid; grid-template-columns: repeat(4, minmax(120px, 1fr)); gap: 12px; margin: 18px 0; }}
    .metric {{ border: 1px solid #dadce0; border-radius: 6px; padding: 10px 12px; }}
    .metric b {{ display: block; font-size: 20px; margin-top: 3px; }}
    .figure-wrap {{ border: 1px solid #dadce0; border-radius: 6px; padding: 12px; overflow-x: auto; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 10px; font-size: 13px; }}
    th, td {{ border: 1px solid #dadce0; padding: 7px 8px; text-align: left; vertical-align: top; }}
    th {{ background: #f8f9fa; }}
    code {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 12px; }}
    .note {{ color: #5f6368; font-size: 13px; }}
  </style>
</head>
<body>
  <h1>SNP-aware RNA Primer Design Report</h1>
  <div class="meta">Gene <code>{escape(gene_name)}</code>, transcript <code>{escape(selected_transcript.transcript_id)}</code>, {escape(template.chrom)}:{selected_transcript.start}-{selected_transcript.end} ({escape(template.strand)})</div>
  <div class="summary">
    <div class="metric">Candidate pairs<b>{candidate_count}</b></div>
    <div class="metric">Top pairs reported<b>{len(top_rows)}</b></div>
    <div class="metric">Common junctions<b>{len(common_junctions)}</b></div>
    <div class="metric">Selected filter tier<b>{escape(str(filter_info["selected_tier"]))}</b></div>
  </div>
  <p class="note">Crossing-primer search was attempted first. Fallback junction-spanning search used: <b>{fallback_used}</b>.</p>
  <p class="note">Product length is calculated on the spliced cDNA template. Genomic span is reported separately and is not the PCR product length. Junction-crossing primers are split across exon arms and may not produce one contiguous genomic BLAT/BLAST hit.</p>
  <h2>Post-Primer3 Filters</h2>
  {_filter_summary_table(filter_info)}
  {_best_pair_block(best)}
  <h2>cDNA Coordinate Diagram</h2>
  <div class="figure-wrap">{figure}</div>
  <h2>Common Junctions</h2>
  {_junction_table(common_junctions)}
  <h2>Top Primer Candidates</h2>
  {_candidate_table(top_rows)}
</body>
</html>
"""
    path.write_text(html)


def _best_pair_block(best: dict | None) -> str:
    if not best:
        return "<p><b>No orderable primer pair passed the hard selectable filters.</b> Inspect <code>SNPaware_RNA_primer_candidates.csv</code> for failure reasons.</p>"
    return f"""
  <h2>Best Ranked Pair</h2>
  <table>
    <tr><th>Design mode</th><td>{escape(str(best["Design_Mode"]))}</td></tr>
    <tr><th>Left primer</th><td><code>{escape(str(best["Left_Primer_Seq"]))}</code><br>{escape(str(best["Left_Primer_Genomic_Segments"]))}</td></tr>
    <tr><th>Right primer</th><td><code>{escape(str(best["Right_Primer_Seq"]))}</code><br>{escape(str(best["Right_Primer_Genomic_Segments"]))}</td></tr>
    <tr><th>cDNA amplicon</th><td>{best["Amplicon_cDNA_Start"]}-{best["Amplicon_cDNA_End"]} ({best.get("Amplicon_cDNA_Size", best["Amplicon_Size"])} bp cDNA product)</td></tr>
    <tr><th>Genomic span</th><td>{escape(str(best.get("Amplicon_Genomic_Span", "")))} bp genomic span, not PCR product length</td></tr>
    <tr><th>Genomic segments</th><td>{escape(str(best["Amplicon_Genomic_Segments"]))}</td></tr>
    <tr><th>Target junction</th><td>{escape(str(best.get("Target_Junction_Genomic", "")))}</td></tr>
    <tr><th>Primary crossing junction</th><td>{escape(str(best.get("Primary_Crossing_Junction", "")))}</td></tr>
    <tr><th>Crossing primer side</th><td>{escape(str(best.get("Crossing_Primer_Side", "")))}</td></tr>
    <tr><th>Left primer junction arms</th><td>5-prime/3-prime: {escape(str(best.get("Left_Primer_Junction_5p_Arm_Bases", "")))}/{escape(str(best.get("Left_Primer_Junction_3p_Arm_Bases", "")))} bases</td></tr>
    <tr><th>Right primer junction arms</th><td>5-prime/3-prime: {escape(str(best.get("Right_Primer_Junction_5p_Arm_Bases", "")))}/{escape(str(best.get("Right_Primer_Junction_3p_Arm_Bases", "")))} bases</td></tr>
    <tr><th>Crossing junctions</th><td>{escape(str(best.get("Crossing_Junctions", "")))}</td></tr>
    <tr><th>Spanned junctions</th><td>{escape(str(best.get("Spanned_Junctions", "")))}</td></tr>
    <tr><th>Informative SNPs</th><td>{best["SNP_Count_In_Amplicon"]}: {escape(str(best["SNP_Positions_In_Amplicon"]))}</td></tr>
    <tr><th>Primer overlaps SNP</th><td>{escape(str(best["Primer_Overlaps_SNP"]))}</td></tr>
    <tr><th>Specificity pass</th><td>{escape(str(best.get("Primer_Pair_Specificity_Pass", "")))}</td></tr>
    <tr><th>Filter tier</th><td>{escape(str(best.get("Filter_Tier", "")))}</td></tr>
    <tr><th>Relaxed criteria used</th><td>{escape(str(best.get("Relaxed_Criteria_Used", "")))}</td></tr>
    <tr><th>Worst thermodynamic dG</th><td>{escape(str(best.get("Worst_Thermo_dG_kcal", "")))} kcal/mol</td></tr>
    <tr><th>Hard-filter failures</th><td>{escape(str(best.get("Hard_Filter_Failure_Reasons", "")))}</td></tr>
  </table>
"""


def _filter_summary_table(filter_info: dict) -> str:
    counts = filter_info.get("tier_counts", {})
    explanations = filter_info.get("tier_explanations", {})
    rows = []
    for tier in ["strict_all", "relaxed_amplicon", "relaxed_specificity", "relaxed_gc_tm", "failed"]:
        rows.append(
            f"<tr><td>{escape(tier)}</td><td>{counts.get(tier, 0)}</td><td>{escape(str(explanations.get(tier, '')))}</td></tr>"
        )
    return (
        f"<p class='note'>Hard-filter passing candidates: {filter_info.get('hard_pass_count', 0)}. "
        f"Hard-filter failing candidates: {filter_info.get('hard_fail_count', 0)}.</p>"
        "<table><tr><th>Tier</th><th>Candidate count</th><th>Explanation</th></tr>"
        + "\n".join(rows)
        + "</table>"
    )


def _junction_table(junctions: list[CommonJunction]) -> str:
    body = []
    for idx, junction in enumerate(junctions, start=1):
        body.append(
            "<tr>"
            f"<td>{idx}</td>"
            f"<td>{escape(junction.chrom)}:{junction.upstream_exon_end_genomic}|{junction.downstream_exon_start_genomic}</td>"
            f"<td>{escape(junction.strand)}</td>"
            f"<td>{junction.cdna_left}|{junction.cdna_right}</td>"
            f"<td>{junction.upstream_exon_number}|{junction.downstream_exon_number}</td>"
            f"<td>{junction.transcript_count}</td>"
            "</tr>"
        )
    return (
        "<table><tr><th>#</th><th>Genomic junction</th><th>Strand</th><th>cDNA junction</th><th>Exons</th><th>Curated transcript count</th></tr>"
        + "\n".join(body)
        + "</table>"
    )


def _candidate_table(rows: list[dict]) -> str:
    if not rows:
        return "<p>No top candidates available.</p>"
    body = []
    for idx, row in enumerate(rows, start=1):
        body.append(
            "<tr>"
            f"<td>{idx}</td>"
            f"<td>{escape(str(row.get('Filter_Tier', '')))}</td>"
            f"<td>{escape(str(row.get('Design_Mode', '')))}</td>"
            f"<td>{row['Amplicon_cDNA_Start']}-{row['Amplicon_cDNA_End']}</td>"
            f"<td>{row.get('Amplicon_cDNA_Size', row['Amplicon_Size'])}</td>"
            f"<td>{escape(str(row.get('Amplicon_Genomic_Span', '')))}</td>"
            f"<td>{row['SNP_Count_In_Amplicon']}</td>"
            f"<td><code>{escape(str(row['Left_Primer_Seq']))}</code></td>"
            f"<td><code>{escape(str(row['Right_Primer_Seq']))}</code></td>"
            f"<td>{escape(str(row.get('Primary_Crossing_Junction', '')))}</td>"
            f"<td>{escape(str(row.get('Primer_Pair_Specificity_Pass', '')))}</td>"
            f"<td>{escape(str(row.get('Relaxed_Criteria_Used', '')))}</td>"
            "</tr>"
        )
    return (
        "<table><tr><th>Rank</th><th>Tier</th><th>Mode</th><th>cDNA amplicon</th><th>cDNA bp</th><th>Genomic span</th><th>SNPs</th><th>Left primer</th><th>Right primer</th><th>Primary crossing junction</th><th>Specificity pass</th><th>Relaxed</th></tr>"
        + "\n".join(body)
        + "</table>"
    )


def _build_cdna_figure(
    template: TranscriptTemplate,
    junctions: list[CommonJunction],
    rows: list[dict],
) -> str:
    if not rows:
        return "<svg width='900' height='120' role='img'><text x='20' y='60'>No primer candidates generated.</text></svg>"

    snps = _parse_snps_from_rows(rows)
    width = 1180
    margin_left = 110
    margin_right = 45
    axis_width = width - margin_left - margin_right
    row_height = 60
    height = 220 + len(rows) * row_height

    def xcoord(pos: int) -> float:
        if template.length <= 1:
            return margin_left
        return margin_left + ((pos - 1) / (template.length - 1)) * axis_width

    parts = [
        f"<svg width='{width}' height='{height}' viewBox='0 0 {width} {height}' role='img' aria-label='SNP-aware RNA cDNA coordinate diagram'>",
        "<rect x='0' y='0' width='100%' height='100%' fill='white'/>",
        "<style>.axis{stroke:#3c4043;stroke-width:1}.exon{fill:#e8f0fe;stroke:#1a73e8;stroke-width:1}.amp{stroke:#5f6368;stroke-width:2}.primerL{fill:#1a73e8}.primerR{fill:#d93025}.junction{stroke:#188038;stroke-width:2}.primary{stroke:#7b1fa2;stroke-width:3}.snp{stroke:#f9ab00;stroke-width:2}.label{font:12px sans-serif;fill:#202124}.small{font:10px sans-serif;fill:#5f6368}.tiny{font:9px sans-serif;fill:#5f6368}</style>",
    ]

    axis_y = 44
    parts.append(f"<line class='axis' x1='{margin_left}' y1='{axis_y}' x2='{width - margin_right}' y2='{axis_y}'/>")
    for tick in _ticks(1, template.length, 6):
        x = xcoord(tick)
        parts.append(f"<line class='axis' x1='{x:.1f}' y1='{axis_y - 5}' x2='{x:.1f}' y2='{axis_y + 5}'/>")
        parts.append(f"<text class='small' x='{x:.1f}' y='{axis_y - 10}' text-anchor='middle'>cDNA {tick}</text>")
    parts.append(f"<text class='label' x='12' y='{axis_y + 4}'>cDNA axis</text>")

    exon_y = 86
    parts.append(f"<text class='label' x='12' y='{exon_y + 5}'>Exons</text>")
    for idx, segment in enumerate(template.segments, start=1):
        x1 = xcoord(segment.cdna_start)
        x2 = xcoord(segment.cdna_end)
        parts.append(
            f"<rect class='exon' x='{x1:.1f}' y='{exon_y - 10}' width='{max(3, x2 - x1):.1f}' height='20'>"
            f"<title>Exon {idx}: cDNA {segment.cdna_start}-{segment.cdna_end}; {segment.exon.chrom}:{segment.exon.start}-{segment.exon.end}</title></rect>"
        )

    junction_y1 = 116
    junction_y2 = height - 24
    parts.append(f"<text class='label' x='12' y='{junction_y1 + 4}'>Junctions</text>")
    primary_junctions = {row.get("Primary_Crossing_Junction", "") for row in rows if row.get("Primary_Crossing_Junction")}
    for junction in junctions:
        label = _format_junction(junction)
        cls = "primary" if label in primary_junctions else "junction"
        x = xcoord(junction.cdna_right)
        parts.append(f"<line class='{cls}' x1='{x:.1f}' y1='{junction_y1}' x2='{x:.1f}' y2='{junction_y2}'><title>{escape(label)}</title></line>")
        parts.append(f"<text class='tiny' x='{x:.1f}' y='{junction_y1 - 6}' text-anchor='middle'>{junction.cdna_left}|{junction.cdna_right}</text>")

    snp_y1 = 142
    parts.append(f"<text class='label' x='12' y='{snp_y1 + 4}'>SNPs</text>")
    for pos, label in snps:
        x = xcoord(pos)
        parts.append(f"<line class='snp' x1='{x:.1f}' y1='{snp_y1}' x2='{x:.1f}' y2='{junction_y2}'><title>{escape(label)}</title></line>")

    start_y = 180
    for idx, row in enumerate(rows, start=1):
        y = start_y + (idx - 1) * row_height
        amp_start = int(row["Amplicon_cDNA_Start"])
        amp_end = int(row["Amplicon_cDNA_End"])
        left_start = int(row["Left_Primer_cDNA_Start"])
        left_end = int(row["Left_Primer_cDNA_End"])
        right_start = int(row["Right_Primer_cDNA_Start"])
        right_end = int(row["Right_Primer_cDNA_End"])
        parts.append(f"<text class='label' x='12' y='{y + 4}'>Pair {idx}</text>")
        parts.append(f"<line class='amp' x1='{xcoord(amp_start):.1f}' y1='{y}' x2='{xcoord(amp_end):.1f}' y2='{y}'><title>cDNA product {row.get('Amplicon_cDNA_Size', row['Amplicon_Size'])} bp; genomic span {escape(str(row.get('Amplicon_Genomic_Span', '')))} bp</title></line>")
        parts.append(f"<rect class='primerL' x='{xcoord(left_start):.1f}' y='{y - 7}' width='{max(5, xcoord(left_end) - xcoord(left_start)):.1f}' height='14'><title>Left primer cDNA {left_start}-{left_end}; {escape(str(row.get('Left_Primer_Genomic_Segments', '')))}</title></rect>")
        parts.append(f"<rect class='primerR' x='{xcoord(right_start):.1f}' y='{y - 7}' width='{max(5, xcoord(right_end) - xcoord(right_start)):.1f}' height='14'><title>Right primer cDNA {right_start}-{right_end}; {escape(str(row.get('Right_Primer_Genomic_Segments', '')))}</title></rect>")
        primary = row.get("Primary_Crossing_Junction", "")
        if primary:
            match = re.search(r"cDNA(\d+)\|(\d+)", str(primary))
            if match:
                x = xcoord(int(match.group(2)))
                parts.append(f"<circle cx='{x:.1f}' cy='{y}' r='5' fill='#7b1fa2'><title>Primary crossing junction: {escape(str(primary))}</title></circle>")
        parts.append(
            f"<text class='small' x='{xcoord(amp_end):.1f}' y='{y + 22}' text-anchor='end'>"
            f"{row.get('Amplicon_cDNA_Size', row['Amplicon_Size'])} bp cDNA, genomic span {escape(str(row.get('Amplicon_Genomic_Span', '')))} bp, {row['SNP_Count_In_Amplicon']} SNPs</text>"
        )

    parts.append(f"<text class='small' x='{margin_left}' y='{height - 7}'>Blue: left primer. Red: right primer. Green: common junctions. Purple: primary crossed junction. Orange: informative SNPs.</text>")
    parts.append("</svg>")
    return "\n".join(parts)


def _parse_snps_from_rows(rows: list[dict]) -> list[tuple[int, str]]:
    snps: dict[int, str] = {}
    for row in rows:
        for item in str(row.get("SNP_Positions_In_Amplicon", "")).split(";"):
            if not item:
                continue
            match = re.search(r":cDNA(\d+):", item)
            if match:
                snps[int(match.group(1))] = item
    return sorted(snps.items())


def _format_junction(junction: CommonJunction) -> str:
    return (
        f"{junction.chrom}:{junction.upstream_exon_end_genomic}|"
        f"{junction.downstream_exon_start_genomic}:{junction.strand}:"
        f"cDNA{junction.cdna_left}|{junction.cdna_right}"
    )


def _ticks(start: int, end: int, count: int) -> list[int]:
    if count <= 1 or end <= start:
        return [start, end]
    return [round(start + (end - start) * i / (count - 1)) for i in range(count)]
