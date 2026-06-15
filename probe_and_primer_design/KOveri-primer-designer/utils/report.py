"""HTML report generation for KOveri outputs."""

from __future__ import annotations

from html import escape
from pathlib import Path

from .models import GuideHit, TargetRegion


def write_html_report(
    path: str | Path,
    target: TargetRegion,
    guide_hits: list[GuideHit],
    top_rows: list[dict],
    candidate_count: int,
    config: dict,
    filter_info: dict,
) -> None:
    """Write a consolidated portable HTML report with an inline coordinate figure."""
    path = Path(path)
    figure = _build_svg_figure(target, guide_hits, top_rows[:5])
    best = top_rows[0] if top_rows else None
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>KOveri Primer Design Report</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 28px; color: #202124; }}
    h1, h2 {{ margin: 0 0 12px; }}
    h1 {{ font-size: 26px; }}
    h2 {{ font-size: 18px; margin-top: 28px; }}
    .meta {{ color: #5f6368; margin-bottom: 20px; }}
    .summary {{ display: grid; grid-template-columns: repeat(4, minmax(120px, 1fr)); gap: 12px; margin: 18px 0; }}
    .metric {{ border: 1px solid #dadce0; border-radius: 6px; padding: 10px 12px; }}
    .metric b {{ display: block; font-size: 20px; margin-top: 3px; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 10px; font-size: 13px; }}
    th, td {{ border: 1px solid #dadce0; padding: 7px 8px; text-align: left; vertical-align: top; }}
    th {{ background: #f8f9fa; }}
    code {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 12px; }}
    .figure-wrap {{ border: 1px solid #dadce0; border-radius: 6px; padding: 12px; overflow-x: auto; }}
    .note {{ color: #5f6368; font-size: 13px; }}
  </style>
</head>
<body>
  <h1>KOveri Primer Design Report</h1>
  <div class="meta">Target interval: <code>{escape(target.chrom)}:{target.start}-{target.end}</code></div>
  <div class="summary">
    <div class="metric">Candidate pairs<b>{candidate_count}</b></div>
    <div class="metric">Top pairs reported<b>{len(top_rows)}</b></div>
    <div class="metric">Minimum SNP goal<b>{config["min_informative_snps"]}</b></div>
    <div class="metric">Selected filter tier<b>{escape(str(filter_info["selected_tier"]))}</b></div>
  </div>
  <h2>Post-Primer3 Filters</h2>
  {_filter_summary_table(filter_info)}
  {_best_pair_block(best)}
  <h2>Coordinate Figure</h2>
  <div class="figure-wrap">{figure}</div>
  <p class="note">Coordinates are 1-indexed inclusive. The BED output uses standard 0-based starts.</p>
  <h2>Selected Guide Hits</h2>
  {_guide_table(guide_hits, target)}
  <h2>Top Primer Candidates</h2>
  {_candidate_table(top_rows)}
</body>
</html>
"""
    path.write_text(html)


def _best_pair_block(best: dict | None) -> str:
    if not best:
        return "<p><b>No orderable primer pair passed the hard selectable filters.</b> Inspect <code>KOveri_primer_candidates.csv</code> for hard-filter failure reasons.</p>"
    return f"""
  <h2>Best Ranked Pair</h2>
  <table>
    <tr><th>Left primer</th><td><code>{escape(str(best["Left_Primer_Seq"]))}</code> ({escape(str(best["Chromosome"]))}:{best["Left_Primer_Start"]}-{best["Left_Primer_End"]})</td></tr>
    <tr><th>Right primer</th><td><code>{escape(str(best["Right_Primer_Seq"]))}</code> ({escape(str(best["Chromosome"]))}:{best["Right_Primer_Start"]}-{best["Right_Primer_End"]})</td></tr>
    <tr><th>Amplicon</th><td>{escape(str(best["Chromosome"]))}:{best["Amplicon_Start"]}-{best["Amplicon_End"]} ({best["Amplicon_Size"]} bp)</td></tr>
    <tr><th>Informative SNPs</th><td>{best["SNP_Count_In_Amplicon"]}: {escape(str(best["SNP_Positions_In_Amplicon"]))}</td></tr>
    <tr><th>Primer overlaps SNP</th><td>{escape(str(best["Primer_Overlaps_SNP"]))}</td></tr>
    <tr><th>Specificity pass</th><td>{escape(str(best.get("Primer_Pair_Specificity_Pass", "")))}</td></tr>
    <tr><th>Filter tier</th><td>{escape(str(best.get("Filter_Tier", "")))}</td></tr>
    <tr><th>Relaxed criteria used</th><td>{escape(str(best.get("Relaxed_Criteria_Used", "")))}</td></tr>
    <tr><th>Worst thermodynamic dG</th><td>{escape(str(best.get("Worst_Thermo_dG_kcal", "")))} kcal/mol</td></tr>
    <tr><th>Hard-filter failures</th><td>{escape(str(best.get("Hard_Filter_Failure_Reasons", "")))}</td></tr>
    <tr><th>Soft relaxed reasons</th><td>{escape(str(best.get("Soft_Filter_Relaxed_Reasons", "")))}</td></tr>
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
        "<p class='note'>Top outputs use only the strictest non-failed tier with at least one hard-filter-passing candidate. Relaxation order: amplicon size, specificity, then GC/Tm. Primers that fail hard filters, including missing 3-prime G/C clamp, are not orderable.</p>"
    )


def _build_svg_figure(
    target: TargetRegion,
    guide_hits: list[GuideHit],
    rows: list[dict],
) -> str:
    if not rows:
        return "<svg width='900' height='120' role='img'><text x='20' y='60'>No primer candidates generated.</text></svg>"

    region_start = min(
        [target.start]
        + [hit.start for hit in guide_hits]
        + [int(row["Amplicon_Start"]) for row in rows]
    )
    region_end = max(
        [target.end]
        + [hit.end for hit in guide_hits]
        + [int(row["Amplicon_End"]) for row in rows]
    )
    snps = _parse_snps_from_rows(rows)
    if snps:
        region_start = min(region_start, min(pos for pos, _ in snps))
        region_end = max(region_end, max(pos for pos, _ in snps))

    width = 1100
    margin_left = 120
    margin_right = 40
    axis_width = width - margin_left - margin_right
    row_height = 56
    height = 175 + len(rows) * row_height

    def xcoord(pos: int) -> float:
        if region_end == region_start:
            return margin_left
        return margin_left + ((pos - region_start) / (region_end - region_start)) * axis_width

    parts = [
        f"<svg width='{width}' height='{height}' viewBox='0 0 {width} {height}' role='img' aria-label='KOveri primer coordinate figure'>",
        "<rect x='0' y='0' width='100%' height='100%' fill='white'/>",
        "<style>.axis{stroke:#3c4043;stroke-width:1}.amp{stroke:#5f6368;stroke-width:2}.primerL{fill:#1a73e8}.primerR{fill:#d93025}.guide{fill:#188038}.snp{stroke:#f9ab00;stroke-width:2}.label{font:12px sans-serif;fill:#202124}.small{font:10px sans-serif;fill:#5f6368}</style>",
    ]

    axis_y = 42
    parts.append(f"<line class='axis' x1='{margin_left}' y1='{axis_y}' x2='{width - margin_right}' y2='{axis_y}'/>")
    for tick in _ticks(region_start, region_end, 5):
        x = xcoord(tick)
        parts.append(f"<line class='axis' x1='{x:.1f}' y1='{axis_y - 5}' x2='{x:.1f}' y2='{axis_y + 5}'/>")
        parts.append(f"<text class='small' x='{x:.1f}' y='{axis_y - 10}' text-anchor='middle'>{tick}</text>")
    parts.append(f"<text class='label' x='12' y='{axis_y + 4}'>Genomic coords</text>")

    guide_y = 88
    parts.append(f"<line class='amp' x1='{xcoord(target.start):.1f}' y1='{guide_y}' x2='{xcoord(target.end):.1f}' y2='{guide_y}'/>")
    parts.append(f"<text class='label' x='12' y='{guide_y + 4}'>gRNAs</text>")
    for hit in guide_hits:
        x1 = xcoord(hit.start)
        x2 = xcoord(hit.end)
        x = min(x1, x2)
        w = max(4, abs(x2 - x1))
        parts.append(f"<rect class='guide' x='{x:.1f}' y='{guide_y - 8}' width='{w:.1f}' height='16' rx='2'/>")
        strand_mark = ">" if hit.strand == "+" else "<"
        parts.append(
            f"<text class='small' x='{(x1 + x2) / 2:.1f}' y='{guide_y + 25}' text-anchor='middle'>{escape(hit.guide_id)} {strand_mark}</text>"
        )

    snp_y1 = 118
    snp_y2 = height - 24
    parts.append(f"<text class='label' x='12' y='{snp_y1 + 4}'>SNPs</text>")
    for pos, label in snps:
        x = xcoord(pos)
        parts.append(f"<line class='snp' x1='{x:.1f}' y1='{snp_y1}' x2='{x:.1f}' y2='{snp_y2}'/>")
        parts.append(f"<text class='small' x='{x:.1f}' y='{snp_y1 - 6}' text-anchor='middle'>{pos}</text>")

    start_y = 150
    for idx, row in enumerate(rows, start=1):
        y = start_y + (idx - 1) * row_height
        amp_start = int(row["Amplicon_Start"])
        amp_end = int(row["Amplicon_End"])
        left_start = int(row["Left_Primer_Start"])
        left_end = int(row["Left_Primer_End"])
        right_start = int(row["Right_Primer_Start"])
        right_end = int(row["Right_Primer_End"])
        parts.append(f"<text class='label' x='12' y='{y + 4}'>Pair {idx}</text>")
        parts.append(f"<line class='amp' x1='{xcoord(amp_start):.1f}' y1='{y}' x2='{xcoord(amp_end):.1f}' y2='{y}'/>")
        parts.append(
            f"<rect class='primerL' x='{xcoord(left_start):.1f}' y='{y - 7}' width='{max(5, xcoord(left_end) - xcoord(left_start)):.1f}' height='14' rx='2'/>"
        )
        parts.append(
            f"<rect class='primerR' x='{xcoord(right_start):.1f}' y='{y - 7}' width='{max(5, xcoord(right_end) - xcoord(right_start)):.1f}' height='14' rx='2'/>"
        )
        parts.append(
            f"<text class='small' x='{xcoord(amp_end):.1f}' y='{y + 20}' text-anchor='end'>{row['Amplicon_Size']} bp, {row['SNP_Count_In_Amplicon']} SNPs</text>"
        )

    parts.append("<text class='small' x='120' y='" + str(height - 6) + "'>Blue: left primer. Red: right primer. Green: gRNAs. Orange: B6/Cast informative SNPs.</text>")
    parts.append("</svg>")
    return "\n".join(parts)


def _parse_snps_from_rows(rows: list[dict]) -> list[tuple[int, str]]:
    snps: dict[int, str] = {}
    for row in rows:
        for item in str(row.get("SNP_Positions_In_Amplicon", "")).split(";"):
            if not item:
                continue
            fields = item.split(":")
            if len(fields) < 3:
                continue
            try:
                snps[int(fields[1])] = item
            except ValueError:
                continue
    return sorted(snps.items())


def _ticks(start: int, end: int, count: int) -> list[int]:
    if count <= 1 or end <= start:
        return [start, end]
    return [round(start + (end - start) * i / (count - 1)) for i in range(count)]


def _guide_table(guide_hits: list[GuideHit], target: TargetRegion) -> str:
    rows = []
    for hit in guide_hits:
        selected = target.chrom == hit.chrom and target.start <= hit.start and hit.end <= target.end
        rows.append(
            "<tr>"
            f"<td>{escape(hit.guide_id)}</td>"
            f"<td><code>{escape(hit.sequence)}</code></td>"
            f"<td>{escape(hit.chrom)}:{hit.start}-{hit.end}</td>"
            f"<td>{escape(hit.strand)}</td>"
            f"<td>{selected}</td>"
            "</tr>"
        )
    return (
        "<table><tr><th>Guide</th><th>Sequence</th><th>Coordinate</th><th>Strand</th><th>Selected</th></tr>"
        + "\n".join(rows)
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
            f"<td>{escape(str(row['Chromosome']))}:{row['Amplicon_Start']}-{row['Amplicon_End']}</td>"
            f"<td>{row['Amplicon_Size']}</td>"
            f"<td>{row['SNP_Count_In_Amplicon']}</td>"
            f"<td><code>{escape(str(row['Left_Primer_Seq']))}</code></td>"
            f"<td><code>{escape(str(row['Right_Primer_Seq']))}</code></td>"
            f"<td>{escape(str(row['Primer_Overlaps_SNP']))}</td>"
            f"<td>{escape(str(row.get('Primer_Pair_Specificity_Pass', '')))}</td>"
            f"<td>{escape(str(row.get('Worst_Thermo_dG_kcal', '')))}</td>"
            f"<td>{escape(str(row.get('Relaxed_Criteria_Used', '')))}</td>"
            "</tr>"
        )
    return (
        "<table><tr><th>Rank</th><th>Tier</th><th>Amplicon</th><th>bp</th><th>SNPs</th><th>Left primer</th><th>Right primer</th><th>Primer on SNP</th><th>Specificity pass</th><th>Worst dG</th><th>Relaxed</th></tr>"
        + "\n".join(body)
        + "</table>"
    )
