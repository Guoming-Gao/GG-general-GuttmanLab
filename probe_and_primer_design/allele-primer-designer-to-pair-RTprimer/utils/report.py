"""HTML report generation for allele primer outputs."""

from __future__ import annotations

from html import escape
from pathlib import Path

from .models import TranscriptModel, TranscriptTemplate


def write_html_report(
    path: str | Path,
    gene_name: str,
    selected_transcript: TranscriptModel,
    template: TranscriptTemplate,
    snp_rows: list[dict],
    candidate_rows: list[dict],
    top_insert_rows: list[dict],
    top_oligo_rows: list[dict],
    config: dict,
) -> None:
    path = Path(path)
    best = top_insert_rows[0] if top_insert_rows else None
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Allele Primer Design Report</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 28px; color: #202124; }}
    h1 {{ font-size: 26px; margin: 0 0 8px; }}
    h2 {{ font-size: 18px; margin: 26px 0 10px; }}
    .meta, .note {{ color: #5f6368; }}
    .summary {{ display: grid; grid-template-columns: repeat(5, minmax(120px, 1fr)); gap: 12px; margin: 18px 0; }}
    .metric {{ border: 1px solid #dadce0; border-radius: 6px; padding: 10px 12px; }}
    .metric b {{ display: block; font-size: 20px; margin-top: 3px; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 10px; font-size: 13px; }}
    th, td {{ border: 1px solid #dadce0; padding: 7px 8px; text-align: left; vertical-align: top; }}
    th {{ background: #f8f9fa; }}
    code {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 12px; }}
  </style>
</head>
<body>
  <h1>Allele Primer Design Report</h1>
  <div class="meta">Gene <code>{escape(gene_name)}</code>, transcript <code>{escape(selected_transcript.transcript_id)}</code>, {escape(template.chrom)}:{selected_transcript.start}-{selected_transcript.end} ({escape(template.strand)}), handle <code>{escape(str(config.get('handle', '2PUNI')))}</code></div>
  <div class="summary">
    <div class="metric">Informative SNPs<b>{len(snp_rows)}</b></div>
    <div class="metric">Candidate inserts<b>{len(candidate_rows)}</b></div>
    <div class="metric">Passing inserts<b>{sum(row.get("Filter_Pass") is True for row in candidate_rows)}</b></div>
    <div class="metric">Top inserts<b>{len(top_insert_rows)}</b></div>
    <div class="metric">Top oligos<b>{len(top_oligo_rows)}</b></div>
  </div>
  <p class="note">All coordinates labeled cDNA are spliced mature RNA coordinates in RNA 5-prime to 3-prime order. The SNP window is checked only from the primer insert 3-prime end toward the transcript/poly-dT 3-prime end.</p>
  {_best_block(best)}
  <h2>Top Inserts</h2>
  {_top_insert_table(top_insert_rows)}
  <h2>Top Order Oligos</h2>
  {_top_oligo_table(top_oligo_rows)}
  <h2>Closest Informative SNPs</h2>
  {_snp_table(snp_rows[:20])}
</body>
</html>
"""
    path.write_text(html)


def write_consolidated_html_report(
    path: str | Path,
    batch_rows: list[dict],
    top_oligo_rows: list[dict],
    config: dict,
) -> None:
    """Write a batch-level report containing statuses and all order-ready oligos."""
    successes = sum(row.get("Status") == "success" for row in batch_rows)
    failures = len(batch_rows) - successes
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Consolidated Allele Primer Design Report</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 28px; color: #202124; }}
    h1 {{ font-size: 26px; margin: 0 0 8px; }} h2 {{ font-size: 18px; margin: 26px 0 10px; }}
    .meta {{ color: #5f6368; }} .summary {{ display: flex; gap: 12px; margin: 18px 0; }}
    .metric {{ border: 1px solid #dadce0; border-radius: 6px; padding: 10px 14px; }}
    .metric b {{ display: block; font-size: 20px; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 10px; font-size: 13px; }}
    th, td {{ border: 1px solid #dadce0; padding: 7px 8px; text-align: left; vertical-align: top; }}
    th {{ background: #f8f9fa; }} code {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 12px; }}
  </style>
</head>
<body>
  <h1>Consolidated Allele Primer Design Report</h1>
  <div class="meta">Handle <code>{escape(str(config.get('handle', '2PUNI')))}</code>; SNP window toward poly-dT {config['snp_window_toward_polydt']} nt; BLAST specificity {escape(str(config.get('blast_specificity', True)))}</div>
  <div class="summary"><div class="metric">Successful genes<b>{successes}</b></div><div class="metric">Failed genes<b>{failures}</b></div><div class="metric">Order-ready oligos<b>{len(top_oligo_rows)}</b></div></div>
  <h2>Gene Summary</h2>
  {_batch_table(batch_rows)}
  <h2>Top Order-Ready Oligos</h2>
  {_consolidated_oligo_table(top_oligo_rows)}
</body>
</html>
"""
    Path(path).write_text(html)


def _batch_table(rows: list[dict]) -> str:
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{escape(str(row.get(key, '')))}</td>" for key in ("Gene", "Status", "Selected_Transcript", "Strand", "Informative_SNPs", "Candidate_Inserts", "Top_Inserts", "Top_Oligos", "Error")) + "</tr>")
    return "<table><tr><th>Gene</th><th>Status</th><th>Transcript</th><th>Strand</th><th>SNPs</th><th>Candidates</th><th>Top inserts</th><th>Top oligos</th><th>Error</th></tr>" + "\n".join(body) + "</table>"


def _consolidated_oligo_table(rows: list[dict]) -> str:
    if not rows:
        return "<p>No order-ready oligos were found.</p>"
    body = []
    for row in rows:
        body.append(
            "<tr>"
            f"<td>{escape(str(row['Gene']))}</td><td>{escape(str(row['Transcript_ID']))}</td><td>{escape(str(row['Strand']))}</td>"
            f"<td>{row['Insert_Rank']}</td><td><code>{escape(str(row['Primer_Insert_Seq']))}</code></td>"
            f"<td>{escape(str(row['Handle_Name']))}</td><td><code>{escape(str(row['Final_Order_Sequence']))}</code></td>"
            f"<td>{row['SNP_cDNA']}</td><td>{row['SNP_Genomic']}</td><td>{row['SNP_Distance_To_RNA_3p']}</td>"
            f"<td>{row['Primer_To_SNP_Distance_Toward_PolyDT']}</td><td>{row['Primer_Tm']}</td><td>{row['Primer_GC']}</td>"
            f"<td>{escape(str(row.get('Primer_Specificity_Pass', '')))}</td></tr>"
        )
    return "<table><tr><th>Gene</th><th>Transcript</th><th>Strand</th><th>Rank</th><th>Insert</th><th>Handle</th><th>Final order sequence</th><th>SNP cDNA</th><th>SNP genomic</th><th>SNP to RNA 3-prime</th><th>Primer to SNP</th><th>Tm</th><th>GC%</th><th>Specificity</th></tr>" + "\n".join(body) + "</table>"


def _best_block(best: dict | None) -> str:
    if not best:
        return "<p><b>No passing primer insert was found.</b> Inspect <code>allele_primer_candidate_inserts.csv</code> for failure reasons.</p>"
    return f"""
  <h2>Best Insert</h2>
  <table>
    <tr><th>Primer insert</th><td><code>{escape(str(best["Primer_Insert_Seq"]))}</code></td></tr>
    <tr><th>Primer cDNA interval</th><td>{best["Primer_cDNA_Start"]}-{best["Primer_cDNA_End"]}</td></tr>
    <tr><th>Primer genomic segment</th><td>{escape(str(best["Primer_Genomic_Segments"]))}</td></tr>
    <tr><th>Target SNP</th><td>cDNA {best["SNP_cDNA"]}, genomic {best["Chromosome"]}:{best["SNP_Genomic"]}, B6 {escape(str(best["SNP_B6_GT"]))}, Cast {escape(str(best["SNP_CAST_GT"]))}</td></tr>
    <tr><th>SNP distance to RNA 3-prime end</th><td>{best["SNP_Distance_To_RNA_3p"]} spliced nt</td></tr>
    <tr><th>Primer to SNP toward poly-dT</th><td>{best["Primer_To_SNP_Distance_Toward_PolyDT"]} nt</td></tr>
    <tr><th>Tm / GC</th><td>{best["Primer_Tm"]} C / {best["Primer_GC"]}%</td></tr>
  </table>
"""


def _top_insert_table(rows: list[dict]) -> str:
    if not rows:
        return "<p>No passing inserts.</p>"
    body = []
    for idx, row in enumerate(rows, start=1):
        body.append(
            "<tr>"
            f"<td>{idx}</td>"
            f"<td><code>{escape(str(row['Primer_Insert_Seq']))}</code></td>"
            f"<td>{row['SNP_cDNA']}</td>"
            f"<td>{row['SNP_Genomic']}</td>"
            f"<td>{row['SNP_Distance_To_RNA_3p']}</td>"
            f"<td>{row['Primer_To_SNP_Distance_Toward_PolyDT']}</td>"
            f"<td>{escape(str(row['Primer_Genomic_Segments']))}</td>"
            f"<td>{escape(str(row.get('Primer_Specificity_Pass', '')))}</td>"
            "</tr>"
        )
    return (
        "<table><tr><th>Rank</th><th>Insert</th><th>SNP cDNA</th><th>SNP genomic</th><th>SNP to RNA 3-prime</th><th>Primer to SNP</th><th>Primer genomic segment</th><th>Specificity</th></tr>"
        + "\n".join(body)
        + "</table>"
    )


def _top_oligo_table(rows: list[dict]) -> str:
    if not rows:
        return "<p>No order oligos.</p>"
    body = []
    for idx, row in enumerate(rows, start=1):
        body.append(
            "<tr>"
            f"<td>{idx}</td>"
            f"<td>{escape(str(row['Handle_Name']))}</td>"
            f"<td><code>{escape(str(row['Final_Order_Sequence']))}</code></td>"
            f"<td>{row['SNP_cDNA']}</td>"
            f"<td>{row['SNP_Distance_To_RNA_3p']}</td>"
            "</tr>"
        )
    return (
        "<table><tr><th>#</th><th>Handle</th><th>Final order sequence</th><th>SNP cDNA</th><th>SNP to RNA 3-prime</th></tr>"
        + "\n".join(body)
        + "</table>"
    )


def _snp_table(rows: list[dict]) -> str:
    if not rows:
        return "<p>No informative SNPs found in selected transcript.</p>"
    body = []
    for row in rows:
        body.append(
            "<tr>"
            f"<td>{row['SNP_cDNA']}</td>"
            f"<td>{row['SNP_Genomic']}</td>"
            f"<td>{row['SNP_Distance_To_RNA_3p']}</td>"
            f"<td>{escape(str(row['SNP_B6_GT']))}</td>"
            f"<td>{escape(str(row['SNP_CAST_GT']))}</td>"
            "</tr>"
        )
    return (
        "<table><tr><th>cDNA</th><th>Genomic</th><th>Distance to RNA 3-prime</th><th>B6</th><th>Cast</th></tr>"
        + "\n".join(body)
        + "</table>"
    )
