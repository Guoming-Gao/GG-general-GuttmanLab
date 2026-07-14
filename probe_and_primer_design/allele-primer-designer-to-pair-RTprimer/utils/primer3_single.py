"""Primer3 wrapper for one RNA-sense gene-specific primer insert."""

from __future__ import annotations

import subprocess


def build_excluded_regions(excluded_positions: list[int], template_start: int, template_end: int) -> str:
    intervals = []
    for pos in sorted(set(excluded_positions)):
        if template_start <= pos <= template_end:
            intervals.append(f"{pos - template_start},1")
    return " ".join(intervals)


def run_primer3_left_primer_design(
    sequence_id: str,
    template: str,
    template_start_cdna: int,
    excluded_positions_cdna: list[int],
    config: dict,
    primer3_core: str,
) -> list[dict]:
    """Run Primer3 on a short RNA-sense window and return left primers.

    Primer3's left primer sequence is reported in the same 5'->3' orientation
    as the supplied template. Because the template is mature RNA/cDNA in gene
    direction, this is the orderable gene-specific primer insert.
    """
    template_end_cdna = template_start_cdna + len(template) - 1
    excluded = build_excluded_regions(excluded_positions_cdna, template_start_cdna, template_end_cdna)

    lines = [
        f"SEQUENCE_ID={sequence_id}",
        f"SEQUENCE_TEMPLATE={template}",
        "PRIMER_TASK=pick_left_only",
        "PRIMER_PICK_LEFT_PRIMER=1",
        "PRIMER_PICK_RIGHT_PRIMER=0",
        "PRIMER_PICK_INTERNAL_OLIGO=0",
        f"PRIMER_NUM_RETURN={config['primer3_num_return']}",
        f"PRIMER_MIN_SIZE={config['primer_min_size']}",
        f"PRIMER_OPT_SIZE={config['primer_opt_size']}",
        f"PRIMER_MAX_SIZE={config['primer_max_size']}",
        f"PRIMER_MIN_TM={config['primer_min_tm']}",
        f"PRIMER_OPT_TM={config['primer_opt_tm']}",
        f"PRIMER_MAX_TM={config['primer_max_tm']}",
        f"PRIMER_MIN_GC={config['primer_min_gc']}",
        f"PRIMER_MAX_GC={config['primer_max_gc']}",
        f"PRIMER_GC_CLAMP={config.get('primer_gc_clamp', 1)}",
        f"PRIMER_MAX_POLY_X={config['primer_max_poly_x']}",
        f"PRIMER_MAX_SELF_ANY={config['primer_max_self_any']}",
        f"PRIMER_MAX_SELF_END={config['primer_max_self_end']}",
    ]
    if excluded:
        lines.append(f"SEQUENCE_EXCLUDED_REGION={excluded}")
    primer3_input = "\n".join(lines) + "\n=\n"

    result = subprocess.run(
        [primer3_core],
        input=primer3_input,
        capture_output=True,
        text=True,
        check=True,
    )
    data = {}
    for raw_line in result.stdout.splitlines():
        if "=" in raw_line:
            key, value = raw_line.split("=", 1)
            data[key] = value

    primers: list[dict] = []
    count = int(data.get("PRIMER_LEFT_NUM_RETURNED", "0") or "0")
    for idx in range(count):
        seq = data.get(f"PRIMER_LEFT_{idx}_SEQUENCE")
        if not seq:
            continue
        left_pos, left_len = _parse_pos_len(data[f"PRIMER_LEFT_{idx}"])
        cdna_start = template_start_cdna + left_pos
        cdna_end = cdna_start + left_len - 1
        primers.append(
            {
                "primer3_index": idx,
                "sequence": seq,
                "cdna_start": cdna_start,
                "cdna_end": cdna_end,
                "length": left_len,
                "tm": float(data.get(f"PRIMER_LEFT_{idx}_TM", "nan")),
                "gc": float(data.get(f"PRIMER_LEFT_{idx}_GC_PERCENT", "nan")),
                "penalty": float(data.get(f"PRIMER_LEFT_{idx}_PENALTY", "nan")),
                "self_any": _float_or_blank(data.get(f"PRIMER_LEFT_{idx}_SELF_ANY_TH")),
                "self_end": _float_or_blank(data.get(f"PRIMER_LEFT_{idx}_SELF_END_TH")),
                "hairpin": _float_or_blank(data.get(f"PRIMER_LEFT_{idx}_HAIRPIN_TH")),
            }
        )
    return primers


def _parse_pos_len(value: str) -> tuple[int, int]:
    pos, length = value.split(",", 1)
    return int(pos), int(length)


def _float_or_blank(value: str | None) -> float | str:
    if value in (None, ""):
        return ""
    try:
        return float(value)
    except ValueError:
        return ""

