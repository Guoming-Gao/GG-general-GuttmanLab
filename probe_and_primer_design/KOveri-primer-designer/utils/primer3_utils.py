"""Primer3 wrapper and parsing."""

from __future__ import annotations

import subprocess

from .models import SnpRecord


def build_excluded_regions(snps: list[SnpRecord], template_start: int, template_end: int) -> str:
    intervals = []
    for snp in snps:
        if template_start <= snp.pos <= template_end:
            intervals.append(f"{snp.pos - template_start},1")
    return " ".join(intervals)


def run_primer3_pair_design(
    sequence_id: str,
    template: str,
    target_start_abs: int,
    target_end_abs: int,
    template_start_abs: int,
    excluded_snps: list[SnpRecord],
    config: dict,
    primer3_core: str,
) -> list[dict]:
    """Run Primer3 pair design and return parsed candidate pairs."""
    target_rel_start = target_start_abs - template_start_abs
    target_len = target_end_abs - target_start_abs + 1
    product_range = f"{config['min_amplicon_size']}-{min(config['max_amplicon_size'], len(template))}"
    excluded = build_excluded_regions(excluded_snps, template_start_abs, template_start_abs + len(template) - 1)

    lines = [
        f"SEQUENCE_ID={sequence_id}",
        f"SEQUENCE_TEMPLATE={template}",
        f"SEQUENCE_TARGET={target_rel_start},{target_len}",
        "PRIMER_TASK=generic",
        "PRIMER_PICK_LEFT_PRIMER=1",
        "PRIMER_PICK_RIGHT_PRIMER=1",
        "PRIMER_PICK_INTERNAL_OLIGO=0",
        f"PRIMER_PRODUCT_SIZE_RANGE={product_range}",
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
        f"PRIMER_PAIR_MAX_COMPL_ANY={config['primer_pair_max_compl_any']}",
        f"PRIMER_PAIR_MAX_COMPL_END={config['primer_pair_max_compl_end']}",
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

    pairs: list[dict] = []
    pair_count = int(data.get("PRIMER_PAIR_NUM_RETURNED", "0") or "0")
    for idx in range(pair_count):
        left_seq = data.get(f"PRIMER_LEFT_{idx}_SEQUENCE")
        right_seq = data.get(f"PRIMER_RIGHT_{idx}_SEQUENCE")
        if not left_seq or not right_seq:
            continue
        left_pos, left_len = _parse_pos_len(data[f"PRIMER_LEFT_{idx}"])
        right_3p, right_len = _parse_pos_len(data[f"PRIMER_RIGHT_{idx}"])
        right_start = right_3p - right_len + 1
        pairs.append(
            {
                "primer3_index": idx,
                "left_sequence": left_seq,
                "right_sequence": right_seq,
                "left_rel_start": left_pos,
                "left_rel_end": left_pos + left_len - 1,
                "right_rel_start": right_start,
                "right_rel_end": right_3p,
                "left_genomic_start": template_start_abs + left_pos,
                "left_genomic_end": template_start_abs + left_pos + left_len - 1,
                "right_genomic_start": template_start_abs + right_start,
                "right_genomic_end": template_start_abs + right_3p,
                "amplicon_start": template_start_abs + left_pos,
                "amplicon_end": template_start_abs + right_3p,
                "amplicon_size": right_3p - left_pos + 1,
                "left_tm": float(data.get(f"PRIMER_LEFT_{idx}_TM", "nan")),
                "right_tm": float(data.get(f"PRIMER_RIGHT_{idx}_TM", "nan")),
                "left_gc": float(data.get(f"PRIMER_LEFT_{idx}_GC_PERCENT", "nan")),
                "right_gc": float(data.get(f"PRIMER_RIGHT_{idx}_GC_PERCENT", "nan")),
                "pair_penalty": float(data.get(f"PRIMER_PAIR_{idx}_PENALTY", "nan")),
                "left_penalty": float(data.get(f"PRIMER_LEFT_{idx}_PENALTY", "nan")),
                "right_penalty": float(data.get(f"PRIMER_RIGHT_{idx}_PENALTY", "nan")),
                "product_size": int(data.get(f"PRIMER_PAIR_{idx}_PRODUCT_SIZE", "0")),
            }
        )
    return pairs


def _parse_pos_len(value: str) -> tuple[int, int]:
    pos, length = value.split(",", 1)
    return int(pos), int(length)
