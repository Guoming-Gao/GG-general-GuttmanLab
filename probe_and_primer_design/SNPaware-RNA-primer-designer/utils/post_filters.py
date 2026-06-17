"""Post-Primer3 genotyping primer filters."""

from __future__ import annotations

import math
import re
import subprocess

from .sequence_utils import reverse_complement


TIER_ORDER = ["strict_all", "relaxed_amplicon", "relaxed_specificity", "relaxed_gc_tm"]


def annotate_post_filters(rows: list[dict], config: dict, ntthal: str) -> dict:
    """Annotate candidate rows with post-design filter metrics and select the best tier."""
    thermo_cache: dict[tuple[str, str, str | None], float | None] = {}
    tier_explanations = _tier_explanations(config)

    for row in rows:
        _annotate_row(row, config, ntthal, thermo_cache, tier_explanations)

    counts = {tier: sum(row["Filter_Tier"] == tier for row in rows) for tier in TIER_ORDER}
    counts["failed"] = sum(row["Filter_Tier"] == "failed" for row in rows)
    hard_pass_count = sum(row["Hard_Filter_Pass"] is True for row in rows)
    selected_tier = next((tier for tier in TIER_ORDER if counts[tier] > 0), "failed")
    for row in rows:
        row["Selected_Filter_Tier"] = selected_tier
        row["Selected_For_Top_Output"] = selected_tier != "failed" and row["Filter_Tier"] == selected_tier

    return {
        "selected_tier": selected_tier,
        "tier_counts": counts,
        "hard_pass_count": hard_pass_count,
        "hard_fail_count": len(rows) - hard_pass_count,
        "tier_explanations": tier_explanations,
    }


def annotate_post_filter_row(
    row: dict,
    config: dict,
    ntthal: str,
    thermo_cache: dict[tuple[str, str, str | None], float | None] | None = None,
) -> None:
    """Annotate one candidate row with the same config-aware post filters as batch annotation."""
    _annotate_row(row, config, ntthal, thermo_cache if thermo_cache is not None else {}, _tier_explanations(config))


def _tier_explanations(config: dict) -> dict[str, str]:
    strict_amplicon = f"{config['strict_amplicon_min']}-{config['strict_amplicon_max']} bp"
    configured_amplicon = f"{config['min_amplicon_size']}-{config['max_amplicon_size']} bp"
    strict_gc = f"{config['post_filter_gc_min']}-{config['post_filter_gc_max']}%"
    relaxed_gc = f"{config['relaxed_post_filter_gc_min']}-{config['relaxed_post_filter_gc_max']}%"
    strict_tm = f"{config['post_filter_tm_min']}-{config['post_filter_tm_max']} C"
    relaxed_tm = f"{config['relaxed_post_filter_tm_min']}-{config['relaxed_post_filter_tm_max']} C"
    min_snps = config["min_informative_snps"]
    return {
        "strict_all": (
            "All hard filters and soft preferences passed: strict BLAST specificity, "
            f"{strict_amplicon} cDNA amplicon, at least {min_snps} informative SNPs, "
            f"primer GC {strict_gc}, and primer Tm {strict_tm}."
        ),
        "relaxed_amplicon": (
            "Hard filters passed with strict BLAST specificity; only the strict cDNA "
            f"amplicon preference ({strict_amplicon}) was relaxed within the configured "
            f"allowed range ({configured_amplicon})."
        ),
        "relaxed_specificity": (
            "Hard filters passed within the configured cDNA amplicon range "
            f"({configured_amplicon}); strict BLAST specificity was relaxed."
        ),
        "relaxed_gc_tm": (
            "Last-resort rescue: primer GC/Tm was relaxed from "
            f"{strict_gc}/{strict_tm} to {relaxed_gc}/{relaxed_tm}, while the 3-prime "
            "G/C clamp and other hard safety filters still passed."
        ),
        "failed": "Not orderable: one or more hard selectable filters failed.",
    }


def _annotate_row(
    row: dict,
    config: dict,
    ntthal: str,
    thermo_cache: dict[tuple[str, str, str | None], float | None],
    tier_explanations: dict[str, str],
) -> None:
    left = str(row.get("Left_Primer_Seq", "") or "").upper()
    right = str(row.get("Right_Primer_Seq", "") or "").upper()
    hard_reasons: list[str] = []
    soft_reasons: list[str] = []

    left_len = len(left)
    right_len = len(right)
    left_gc = _float(row.get("Left_Primer_GC"))
    right_gc = _float(row.get("Right_Primer_GC"))
    left_tm = _float(row.get("Left_Primer_Tm"))
    right_tm = _float(row.get("Right_Primer_Tm"))
    tm_delta = _float(row.get("Primer_Tm_Delta"))
    amplicon_size = _int(row.get("Amplicon_Size"))
    snp_count = _int(row.get("SNP_Count_In_Amplicon"))
    min_snps = int(config["min_informative_snps"])

    left_hairpin = _ntthal_dg(left, None, "HAIRPIN", ntthal, thermo_cache)
    right_hairpin = _ntthal_dg(right, None, "HAIRPIN", ntthal, thermo_cache)
    left_self = _ntthal_dg(left, left, "ANY", ntthal, thermo_cache)
    right_self = _ntthal_dg(right, right, "ANY", ntthal, thermo_cache)
    hetero = _ntthal_dg(left, right, "ANY", ntthal, thermo_cache)
    thermo_values = [left_hairpin, right_hairpin, left_self, right_self, hetero]
    thermo_min = min([v for v in thermo_values if v is not None], default=None)

    valid_len = _in_range(left_len, config["valid_primer_length_min"], config["valid_primer_length_max"]) and _in_range(
        right_len, config["valid_primer_length_min"], config["valid_primer_length_max"]
    )
    gc_pass = _in_range(left_gc, config["post_filter_gc_min"], config["post_filter_gc_max"]) and _in_range(
        right_gc, config["post_filter_gc_min"], config["post_filter_gc_max"]
    )
    relaxed_gc_pass = _in_range(left_gc, config["relaxed_post_filter_gc_min"], config["relaxed_post_filter_gc_max"]) and _in_range(
        right_gc, config["relaxed_post_filter_gc_min"], config["relaxed_post_filter_gc_max"]
    )
    tm_pass = _in_range(left_tm, config["post_filter_tm_min"], config["post_filter_tm_max"]) and _in_range(
        right_tm, config["post_filter_tm_min"], config["post_filter_tm_max"]
    )
    relaxed_tm_pass = _in_range(left_tm, config["relaxed_post_filter_tm_min"], config["relaxed_post_filter_tm_max"]) and _in_range(
        right_tm, config["relaxed_post_filter_tm_min"], config["relaxed_post_filter_tm_max"]
    )
    tm_delta_pass = tm_delta <= float(config["strict_tm_delta_max"])
    gc_clamp_pass = _has_3p_gc_clamp(left) and _has_3p_gc_clamp(right)
    no_5p_g_pass = not left.startswith("G") and not right.startswith("G")
    homopolymer_pass = not _has_homopolymer(left, config["max_homopolymer_run"] + 1) and not _has_homopolymer(
        right, config["max_homopolymer_run"] + 1
    )
    quad_g_pass = "G" * (int(config["max_quad_g_run"]) + 1) not in left and "G" * (int(config["max_quad_g_run"]) + 1) not in right
    alternating_pass = not _has_alternating_dinuc(left, config["max_alternating_dinuc_bases"] + 1) and not _has_alternating_dinuc(
        right, config["max_alternating_dinuc_bases"] + 1
    )
    interprimer_run = _max_interprimer_complement_run(left, right)
    interprimer_run_pass = interprimer_run <= int(config["max_interprimer_complement_bases"])
    thermo_pass = all(v is not None and v > float(config["min_thermo_dg_kcal"]) for v in thermo_values)

    hard_chemistry = all(
        [
            tm_delta_pass,
            gc_clamp_pass,
            homopolymer_pass,
            quad_g_pass,
            alternating_pass,
            interprimer_run_pass,
            thermo_pass,
        ]
    )
    strict_chemistry = hard_chemistry and valid_len and gc_pass and tm_pass
    relaxed_gc_tm_chemistry = hard_chemistry and valid_len and relaxed_gc_pass and relaxed_tm_pass
    strict_specificity = row.get("Primer_Pair_Specificity_Pass") is True or row.get("Primer_Pair_Specificity_Pass") == "True"
    relaxed_specificity = row.get("Primer_Pair_Relaxed_Specificity_Pass") is True or row.get("Primer_Pair_Relaxed_Specificity_Pass") == "True"
    if not relaxed_specificity:
        relaxed_specificity = _int(row.get("Left_Primer_PerfectOrNear_Hits")) >= 1 and _int(row.get("Right_Primer_PerfectOrNear_Hits")) >= 1
    strict_amplicon = _in_range(amplicon_size, config["strict_amplicon_min"], config["strict_amplicon_max"])
    relaxed_amplicon = _in_range(amplicon_size, config["min_amplicon_size"], config["max_amplicon_size"])
    target_pass_field = config.get("target_pass_field", "Covers_All_Guides")
    target_failure_reason = config.get("target_failure_reason", "does_not_cover_all_guides")
    base_target_pass = (
        row.get(target_pass_field) is True
        and snp_count >= min_snps
        and row.get("Primer_Overlaps_SNP") is False
    )

    hard_filter_pass = base_target_pass and hard_chemistry and valid_len and relaxed_gc_pass and relaxed_tm_pass and relaxed_amplicon

    if not base_target_pass:
        if row.get(target_pass_field) is not True:
            hard_reasons.append(target_failure_reason)
        if snp_count < min_snps:
            hard_reasons.append(f"fewer_than_{min_snps}_snps")
        if row.get("Primer_Overlaps_SNP") is not False:
            hard_reasons.append("primer_overlaps_snp")
    if not valid_len:
        hard_reasons.append("length_outside_18_30")
    if not gc_pass:
        soft_reasons.append("gc_outside_40_60")
    if not relaxed_gc_pass:
        hard_reasons.append("gc_outside_35_65")
    if not tm_pass:
        soft_reasons.append("tm_outside_55_65")
    if not relaxed_tm_pass:
        hard_reasons.append("tm_outside_52_68")
    if not tm_delta_pass:
        hard_reasons.append("tm_delta_gt_5")
    if not gc_clamp_pass:
        hard_reasons.append("missing_3p_gc_clamp")
    if not no_5p_g_pass:
        soft_reasons.append("starts_with_5p_g")
    if not homopolymer_pass:
        hard_reasons.append("homopolymer_run_ge_4")
    if not quad_g_pass:
        hard_reasons.append("quad_g_ge_4")
    if not alternating_pass:
        hard_reasons.append("alternating_dinuc_repeat_ge_6_bases")
    if not interprimer_run_pass:
        hard_reasons.append("interprimer_complement_run_gt_3")
    if not thermo_pass:
        hard_reasons.append("thermo_dg_le_-9_or_unavailable")
    if not strict_specificity:
        soft_reasons.append("strict_specificity_failed")
    if not relaxed_specificity:
        soft_reasons.append("relaxed_specificity_failed")
    if not strict_amplicon:
        soft_reasons.append(f"outside_strict_amplicon_{config['strict_amplicon_min']}_{config['strict_amplicon_max']}")
    if not relaxed_amplicon:
        hard_reasons.append("outside_configured_amplicon_range")

    if base_target_pass and strict_chemistry and strict_specificity and strict_amplicon:
        tier = "strict_all"
        relaxed_criteria = ""
    elif base_target_pass and strict_chemistry and strict_specificity and relaxed_amplicon:
        tier = "relaxed_amplicon"
        relaxed_criteria = "amplicon"
    elif base_target_pass and strict_chemistry and relaxed_specificity and relaxed_amplicon:
        tier = "relaxed_specificity"
        relaxed_criteria = "specificity" + (",amplicon" if not strict_amplicon else "")
    elif base_target_pass and relaxed_gc_tm_chemistry and relaxed_specificity and relaxed_amplicon:
        tier = "relaxed_gc_tm"
        relaxed_criteria = ",".join(
            reason for reason in ["gc" if not gc_pass else "", "tm" if not tm_pass else "", "specificity" if not strict_specificity else "", "amplicon" if not strict_amplicon else ""] if reason
        )
    else:
        tier = "failed"
        relaxed_criteria = ""

    row.update(
        {
            "Left_Primer_Length": left_len,
            "Right_Primer_Length": right_len,
            "PostFilter_Length_18_30_Pass": valid_len,
            "PostFilter_GC_Pass": gc_pass,
            "PostFilter_Relaxed_GC_Pass": relaxed_gc_pass,
            "PostFilter_Tm_Pass": tm_pass,
            "PostFilter_Relaxed_Tm_Pass": relaxed_tm_pass,
            "PostFilter_Tm_Delta_Pass": tm_delta_pass,
            "PostFilter_3p_GC_Clamp_Pass": gc_clamp_pass,
            "PostFilter_No_5p_G_Pass": no_5p_g_pass,
            "PostFilter_Homopolymer_Pass": homopolymer_pass,
            "PostFilter_QuadG_Pass": quad_g_pass,
            "PostFilter_Alternating_Dinuc_Pass": alternating_pass,
            "PostFilter_Interprimer_Max_Complement_Run": interprimer_run,
            "PostFilter_Interprimer_Run_Pass": interprimer_run_pass,
            "Left_Hairpin_dG_kcal": _round_or_blank(left_hairpin),
            "Right_Hairpin_dG_kcal": _round_or_blank(right_hairpin),
            "Left_SelfDimer_dG_kcal": _round_or_blank(left_self),
            "Right_SelfDimer_dG_kcal": _round_or_blank(right_self),
            "Heterodimer_dG_kcal": _round_or_blank(hetero),
            "Worst_Thermo_dG_kcal": _round_or_blank(thermo_min),
            "PostFilter_Thermo_dG_Pass": thermo_pass,
            "Hard_Filter_Pass": hard_filter_pass,
            "PostFilter_Strict_Chemistry_Pass": strict_chemistry,
            "PostFilter_Relaxed_Chemistry_Pass": relaxed_gc_tm_chemistry,
            "PostFilter_Strict_Specificity_Pass": strict_specificity,
            "PostFilter_Relaxed_Specificity_Pass": relaxed_specificity,
            "PostFilter_Strict_Amplicon_Pass": strict_amplicon,
            "PostFilter_Relaxed_Amplicon_Pass": relaxed_amplicon,
            "PostFilter_Target_SNP_Pass": base_target_pass,
            "Filter_Tier": tier,
            "Soft_Filter_Tier": tier,
            "Selectable_For_Top_Output": tier != "failed",
            "Tier_Explanation": tier_explanations[tier],
            "Relaxed_Criteria_Used": relaxed_criteria,
            "Hard_Filter_Failure_Reasons": ";".join(dict.fromkeys(hard_reasons)),
            "Soft_Filter_Relaxed_Reasons": ";".join(dict.fromkeys(soft_reasons)),
            "PostFilter_Failure_Reasons": ";".join(dict.fromkeys(hard_reasons + soft_reasons)),
        }
    )


def _ntthal_dg(
    seq1: str,
    seq2: str | None,
    mode: str,
    ntthal: str,
    cache: dict[tuple[str, str, str | None], float | None],
) -> float | None:
    key = (mode, seq1, seq2)
    if key in cache:
        return cache[key]
    cmd = [ntthal, "-a", mode, "-s1", seq1]
    if seq2 is not None:
        cmd.extend(["-s2", seq2])
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    text = f"{result.stdout}\n{result.stderr}"
    if "No secondary structure could be calculated" in text:
        cache[key] = 0.0
        return 0.0
    match = re.search(r"dG\s*=\s*([-0-9.]+)", text)
    if not match:
        cache[key] = None
        return None
    cache[key] = float(match.group(1)) / 1000.0
    return cache[key]


def _has_3p_gc_clamp(seq: str) -> bool:
    return bool(seq) and seq[-1] in {"G", "C"}


def _has_homopolymer(seq: str, length: int) -> bool:
    return re.search(rf"([ACGT])\1{{{length - 1},}}", seq) is not None


def _has_alternating_dinuc(seq: str, min_bases: int) -> bool:
    repeats = math.ceil(min_bases / 2)
    for i in range(len(seq) - 1):
        motif = seq[i : i + 2]
        if len(set(motif)) != 2:
            continue
        if motif * repeats in seq:
            return True
    return False


def _max_interprimer_complement_run(left: str, right: str) -> int:
    rc_right = reverse_complement(right)
    best = 0
    for offset in range(-len(rc_right) + 1, len(left)):
        run = 0
        for i, base in enumerate(left):
            j = i - offset
            if 0 <= j < len(rc_right) and base == rc_right[j]:
                run += 1
                best = max(best, run)
            else:
                run = 0
    return best


def _in_range(value: float | int | None, lower: float | int, upper: float | int) -> bool:
    return value is not None and float(lower) <= float(value) <= float(upper)


def _float(value) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(parsed):
        return None
    return parsed


def _int(value) -> int:
    try:
        if value == "":
            return 0
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def _round_or_blank(value: float | None) -> float | str:
    if value is None:
        return ""
    return round(value, 3)
