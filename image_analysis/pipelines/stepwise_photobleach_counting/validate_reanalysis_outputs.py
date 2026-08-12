"""Acceptance checks for the consolidated cluster-aware PBSA result tree."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from cluster_reanalysis import STREAMS, load_reanalysis_config, roots


def validate(cfg: dict, require_drift: bool = False) -> dict[str, int | bool]:
    path = roots(cfg); ledger = pd.read_csv(path["tables"] / "complete_seed_ledger.csv")
    mapping = pd.read_csv(path["tables"] / "seed_to_analysis_unit_mapping.csv")
    assert len(ledger) == 14712, len(ledger)
    assert len(mapping) == 14712 and mapping.seed_level_unit_id.notna().all() and mapping.grouped_unit_id.notna().all()
    assert mapping[["dataset", "fov", "seed_id"]].duplicated().sum() == 0
    failed = ledger.legacy_growth_status.isin(["seed_below_local_threshold", "region_too_small"])
    assert ledger.loc[failed, "seed_level_unit_id"].notna().all()
    expected_counts = {"seed_level_cluster": 14711, "grouped_cluster": 14273, "extended_cluster_candidate": 917}
    filenames = {"seed_level_cluster": "seed_level_counts.csv", "grouped_cluster": "grouped_cluster_counts.csv", "extended_cluster_candidate": "extended_cluster_counts.csv"}
    for stream in STREAMS:
        table = pd.read_csv(path["tables"] / filenames[stream])
        assert len(table) == expected_counts[stream], (stream, len(table))
        assert table[["dataset", "fov", "analysis_unit_id"]].duplicated().sum() == 0
        assert {"pbsa_flag_meaning", "validated_range_status", "calibration_status"}.issubset(table.columns)
    thresholds = pd.read_csv(path["tables"] / "quickpbsa_thresholds.csv")
    ready = thresholds[(thresholds.selected_representation.astype(bool)) & (thresholds.calibration_status == "production_ready")]
    assert len(ready) == 8 and set(ready.analysis_stream) == {"seed_level_cluster", "grouped_cluster"}
    canonical = list(path["experiments"].glob("*/quickpbsa_counts/*/*/*__quickpbsa_counts.csv"))
    assert len(canonical) == 81, len(canonical)
    assert (path["report"] / "pipeline_audit_and_reanalysis_report.pdf").exists()
    assert len(list((path["figures"] / "seed_level_clusters").glob("*.png"))) >= 3
    assert len(list((path["figures"] / "grouped_clusters").glob("*.png"))) >= 1
    if require_drift:
        assert len(list(path["experiments"].glob("*/drift_correction/corrected_tiffs/*.tif"))) == 27
        assert len(list(path["experiments"].glob("*/input_QC/input_manifest.csv"))) == 2
        forbidden = [item for item in path["root"].rglob("*") if item.name in {".cache", "__pycache__", ".matplotlib"}]
        assert not forbidden, forbidden
        empty = [item for item in path["root"].rglob("*") if item.is_dir() and not any(item.iterdir())]
        assert not empty, empty
    result = {"seeds": len(ledger), "failed_growth_seeds_retained": int(failed.sum()),
              "grouped_traces": expected_counts["grouped_cluster"], "extended_traces": expected_counts["extended_cluster_candidate"],
              "production_calibrations": len(ready), "canonical_per_FOV_tables": len(canonical), "require_drift": require_drift}
    provenance = path["technical"] / "provenance"; provenance.mkdir(parents=True, exist_ok=True)
    (provenance / "acceptance_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--config", default="config_reanalysis.yaml"); parser.add_argument("--require-drift", action="store_true")
    args = parser.parse_args(); cfg = load_reanalysis_config(args.config); print(validate(cfg, args.require_drift))


if __name__ == "__main__": main()
