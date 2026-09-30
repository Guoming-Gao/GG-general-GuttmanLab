"""Check the SPEN design against an optional R output or source checkout.

Usage: python -m oligostan.validate_spen --spen-fasta /path/Spen_mm10.fa
       [--r-reference /path/Oligostan_Spen_ALL.tsv] [--report /path/report.json]
"""

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from .config import DEFAULT_SETTINGS
from .mouse_smifish import load_run_config
from .oligostan_core import get_probes_from_rna_dg37, process_probes_for_output

HUMAN_FIXTURE = Path(__file__).resolve().parent / "tests/data/humanRNU1_1.fa"
EXPECTED_HUMAN = {
    "AAAACCACCTTCGTGATCATGGTATCTCC", "GAGTGCAATGGATAAGCCTCGCCCTG",
    "TTTGGGGAAATCGCAGGGGTCAGCACATC", "CTACCACAAATTATGCAGTCGAGTTTCCCA",
    "CAGGGGAAAGCGCGAACGCAGTCCCC",
}


def design(fasta):
    records = list(SeqIO.parse(fasta, "fasta"))
    if len(records) != 1:
        raise ValueError("Expected exactly one FASTA record")
    antisense = str(records[0].seq.reverse_complement()).upper()
    params = {**DEFAULT_SETTINGS, "use_dustmasker": False}
    probes = get_probes_from_rna_dg37(
        antisense, params["taille_sonde_min"], params["taille_sonde_max"],
        params["fixed_dg37_value"], params["score_min"],
        params["distance_min_inter_sonde"],
    ) or []
    return process_probes_for_output(
        probes, {"name": Path(fasta).stem, "sequence": antisense},
        params["fixed_dg37_value"], **params,
    )


def source_design(fasta, source_checkout):
    script = '''import json,sys
from pathlib import Path
from Bio import SeqIO
sys.path.insert(0, sys.argv[1])
from config import DEFAULT_SETTINGS
from oligostan_core import get_probes_from_rna_dg37,process_probes_for_output
records=list(SeqIO.parse(sys.argv[2],"fasta"))
assert len(records)==1
seq=str(records[0].seq.reverse_complement()).upper()
p={**DEFAULT_SETTINGS,"use_dustmasker":False}
raw=get_probes_from_rna_dg37(seq,p["taille_sonde_min"],p["taille_sonde_max"],p["fixed_dg37_value"],p["score_min"],p["distance_min_inter_sonde"]) or []
print(json.dumps(process_probes_for_output(raw,{"name":Path(sys.argv[2]).stem,"sequence":seq},p["fixed_dg37_value"],**p)))'''
    result = subprocess.run(
        [sys.executable, "-c", script, str(source_checkout), str(fasta)],
        check=True, text=True, capture_output=True,
    )
    return json.loads(result.stdout)


def compare(left, right):
    mismatches = []
    if len(left) != len(right):
        mismatches.append({"field": "row_count", "left": len(left), "right": len(right)})
    for i, (a, b) in enumerate(zip(left, right), start=1):
        for key in sorted(set(a) & set(b)):
            x, y = a[key], b[key]
            if isinstance(x, (int, float)) and isinstance(y, (int, float)):
                equal = abs(x-y) <= 1e-9
            else:
                equal = x == y
            if not equal:
                mismatches.append({"row": i, "field": key, "left": x, "right": y})
                if len(mismatches) >= 30:
                    return mismatches
    return mismatches


def compare_r_table(python_rows, reference_path):
    """Compare R all-probes output by locus, ignoring presentation row order."""
    r = pd.read_csv(reference_path, sep="\t")
    p = pd.DataFrame(python_rows)
    keys = ["Seq", "theStartPos", "theEndPos", "ProbeSize"]
    fields = ["dGScore", "dG37", "GCpc", "GCFilter", "aCompFilter",
              "aStackFilter", "cCompFilter", "cStackFilter", "cSpecStackFilter",
              "NbOfPNAS", "PNASFilter", "InsideUTR", "HybFlpX", "HybFlpY", "HybFlpZ"]
    absent = [col for col in keys + fields if col not in r.columns or col not in p.columns]
    if absent:
        return {"match": False, "missing_columns": absent}
    joined = r[keys+fields].merge(p[keys+fields], on=keys, suffixes=("_r", "_python"),
                                  how="outer", indicator=True, validate="one_to_one")
    differences = []
    unmatched = joined[joined["_merge"] != "both"]
    for _, row in unmatched.head(30).iterrows():
        differences.append({"kind": str(row["_merge"]), "key": [row[k] for k in keys]})
    matched = joined[joined["_merge"] == "both"]
    for field in fields:
        a, b = matched[field+"_r"], matched[field+"_python"]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            bad = (a-b).abs() > 1e-9
        else:
            bad = a.astype(str) != b.astype(str)
        for _, row in matched[bad].head(max(0, 30-len(differences))).iterrows():
            differences.append({"field": field, "key": [row[k] for k in keys],
                                "r": row[field+"_r"], "python": row[field+"_python"]})
        if len(differences) >= 30:
            break
    r_filtered = int(((r.GCFilter == 1) & (r.PNASFilter == 1)).sum())
    p_filtered = int(((p.GCFilter == 1) & (p.PNASFilter == 1)).sum())
    return {"match": not differences and len(r)==len(p) and r_filtered==p_filtered,
            "r_count": len(r), "python_count": len(p),
            "r_filtered_count": r_filtered, "python_filtered_count": p_filtered,
            "differences": differences}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, help="Local JSON config with optional spen_fasta")
    parser.add_argument("--spen-fasta", type=Path, help="SPEN FASTA, or set spen_fasta in local config")
    parser.add_argument("--r-reference", type=Path)
    parser.add_argument("--source-checkout", type=Path, help="Optional second Python checkout for comparison")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    spen_fasta = args.spen_fasta
    if spen_fasta is None:
        spen_fasta = load_run_config(args.config).get("spen_fasta")
    if spen_fasta is None:
        parser.error("provide --spen-fasta or spen_fasta in the local config")
    human = design(HUMAN_FIXTURE)
    human_pass = len(human) == 5 and {p["Seq"] for p in human} == EXPECTED_HUMAN
    integrated = design(spen_fasta)
    source_diff = (compare(integrated, source_design(spen_fasta, args.source_checkout))
                   if args.source_checkout else None)
    report = {"human_fixture_pass": human_pass, "spen_probe_count": len(integrated),
              "source_checkout_match": None if not args.source_checkout else not source_diff,
              "source_differences": source_diff,
              "r_reference_checked": args.r_reference is not None}
    if args.r_reference:
        report["r_comparison"] = compare_r_table(integrated, args.r_reference)
        report["r_reference_sha256"] = hashlib.sha256(args.r_reference.read_bytes()).hexdigest()
    report["spen_fasta_sha256"] = hashlib.sha256(spen_fasta.read_bytes()).hexdigest()
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))
    if not human_pass or source_diff or (args.r_reference and not report["r_comparison"]["match"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
