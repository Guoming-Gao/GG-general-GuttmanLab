#!/usr/bin/env python3
"""Demultiplex a single Nanopore FASTQ carrying Illumina inline dual indexes.

The script recognizes both read orientations from Illumina adapter anchors,
matches observed 8 bp P7/P5 barcode pairs to a sample map, and writes
per-sample FASTQs. Sample maps can use either:

    SampleName,P7_Barcode,P5_Barcode
    SampleName,i5,i7

For i5/i7 maps, i5 is treated as the P5-side barcode and i7 is reverse
complemented into the internal P7 orientation.
"""

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path

P7_PREFIX_FWD = "CAAGCAGAAGACGGCATACGAGAT"
P7_SUFFIX_FWD = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
P5_RC_PREFIX_FWD = "AGAGTGT"

P5_PREFIX_REV = "AATGATACGGCGACCACCGAGATCTACAC"
P5_SUFFIX_REV = "ACACTCTTTCCCT"
P7_RC_PREFIX_REV = "GAACTCCAGTCA"

BARCODE_LEN = 8
OFFSET_WINDOW = range(-2, 3)
RC_TRANS = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq):
    return seq.translate(RC_TRANS)[::-1].upper()


def open_text(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def sanitize_sample_name(name):
    name = name.strip()
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    return name.strip("._-")


def read_fastq(handle):
    while True:
        header = handle.readline()
        if not header:
            return
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not qual:
            raise ValueError("FASTQ ended mid-record")
        yield header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def write_fastq(handle, rec):
    header, seq, plus, qual = rec
    handle.write(f"{header}\n{seq}\n{plus}\n{qual}\n")


def hamming_or_large(a, b):
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(x != y for x, y in zip(a, b))


def candidate_barcodes(seq, expected_start, reverse=False):
    candidates = []
    for offset in OFFSET_WINDOW:
        start = expected_start + offset
        end = start + BARCODE_LEN
        if start < 0 or end > len(seq):
            continue
        raw = seq[start:end]
        barcode = reverse_complement(raw) if reverse else raw
        candidates.append({"barcode": barcode, "raw": raw, "offset": offset})
    return candidates


def best_candidate(candidates, expected):
    best = None
    for cand in candidates:
        dist = hamming_or_large(cand["barcode"], expected)
        item = (dist, abs(cand["offset"]), cand)
        if best is None or item[:2] < best[:2]:
            best = item
    if best is None:
        return None, None
    return best[0], best[2]


def orientation_candidates(seq):
    seq = seq.upper()
    out = []

    start = seq.find(P7_PREFIX_FWD)
    suffix = seq.find(P7_SUFFIX_FWD, start + len(P7_PREFIX_FWD)) if start != -1 else -1
    tail = seq.rfind(P5_RC_PREFIX_FWD)
    if start != -1 and suffix != -1 and tail != -1:
        p7_start = start + len(P7_PREFIX_FWD)
        p5_start = tail + len(P5_RC_PREFIX_FWD)
        out.append({
            "orientation": "forward",
            "p7_candidates": candidate_barcodes(seq, p7_start, reverse=False),
            "p5_candidates": candidate_barcodes(seq, p5_start, reverse=True),
        })

    start = seq.find(P5_PREFIX_REV)
    suffix = seq.find(P5_SUFFIX_REV, start + len(P5_PREFIX_REV)) if start != -1 else -1
    tail = seq.rfind(P7_RC_PREFIX_REV)
    if start != -1 and suffix != -1 and tail != -1:
        p5_start = start + len(P5_PREFIX_REV)
        p7_start = tail + len(P7_RC_PREFIX_REV)
        out.append({
            "orientation": "reverse",
            "p7_candidates": candidate_barcodes(seq, p7_start, reverse=True),
            "p5_candidates": candidate_barcodes(seq, p5_start, reverse=False),
        })

    return out


def load_sample_map(path, allow_blank=False):
    rows = []
    seen_samples = set()
    seen_pairs = set()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or [])
        has_p7_p5 = {"SampleName", "P7_Barcode", "P5_Barcode"}.issubset(fieldnames)
        has_i5_i7 = {"SampleName", "i5", "i7"}.issubset(fieldnames)
        if not has_p7_p5 and not has_i5_i7:
            raise ValueError(
                "Sample map must contain either SampleName,P7_Barcode,P5_Barcode "
                "or SampleName,i5,i7 columns"
            )

        for line_no, row in enumerate(reader, start=2):
            sample = sanitize_sample_name(row.get("SampleName", ""))
            if has_p7_p5:
                p7 = row.get("P7_Barcode", "").strip().upper()
                p5 = row.get("P5_Barcode", "").strip().upper()
            else:
                i5 = row.get("i5", "").strip().upper()
                i7 = row.get("i7", "").strip().upper()
                p7 = reverse_complement(i7) if i7 else ""
                p5 = i5

            if not p7 and not p5 and not sample:
                continue
            if len(p7) != BARCODE_LEN or len(p5) != BARCODE_LEN:
                label = "P7_Barcode/P5_Barcode" if has_p7_p5 else "i5/i7"
                raise ValueError(f"Line {line_no}: {label} must both be 8 bp")
            pair = (p7, p5)
            if pair in seen_pairs:
                raise ValueError(f"Line {line_no}: duplicate effective barcode pair P7={p7},P5={p5}")
            seen_pairs.add(pair)
            if sample:
                if sample in seen_samples:
                    raise ValueError(f"Line {line_no}: duplicate sample name after sanitizing: {sample}")
                seen_samples.add(sample)
            elif not allow_blank:
                raise ValueError(
                    f"Line {line_no}: SampleName is blank. Fill the sample map before demultiplexing, "
                    "or pass --allow-blank-samples to audit barcode detection only."
                )
            rows.append({"sample": sample, "p7": p7, "p5": p5})
    if not rows:
        raise ValueError("Sample map has no barcode rows")
    return rows


def match_read(seq, rows, max_mismatches):
    candidates = []
    orients = orientation_candidates(seq)
    for orient in orients:
        for row in rows:
            p7_dist, p7_cand = best_candidate(orient["p7_candidates"], row["p7"])
            p5_dist, p5_cand = best_candidate(orient["p5_candidates"], row["p5"])
            if p7_cand is None or p5_cand is None:
                continue
            total = p7_dist + p5_dist
            if total <= max_mismatches:
                candidates.append({
                    "distance": total,
                    "offset_penalty": abs(p7_cand["offset"]) + abs(p5_cand["offset"]),
                    "row": row,
                    "orientation": orient["orientation"],
                    "observed_p7": p7_cand["barcode"],
                    "observed_p5": p5_cand["barcode"],
                    "p7_offset": p7_cand["offset"],
                    "p5_offset": p5_cand["offset"],
                })

    if not candidates:
        return None, "unmatched_barcode" if orients else "unrecognized_structure"

    candidates.sort(key=lambda item: (item["distance"], item["offset_penalty"]))
    best = candidates[0]
    tied = [
        c for c in candidates[1:]
        if c["distance"] == best["distance"] and c["offset_penalty"] == best["offset_penalty"]
    ]
    if tied:
        return best, "ambiguous_barcode"
    return best, "matched"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq", required=True, help="Input multiplexed FASTQ or FASTQ.GZ")
    parser.add_argument("--sample-map", required=True, help="CSV with SampleName plus either i5/i7 or P7_Barcode/P5_Barcode")
    parser.add_argument("--output-dir", required=True, help="Directory for per-sample FASTQs")
    parser.add_argument("--summary", default=None, help="Summary CSV path; default is output-dir/demux_summary.csv")
    parser.add_argument("--unassigned-fastq", default=None, help="FASTQ for unassigned/ambiguous reads")
    parser.add_argument("--max-mismatches", type=int, default=1, help="Maximum total barcode mismatches across P7+P5 after local offset search")
    parser.add_argument("--allow-blank-samples", action="store_true", help="Audit barcode detection without writing blank sample outputs")
    args = parser.parse_args()

    rows = load_sample_map(args.sample_map, allow_blank=args.allow_blank_samples)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = Path(args.summary) if args.summary else out_dir / "demux_summary.csv"
    unassigned_path = Path(args.unassigned_fastq) if args.unassigned_fastq else out_dir / "unassigned.fastq"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    unassigned_path.parent.mkdir(parents=True, exist_ok=True)

    sample_handles = {}
    try:
        for row in rows:
            if not row["sample"]:
                continue
            sample_handles[row["sample"]] = open(out_dir / f"{row['sample']}.fastq", "w")

        unassigned_handle = open(unassigned_path, "w")
        counts = Counter()
        sample_counts = Counter()
        orientation_counts = Counter()
        mismatch_counts = Counter()
        canonical_pair_counts = Counter()
        observed_pair_counts = Counter()
        offset_counts = Counter()
        status_by_pair = defaultdict(Counter)

        with open_text(args.fastq, "rt") as handle:
            for rec in read_fastq(handle):
                counts["total_reads"] += 1
                best, status = match_read(rec[1], rows, args.max_mismatches)

                if status != "matched":
                    counts[status] += 1
                    write_fastq(unassigned_handle, rec)
                    if best:
                        canonical = (best["row"]["p7"], best["row"]["p5"])
                        status_by_pair[canonical][status] += 1
                    continue

                row = best["row"]
                orientation_counts[best["orientation"]] += 1
                mismatch_counts[best["distance"]] += 1
                canonical = (row["p7"], row["p5"])
                observed = (best["observed_p7"], best["observed_p5"])
                canonical_pair_counts[canonical] += 1
                observed_pair_counts[observed] += 1
                offset_counts[(best["p7_offset"], best["p5_offset"])] += 1
                status_by_pair[canonical][status] += 1

                if not row["sample"]:
                    counts["matched_blank_sample"] += 1
                    write_fastq(unassigned_handle, rec)
                    continue

                counts["assigned"] += 1
                sample_counts[row["sample"]] += 1
                write_fastq(sample_handles[row["sample"]], rec)
    finally:
        for handle in sample_handles.values():
            handle.close()
        if "unassigned_handle" in locals():
            unassigned_handle.close()

    with open(summary_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Section", "Key", "Value"])
        for key in [
            "total_reads",
            "assigned",
            "matched_blank_sample",
            "unrecognized_structure",
            "unmatched_barcode",
            "ambiguous_barcode",
        ]:
            writer.writerow(["overall", key, counts.get(key, 0)])
        for key, value in sorted(orientation_counts.items()):
            writer.writerow(["orientation", key, value])
        for key, value in sorted(mismatch_counts.items()):
            writer.writerow(["barcode_mismatches", key, value])
        for key, value in sorted(offset_counts.items()):
            writer.writerow(["barcode_offsets", f"P7={key[0]};P5={key[1]}", value])
        for sample, value in sorted(sample_counts.items()):
            writer.writerow(["sample", sample, value])
        for (p7, p5), value in canonical_pair_counts.most_common():
            statuses = ";".join(f"{k}:{v}" for k, v in sorted(status_by_pair[(p7, p5)].items()))
            writer.writerow(["canonical_pair", f"{p7},{p5},{statuses}", value])
        for (p7, p5), value in observed_pair_counts.most_common(100):
            writer.writerow(["observed_pair", f"{p7},{p5}", value])

    print(f"Total reads: {counts.get('total_reads', 0)}")
    print(f"Assigned reads: {counts.get('assigned', 0)}")
    print(f"Blank-sample matches: {counts.get('matched_blank_sample', 0)}")
    print(f"Unrecognized structure: {counts.get('unrecognized_structure', 0)}")
    print(f"Unmatched barcode: {counts.get('unmatched_barcode', 0)}")
    print(f"Ambiguous barcode: {counts.get('ambiguous_barcode', 0)}")
    print(f"Summary: {summary_path}")
    print(f"Unassigned/audit FASTQ: {unassigned_path}")


if __name__ == "__main__":
    main()
