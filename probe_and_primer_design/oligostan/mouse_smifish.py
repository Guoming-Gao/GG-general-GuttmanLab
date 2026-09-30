"""Mouse exon/intron smiFISH candidates and local BLAST verification.

Coordinates in exported tables are mm10/GRCm38, one-based and inclusive.
"""

from __future__ import annotations

import gzip
import json
import math
import re
import subprocess
import tempfile
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

from .config import DEFAULT_SETTINGS, FLAP_SEQUENCES
from .minimum_span import expected_sets, select_minimum_span_sets
from .oligostan_core import get_probes_from_rna_dg37, process_probes_for_output

GENES = ("Spen", "Arid5b", "Jarid2", "Sfmbt2", "Sfmbt1", "Mdm4", "Malat1")
DEFAULT_CONFIG = Path(__file__).resolve().parent.parent / ".smifish-local.json"
CONFIG_KEYS = ("gtf", "fasta", "blast_db", "blastn", "output_parent")
ATTR = re.compile(r'(\w+) "([^"]+)"')
MIN_PROBES_PER_SET = 30
QUALITY_TIERS = (
    {"name": "baseline", "min_gc": DEFAULT_SETTINGS["min_gc"],
     "max_gc": DEFAULT_SETTINGS["max_gc"],
     "rules": tuple(DEFAULT_SETTINGS["pnas_filter_option"])},
    {"name": "relax_a_content", "min_gc": 0.35, "max_gc": 0.65, "rules": (2, 4)},
    {"name": "relax_a_runs", "min_gc": 0.30, "max_gc": 0.70, "rules": (4,)},
    {"name": "gc_only", "min_gc": 0.30, "max_gc": 0.70, "rules": ()},
)


def progress():
    return Progress(SpinnerColumn(), TextColumn("{task.description}"), BarColumn(),
                    TextColumn("{task.completed}/{task.total}"), TimeElapsedColumn())


def load_run_config(config_path=None):
    """Load machine-specific reference paths from an untracked JSON file."""
    path = Path(config_path) if config_path is not None else DEFAULT_CONFIG
    if not path.is_file():
        raise FileNotFoundError(f"Create {path} using mouse_smifish.example.json as a template")
    raw = json.loads(path.read_text())
    missing = [key for key in CONFIG_KEYS if not isinstance(raw.get(key), str) or not raw[key].strip()]
    if missing:
        raise ValueError(f"Missing paths in {path}: {', '.join(missing)}")
    config = {key: Path(raw[key]).expanduser() for key in CONFIG_KEYS}
    if raw.get("spen_fasta"):
        config["spen_fasta"] = Path(raw["spen_fasta"]).expanduser()
    expected = {
        "gtf": config["gtf"], "fasta": config["fasta"],
        "fasta.fai": Path(str(config["fasta"]) + ".fai"),
        "blast_db.njs": Path(str(config["blast_db"]) + ".njs"),
        "blastn": config["blastn"], "output_parent": config["output_parent"],
    }
    if "spen_fasta" in config:
        expected["spen_fasta"] = config["spen_fasta"]
    absent = [f"{key}: {value}" for key, value in expected.items() if not value.exists()]
    if absent:
        raise FileNotFoundError("Missing configured inputs:\n" + "\n".join(absent))
    if not config["blastn"].is_file() or not config["output_parent"].is_dir():
        raise ValueError("blastn must be a file and output_parent must be a directory")
    return config


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(set(intervals)):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [tuple(x) for x in merged]


def subtract_intervals(parent, obstacles):
    pieces = []
    cursor = parent[0]
    for start, end in obstacles:
        if end < cursor or start > parent[1]:
            continue
        if start > cursor:
            pieces.append((cursor, min(start - 1, parent[1])))
        cursor = max(cursor, end + 1)
        if cursor > parent[1]:
            break
    if cursor <= parent[1]:
        pieces.append((cursor, parent[1]))
    return pieces


def load_models(gtf_path=None, genes=GENES):
    """Select the longest spliced RefGene transcript for each symbol.

    Intron intervals exclude every exon annotated for that gene, including
    exons from other isoforms. This prevents alternative exons entering an
    intron probe set.
    """
    gtf_path = Path(gtf_path) if gtf_path is not None else load_run_config()["gtf"]
    if not gtf_path.is_file():
        raise FileNotFoundError(gtf_path)
    wanted = {g.lower() for g in genes}
    raw = defaultdict(lambda: {"transcripts": defaultdict(list), "all_exons": []})
    with progress() as bar:
        task = bar.add_task("Reading mm10 GTF", total=1)
        opener = gzip.open if gtf_path.suffix == ".gz" else open
        with opener(gtf_path, "rt") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2] != "exon":
                    continue
                attrs = dict(ATTR.findall(fields[8]))
                symbol, tid = attrs.get("gene_name", ""), attrs.get("transcript_id", "")
                if symbol.lower() not in wanted or not tid:
                    continue
                item = raw[symbol]
                item["chrom"], item["strand"] = fields[0], fields[6]
                item["gene_id"] = attrs.get("gene_id", "")
                interval = (int(fields[3]), int(fields[4]))
                item["transcripts"][tid].append(interval)
                item["all_exons"].append(interval)
        bar.advance(task)
    models = {}
    for requested in genes:
        matches = [name for name in raw if name.lower() == requested.lower()]
        if len(matches) != 1:
            raise ValueError(f"Expected one GTF gene for {requested}, found {matches}")
        symbol = matches[0]
        item = raw[symbol]
        transcripts = {tid: merge_intervals(exons) for tid, exons in item["transcripts"].items()}
        tid = max(transcripts, key=lambda t: (sum(e-s+1 for s,e in transcripts[t]), t))
        exons = transcripts[tid]
        all_exons = merge_intervals(item["all_exons"])
        introns = []
        for left, right in zip(exons, exons[1:]):
            if left[1]+1 <= right[0]-1:
                introns.extend(subtract_intervals((left[1]+1, right[0]-1), all_exons))
        models[symbol] = {
            "gene": symbol, "gene_id": item["gene_id"], "transcript_id": tid,
            "chrom": item["chrom"], "strand": item["strand"],
            "exon": exons, "intron": [x for x in introns if x[1]-x[0]+1 >= 32],
            "all_gene_exons": all_exons,
        }
    return models


class FastaIndex:
    def __init__(self, fasta=None):
        self.fasta = Path(fasta) if fasta is not None else load_run_config()["fasta"]
        self.fai = Path(str(self.fasta) + ".fai")
        if not self.fasta.is_file() or not self.fai.is_file():
            raise FileNotFoundError(f"FASTA or .fai missing: {self.fasta}")
        self.index = {}
        for line in self.fai.read_text().splitlines():
            name, length, offset, line_bases, line_width, *_ = line.split("\t")
            self.index[name] = tuple(map(int, (length, offset, line_bases, line_width)))

    def fetch(self, chrom, start, end):
        """Get reference plus-strand bases at one-based inclusive positions."""
        length, offset, line_bases, line_width = self.index[chrom]
        if not 1 <= start <= end <= length:
            raise ValueError((chrom, start, end, length))
        first, last = start - 1, end - 1
        byte_start = offset + (first // line_bases) * line_width + first % line_bases
        byte_end = offset + (last // line_bases) * line_width + last % line_bases
        with self.fasta.open("rb") as handle:
            handle.seek(byte_start)
            sequence = handle.read(byte_end - byte_start + 1).replace(b"\n", b"").replace(b"\r", b"").decode().upper()
        if len(sequence) != end-start+1:
            raise RuntimeError(f"FASTA index length mismatch at {chrom}:{start}-{end}")
        return sequence


def tile_intervals(intervals, tile_size=4000):
    """Interleave tiles across regions so a long first intron cannot dominate."""
    groups = []
    for start, end in intervals:
        tiles = [(s, min(s+tile_size-1, end)) for s in range(start, end+1, tile_size)]
        groups.append([x for x in tiles if x[1]-x[0]+1 >= 32])
    return [tile for depth in range(max(map(len, groups), default=0))
            for group in groups if depth < len(group)
            for tile in [group[depth]]]


def quality_tier(metric):
    flags = {1: metric["aCompFilter"], 2: metric["aStackFilter"],
             3: metric["cCompFilter"], 4: metric["cStackFilter"],
             5: metric["cSpecStackFilter"]}
    for index, tier in enumerate(QUALITY_TIERS):
        if tier["min_gc"] <= metric["GCpc"] <= tier["max_gc"] and all(flags[r] for r in tier["rules"]):
            return index
    return None


def generate_candidates(models, fasta=None, intron_pool=None):
    """Design across all eligible tiles; an optional cap is exploratory only."""
    reference = FastaIndex(fasta)
    params = {**DEFAULT_SETTINGS, "use_dustmasker": False}
    rows = []
    jobs = [(m, region) for m in models.values() for region in ("exon", "intron") if m[region]]
    total_tiles = sum(len(tile_intervals(m[region])) for m, region in jobs)
    with progress() as bar:
        task = bar.add_task("Designing exon and intron candidates", total=total_tiles)
        for model, region in jobs:
            count = 0
            bar.update(task, description=f"Designing {model['gene']} {region}")
            for tile_start, tile_end in tile_intervals(model[region]):
                plus = reference.fetch(model["chrom"], tile_start, tile_end)
                target = plus if model["strand"] == "+" else str(Seq(plus).reverse_complement())
                if set(target) - set("ACGT"):
                    # Ambiguous bases cannot be ordered or BLAST verified reliably.
                    continue
                antisense = str(Seq(target).reverse_complement())
                probes = get_probes_from_rna_dg37(
                    antisense,
                    min_size_probe=params["taille_sonde_min"],
                    max_size_probe=params["taille_sonde_max"],
                    desired_dg=params["fixed_dg37_value"],
                    min_score_value=params["score_min"],
                    inc_betw_prob=params["distance_min_inter_sonde"],
                ) or []
                metrics = process_probes_for_output(
                    probes, {"name": model["gene"], "sequence": antisense},
                    params["fixed_dg37_value"], **params,
                )
                for (size, score, pos, oligo), metric in zip(probes, metrics):
                    tier = quality_tier(metric)
                    if tier is None:
                        continue
                    target_offset = len(target) - pos - size + 1  # zero-based
                    if model["strand"] == "+":
                        start = tile_start + target_offset
                        end = start + size - 1
                    else:
                        end = tile_end - target_offset
                        start = end - size + 1
                    target_seq = target[target_offset:target_offset+size]
                    if oligo != str(Seq(target_seq).reverse_complement()):
                        raise AssertionError((model["gene"], region, start, end))
                    rows.append({
                        "gene": model["gene"], "gene_id": model["gene_id"],
                        "transcript_id": model["transcript_id"], "region": region,
                        "genome_build": "mm10/GRCm38", "chrom": model["chrom"],
                        "start": start, "end": end, "target_strand": model["strand"],
                        "target_seq": target_seq, "probe_seq": oligo,
                        "probe_length": size, "dGScore": score,
                        "dG37": metric["dG37"], "GCpc": metric["GCpc"],
                        "GCFilter": metric["GCFilter"], "PNASFilter": metric["PNASFilter"],
                        "aCompFilter": metric["aCompFilter"],
                        "aStackFilter": metric["aStackFilter"],
                        "cCompFilter": metric["cCompFilter"],
                        "cStackFilter": metric["cStackFilter"],
                        "cSpecStackFilter": metric["cSpecStackFilter"],
                        "NbOfPNAS": metric["NbOfPNAS"],
                        "quality_tier": tier,
                        "quality_tier_name": QUALITY_TIERS[tier]["name"],
                        "HybFlpX": oligo + FLAP_SEQUENCES["X"],
                        "HybFlpY": oligo + FLAP_SEQUENCES["Y"],
                        "HybFlpZ": oligo + FLAP_SEQUENCES["Z"],
                    })
                    count += 1
                bar.advance(task)
                if region == "intron" and intron_pool is not None and count >= intron_pool:
                    break
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise RuntimeError("No candidates passed even the loosest configured quality tier")
    # Repeated oligos across tiles or loci cannot be uniquely assigned.
    frame = frame.drop_duplicates(subset=["probe_seq"], keep=False).reset_index(drop=True)
    frame.insert(0, "probe_id", [f"q{i:06d}" for i in range(1, len(frame)+1)])
    return frame


BLAST_COLUMNS = ("probe_id", "subject", "identity", "alignment_length", "mismatches",
                 "gap_opens", "query_start", "query_end", "subject_start", "subject_end",
                 "evalue", "bitscore", "query_length")


def blast_verify(candidates, blastn=None, database=None,
                 threads=8, min_offtarget_coverage=0.8, min_offtarget_identity=90.0,
                 show_progress=True):
    """Require a full exact intended hit and no strong off-target BLAST hit.

    A strong off-target covers >=80% of the probe at >=90% identity. Shorter or
    more divergent matches remain in the audit table but do not disqualify a probe.
    BLAST searches both strands with dust and soft masking disabled.
    """
    if blastn is None or database is None:
        config = load_run_config()
        blastn = config["blastn"] if blastn is None else blastn
        database = config["blast_db"] if database is None else database
    blastn, database = Path(blastn), Path(database)
    if not blastn.is_file() or not Path(str(database)+".njs").is_file():
        raise FileNotFoundError(f"BLAST executable/database missing: {blastn}, {database}")
    with tempfile.TemporaryDirectory() as tmp:
        input_fasta, output_tsv = Path(tmp)/"queries.fa", Path(tmp)/"hits.tsv"
        with input_fasta.open("w") as out:
            for row in candidates.itertuples(index=False):
                out.write(f">{row.probe_id}\n{row.probe_seq}\n")
        command = [str(blastn), "-task", "blastn-short", "-query", str(input_fasta),
                   "-db", str(database), "-out", str(output_tsv), "-outfmt",
                   "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen",
                   "-word_size", "7", "-dust", "no", "-soft_masking", "false",
                   "-evalue", "1", "-max_target_seqs", "500", "-max_hsps", "50",
                   "-num_threads", str(threads)]
        if show_progress:
            with progress() as bar:
                task = bar.add_task("BLAST verifying candidates", total=1)
                result = subprocess.run(command, text=True, capture_output=True)
                bar.advance(task)
        else:
            result = subprocess.run(command, text=True, capture_output=True)
        if result.returncode:
            raise RuntimeError(f"blastn failed: {result.stderr[-2000:]}")
        if output_tsv.stat().st_size:
            hits = pd.read_csv(output_tsv, sep="\t", names=BLAST_COLUMNS)
        else:
            hits = pd.DataFrame(columns=BLAST_COLUMNS)
    if not hits.empty:
        for col in ("alignment_length", "query_start", "query_end", "subject_start", "subject_end", "query_length"):
            hits[col] = hits[col].astype(int)
        hits = hits.merge(
            candidates[["probe_id", "chrom", "start", "end", "target_strand"]],
            on="probe_id", how="left", validate="many_to_one",
        )
        hits["expected_locus"] = (
            (hits["subject"] == hits["chrom"]) &
            (hits[["subject_start", "subject_end"]].min(axis=1) == hits["start"]) &
            (hits[["subject_start", "subject_end"]].max(axis=1) == hits["end"]) &
            ((hits["subject_start"] < hits["subject_end"]) != (hits["target_strand"] == "+"))
        )
        hits["full_exact"] = (
            (hits["identity"] == 100.0) &
            (hits["alignment_length"] == hits["query_length"]) &
            (hits["query_start"] == 1) &
            (hits["query_end"] == hits["query_length"])
        )
        hits["strong_offtarget"] = (
            ~hits["expected_locus"] &
            (hits["alignment_length"] >= hits["query_length"].map(lambda n: math.ceil(n*min_offtarget_coverage))) &
            (hits["identity"] >= min_offtarget_identity)
        )
    else:
        hits["expected_locus"] = pd.Series(dtype=bool)
        hits["full_exact"] = pd.Series(dtype=bool)
        hits["strong_offtarget"] = pd.Series(dtype=bool)
    passing_target = set(hits.loc[hits.expected_locus & hits.full_exact, "probe_id"])
    failing_offtarget = set(hits.loc[hits.strong_offtarget, "probe_id"])
    verified = candidates.copy()
    verified["blast_target_exact"] = verified.probe_id.isin(passing_target)
    verified["blast_strong_offtarget"] = verified.probe_id.isin(failing_offtarget)
    verified["blast_verified"] = verified.blast_target_exact & ~verified.blast_strong_offtarget
    verified["blast_hit_count"] = verified.probe_id.map(hits.groupby("probe_id").size()).fillna(0).astype(int)
    return verified, hits, command


def adaptive_blast_verify(candidates, batch_size=3000, **kwargs):
    """BLAST every filtered candidate in bounded batches before selection.

    The name is retained for existing callers; there is no tier-based early stop.
    """
    if batch_size < 1 or candidates.empty:
        raise ValueError("BLAST needs candidates and a positive batch_size")
    all_verified, all_hits, commands = [], [], []
    with progress() as bar:
        task = bar.add_task("BLAST verifying complete candidate pool", total=len(candidates))
        for start in range(0, len(candidates), batch_size):
            batch = candidates.iloc[start:start + batch_size]
            verified, hits, command = blast_verify(batch, show_progress=False, **kwargs)
            all_verified.append(verified)
            all_hits.append(hits)
            commands.append(command)
            bar.advance(task, len(batch))
    return pd.concat(all_verified, ignore_index=True), pd.concat(all_hits, ignore_index=True), commands


def select_sets(verified, models, minimum=MIN_PROBES_PER_SET):
    """Compatibility name for the BLAST-first minimum-span add-on."""
    return select_minimum_span_sets(verified, models, probes_per_set=minimum)


def write_outputs(models, candidates, verified, hits, selected, summary, blast_command,
                  output_parent=None, run_name=None, reference_paths=None,
                  intron_pool=None):
    if output_parent is None or reference_paths is None:
        config = load_run_config()
        output_parent = config["output_parent"] if output_parent is None else output_parent
        reference_paths = config if reference_paths is None else reference_paths
    parent = Path(output_parent)
    if not parent.is_dir():
        raise FileNotFoundError(f"Output parent does not exist: {parent}")
    name = run_name or f"SPEN_targets_with_Malat1_control_mm10_{datetime.now():%Y%m%d_%H%M%S}"
    root = parent / name
    root.mkdir(exist_ok=False)
    with progress() as bar:
        task = bar.add_task("Writing audited probe sets", total=len(summary))
        candidates.to_csv(root/"all_design_candidates.csv", index=False)
        verified.to_csv(root/"all_candidates_blast_status.csv", index=False)
        hits.to_csv(root/"blast_hits.tsv", sep="\t", index=False)
        selected.to_csv(root/"all_selected_blast_verified.csv", index=False)
        summary.to_csv(root/"set_summary.csv", index=False)
        for gene, region in expected_sets(models):
            subset = selected[(selected.gene == gene) & (selected.region == region)].sort_values("start")
            subdir = root / gene / region
            subdir.mkdir(parents=True)
            subset.to_csv(subdir/"selected_blast_verified.csv", index=False)
            with (subdir/"selected_probes.fa").open("w") as out:
                for row in subset.itertuples(index=False):
                    out.write(f">{row.probe_id} {gene} {region} {row.chrom}:{row.start}-{row.end}\n{row.probe_seq}\n")
            pd.DataFrame({"Name": [f"{gene}_{region}_{i:02d}" for i in range(1,len(subset)+1)],
                          "Sequence": subset.HybFlpX.to_list(), "Scale": ["25nm"]*len(subset),
                          "Purification": ["STD"]*len(subset)}).to_csv(subdir/"order_FlapX_REVIEW.csv", index=False)
            bar.advance(task)
    manifest = {
        "status": "BLAST verified Python Oligostan candidates; full mouse Spen R parity passed for matched fixed -32 and masking-off settings",
        "genome_build": "mm10/GRCm38", "gtf": str(reference_paths["gtf"]),
        "fasta": str(reference_paths["fasta"]), "blast_database": str(reference_paths["blast_db"]),
        "blast_command": blast_command,
        "blast_acceptance": "100% full-length exact intended genomic hit and no other BLAST alignment covering >=80% of the probe at >=90% identity; blastn-short word size 7, dust off, E-value 1, up to 500 subjects and 50 HSPs per subject",
        "quality_tiers": QUALITY_TIERS,
        "selection_mode": "BLAST-first minimum genomic span of 30 probes from the full verified pool",
        "candidate_pool_complete": intron_pool is None,
        "intron_pool": intron_pool,
        "candidate_count": len(candidates), "blast_tested_count": len(verified),
        "blast_verified_count": int(verified.blast_verified.sum()),
        "selected_probe_count": len(selected), "set_count": len(summary),
        "source_transcripts": {gene: model["transcript_id"] for gene,model in models.items()},
        "reference": "https://bitbucket.org/muellerflorian/fish_quant/src/master/Oligostan/Oligostan.r",
    }
    (root/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    from .coverage_report import write_coverage_report
    write_coverage_report(root, models, selected, summary)
    return root
