from __future__ import annotations

import json
import re
from pathlib import Path
from collections import defaultdict
from dataclasses import asdict, dataclass, field


ATTR = re.compile(r'(\S+) "([^"]*)"')


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    start: int
    end: int
    gene_type: str = "unknown"
    transcript_type: str = "unknown"
    exons: list[list[int]] = field(default_factory=list)
    cds: list[list[int]] = field(default_factory=list)

    @property
    def length(self): return self.end - self.start


def _attrs(s: str) -> dict[str, str]:
    return dict(ATTR.findall(s))


def parse_representative_transcripts(gtf: str) -> dict[str, Transcript]:
    tx = {}
    with open(gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            f = line.rstrip().split("\t")
            if len(f) != 9 or not f[0].startswith("chr"): continue
            a = _attrs(f[8]); typ = f[2]
            tid, gid = a.get("transcript_id"), a.get("gene_id")
            if not tid or not gid: continue
            if tid not in tx:
                tx[tid] = Transcript(tid, gid, a.get("gene_name", gid), f[0], f[6], int(f[3])-1, int(f[4]), a.get("gene_type", "unknown"), a.get("transcript_type", "unknown"))
            t = tx[tid]
            if typ == "transcript": t.start, t.end = int(f[3])-1, int(f[4])
            elif typ == "exon": t.exons.append([int(f[3])-1, int(f[4])])
            elif typ == "CDS": t.cds.append([int(f[3])-1, int(f[4])])
    by_gene = defaultdict(list)
    for t in tx.values():
        if t.exons: by_gene[t.gene_id].append(t)
    reps = {}
    for gid, items in by_gene.items():
        coding = [t for t in items if t.transcript_type == "protein_coding" or t.cds]
        pool = coding or items
        reps[gid] = sorted(pool, key=lambda t: (-t.length, t.transcript_id))[0]
        reps[gid].exons.sort(); reps[gid].cds.sort()
    return reps


def save_models(gtf: str, output_json: str) -> None:
    reps = parse_representative_transcripts(gtf)
    with open(output_json, "w") as fh:
        json.dump({g: asdict(t) for g, t in reps.items()}, fh)


def save_gene_covariates(gtf: str, output_tsv: str) -> None:
    """Write gene-span and union-exon covariates without rebuilding BAM matrices."""
    import pandas as pd
    spans={}; exons=defaultdict(list)
    with open(gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            f=line.rstrip().split("\t")
            if len(f)!=9 or not f[0].startswith("chr") or f[2] not in {"gene","exon"}: continue
            a=_attrs(f[8]); gid=a.get("gene_id")
            if not gid: continue
            interval=(int(f[3])-1,int(f[4]))
            if f[2]=="gene": spans[gid]=(f[0],interval[0],interval[1])
            else: exons[gid].append(interval)
    rows=[]
    for gid,(chrom,start,end) in spans.items():
        merged=[]
        for s,e in sorted(exons.get(gid,[])):
            if merged and s<=merged[-1][1]: merged[-1][1]=max(merged[-1][1],e)
            else: merged.append([s,e])
        exon_bases=sum(e-s for s,e in merged); length=end-start
        rows.append({"gene_id":gid,"chrom":chrom,"gene_span":length,"union_exon_bases":exon_bases,
                     "intron_fraction":1-exon_bases/length if length else float("nan"),"union_exon_count":len(merged)})
    p=Path(output_tsv); p.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(p,sep="\t",index=False)
