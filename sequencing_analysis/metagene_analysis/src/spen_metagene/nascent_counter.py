"""Audited strand-aware gene-body 4sU counter used by the prior analysis."""
from __future__ import annotations

import argparse, json, re, time
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path

import numpy as np
import pysam

COMP = str.maketrans("ACGTN", "TGCAN")


def load_gene_index(gtf, contig_regex=r"^chr"):
    crx=re.compile(contig_regex); genes=[]; seen=set(); per=defaultdict(list)
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"): continue
            f=line.split("\t")
            if len(f)<9 or f[2]!="gene" or not crx.match(f[0]): continue
            a=f[8]; gid=a.split('gene_id "',1)[1].split('"',1)[0]
            if gid in seen: continue
            seen.add(gid); name=a.split('gene_name "',1)[1].split('"',1)[0] if 'gene_name "' in a else gid
            idx=len(genes); genes.append((gid,name,f[0],f[6])); per[(f[0],f[6])].append((int(f[3])-1,int(f[4]),idx))
    out={}
    for k,v in per.items():
        v.sort(); starts=[x[0] for x in v]; ends=[x[1] for x in v]
        out[k]=(starts,ends,[x[2] for x in v],max(e-s for s,e in zip(starts,ends)))
    return out,genes


def assign(index, chrom, strand, pos):
    g=index.get((chrom,strand))
    if g is None: return -1
    starts,ends,ids,maxlen=g; i=bisect_right(starts,pos)-1; floor=pos-maxlen; best=-1; width=None
    while i>=0 and starts[i]>=floor:
        if ends[i]>pos and (width is None or ends[i]-starts[i]<width): best,width=ids[i],ends[i]-starts[i]
        i-=1
    return best


def count(bam_path, gtf, mask_path, label, output, min_bq=20, min_mapq=20, trim=5):
    t0=time.time(); index,genes=load_gene_index(gtf); z=np.load(mask_path)
    mask={k[6:]:np.sort(np.asarray(z[k],dtype=np.int64)) for k in z.files if k.startswith("cast::")}
    all_c=defaultdict(int); nas_c=defaultdict(int); pending={}; totals=defaultdict(int)
    def close(v):
        totals["fragments"]+=1; nas=v["tc"]>0; totals["nascent"]+=int(nas)
        if v["gene"]>=0:
            totals["assigned"]+=1; all_c[v["gene"]]+=1
            if nas: nas_c[v["gene"]]+=1; totals["assigned_nascent"]+=1
    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.mapping_quality<min_mapq or not r.is_paired or not r.is_proper_pair: continue
            if not r.reference_name or not r.reference_name.startswith("chr") or r.query_sequence is None or r.query_qualities is None: continue
            minus=(r.is_read2 and r.is_reverse) or (r.is_read1 and not r.is_reverse); strand="-" if minus else "+"
            try: pairs=r.get_aligned_pairs(matches_only=True,with_seq=True)
            except ValueError: continue
            m=mask.get(r.reference_name); masked=set() if m is None else set(m[np.searchsorted(m,r.reference_start):np.searchsorted(m,r.reference_end)].tolist())
            tc=0; qlen=r.query_length
            for qpos,rpos,refb in pairs:
                if refb is None or r.query_qualities[qpos]<min_bq or qpos<trim or qpos>=qlen-trim or rpos in masked: continue
                rb=refb.upper(); qb=r.query_sequence[qpos]
                if rb not in "ACGT" or qb not in "ACGT": continue
                if minus: rb,qb=rb.translate(COMP),qb.translate(COMP)
                tc += int(rb=="T" and qb=="C")
            gene=assign(index,r.reference_name,strand,(r.reference_start+r.reference_end)//2)
            v=pending.get(r.query_name)
            if v is None: pending[r.query_name]={"tc":tc,"gene":gene}
            else:
                v["tc"]+=tc
                if v["gene"]<0: v["gene"]=gene
                close(v); del pending[r.query_name]
    for v in pending.values(): close(v)
    Path(output).parent.mkdir(parents=True,exist_ok=True)
    with open(output,"w") as fh:
        fh.write("gene_id\tgene_name\tchrom\tstrand\tall_frag\tnascent_frag\n")
        for i,n in sorted(all_c.items()):
            gid,name,chrom,strand=genes[i]; fh.write(f"{gid}\t{name}\t{chrom}\t{strand}\t{n}\t{nas_c.get(i,0)}\n")
    meta={"label":label,"bam":bam_path,**totals,"seconds":round(time.time()-t0,1)}
    with open(str(output).replace(".tsv",".json"),"w") as fh: json.dump(meta,fh,indent=2)

