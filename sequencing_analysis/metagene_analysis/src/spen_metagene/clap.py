from __future__ import annotations

import json
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


FEATURES=["five_prime_utr","cds_exon","three_prime_utr","all_introns","first_intron","later_introns","first_exon","later_exons"]


def load_models(path):
    with open(path) as fh: return json.load(fh)


def build_index(models, flank=2000):
    per=defaultdict(list)
    for gid,t in models.items(): per[(t["chrom"],t["strand"])].append((max(0,t["start"]-flank),t["end"]+flank,gid,t["end"]-t["start"]))
    out={}
    for k,v in per.items():
        v.sort(); out[k]=([x[0] for x in v],[x[1] for x in v],[x[2] for x in v],[x[3] for x in v],max(x[1]-x[0] for x in v))
    return out


def assign(index, chrom, strand, pos):
    g=index.get((chrom,strand))
    if not g:return None
    starts,ends,ids,widths,maxlen=g; i=bisect_right(starts,pos)-1; floor=pos-maxlen; best=None; bw=None
    while i>=0 and starts[i]>=floor:
        if ends[i]>pos and (bw is None or widths[i]<bw): best,bw=ids[i],widths[i]
        i-=1
    return best


def _inside(pos, intervals): return any(a<=pos<b for a,b in intervals)


def feature_membership(t,pos):
    ex=sorted(t["exons"], reverse=t["strand"]=="-"); cds=sorted(t["cds"])
    out=[]
    exon_i=next((i for i,x in enumerate(ex) if x[0]<=pos<x[1]),None)
    if exon_i is not None:
        out.append("first_exon" if exon_i==0 else "later_exons")
        if _inside(pos,cds): out.append("cds_exon")
        elif cds:
            lo=min(x[0] for x in cds); hi=max(x[1] for x in cds)
            is5=(pos<lo) if t["strand"]=="+" else (pos>=hi)
            out.append("five_prime_utr" if is5 else "three_prime_utr")
    else:
        genomic=sorted(t["exons"]); introns=[[genomic[i][1],genomic[i+1][0]] for i in range(len(genomic)-1)]
        ordered=introns if t["strand"]=="+" else list(reversed(introns))
        ii=next((i for i,x in enumerate(ordered) if x[0]<=pos<x[1]),None)
        if ii is not None: out.extend(["all_introns","first_intron" if ii==0 else "later_introns"])
    return out


def count_sample(models_json,bam_path,sample,output_prefix,min_mapq=20,bin_size=100,tss_up=2000,tss_down=20000,scaled_bins=100,flank=2000):
    models=load_models(models_json); index=build_index(models,flank)
    fixed=defaultdict(lambda:np.zeros((tss_up+tss_down)//bin_size,dtype=np.int64))
    scaled=defaultdict(lambda:np.zeros(scaled_bins+2*(flank//bin_size),dtype=np.int64)); feature=defaultdict(lambda:defaultdict(int))
    metrics=defaultdict(int)
    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam.fetch(until_eof=True):
            metrics["records"]+=1
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality<min_mapq: continue
            if r.is_paired and not r.is_read2: continue
            if not r.reference_name or not r.reference_name.startswith("chr"): continue
            strand="-" if r.is_reverse else "+"
            if r.is_paired and r.template_length:
                start=min(r.reference_start,r.next_reference_start); pos=start+abs(r.template_length)//2
            else: pos=(r.reference_start+r.reference_end)//2
            gid=assign(index,r.reference_name,strand,pos)
            if gid is None: continue
            t=models[gid]; tss=t["start"] if strand=="+" else t["end"]-1; rel=pos-tss if strand=="+" else tss-pos
            if -tss_up<=rel<tss_down: fixed[gid][(rel+tss_up)//bin_size]+=1
            length=t["end"]-t["start"]
            if -flank<=rel<length+flank:
                if rel<0: b=(rel+flank)//bin_size
                elif rel<length: b=flank//bin_size+min(scaled_bins-1,int(rel/length*scaled_bins))
                else: b=flank//bin_size+scaled_bins+(rel-length)//bin_size
                if 0<=b<len(scaled[gid]): scaled[gid][b]+=1
            for f in feature_membership(t,pos): feature[gid][f]+=1
            metrics["assigned_fragments"]+=1
    prefix=Path(output_prefix); prefix.parent.mkdir(parents=True,exist_ok=True)
    base={g:{"gene_id":g,"gene_name":t["gene_name"],"chrom":t["chrom"],"strand":t["strand"],"gene_length":t["end"]-t["start"],"intron_count":max(0,len(t["exons"])-1)} for g,t in models.items()}
    for kind,data,n in [("fixed",fixed,(tss_up+tss_down)//bin_size),("scaled",scaled,scaled_bins+2*(flank//bin_size))]:
        rows=[]
        for gid,arr in data.items():
            rows.extend([{**base[gid],"bin":i,"count":int(v)} for i,v in enumerate(arr) if v])
        pd.DataFrame(rows).to_csv(str(prefix)+f".{kind}.tsv.gz",sep="\t",index=False,compression="gzip")
    rows=[]
    for gid,vals in feature.items():
        rows.extend([{**base[gid],"feature":f,"count":int(v)} for f,v in vals.items() if v])
    pd.DataFrame(rows).to_csv(str(prefix)+".features.tsv.gz",sep="\t",index=False,compression="gzip")
    with open(str(prefix)+".metrics.json","w") as fh: json.dump({"sample":sample,"bam":bam_path,**metrics},fh,indent=2)


def combine_sample_tables(inputs, samples, output):
    key=["gene_id","gene_name","chrom","strand","gene_length","intron_count","bin" if ".features." not in inputs[0] else "feature"]
    merged=None
    for f,s in zip(inputs,samples):
        d=pd.read_csv(f,sep="\t").rename(columns={"count":s})
        merged=d if merged is None else merged.merge(d,on=key,how="outer")
    merged.fillna(0).to_csv(output,sep="\t",index=False,compression="gzip")
