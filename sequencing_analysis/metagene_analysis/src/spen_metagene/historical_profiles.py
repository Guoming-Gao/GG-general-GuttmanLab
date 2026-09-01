from __future__ import annotations

import json
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SAMPLES={
    "clap_endospen_1":("nodox","A"), "clap_endospen_2":("nodox","B"),
    "clap_doxendospen_1":("dox","A"), "clap_doxendospen_2":("dox","B"),
}
WINDOW=re.compile(r"^(chr[^:]+):(\d+)-(\d+)$")

def _gene_name(value):
    return re.sub(r"_(?:exon|intron)$","",str(value).split("|")[-1])

def _window(value):
    m=WINDOW.match(str(value))
    return None if m is None else (m.group(1),int(m.group(2)),int(m.group(3)))

def _read_calls(path,method):
    if method=="CA":
        use=["window","feature name","sample count","enrichment (intra)","enrichment (inter)","p-val (intra)","p-val (inter)"]; out=[]
        for d in pd.read_csv(path,sep="\t",usecols=use,chunksize=250000):
            d=d[(d["sample count"]>=5)&(d["enrichment (intra)"]>2)&(d["enrichment (inter)"]>2)&
                (d["p-val (intra)"]<.01)&(d["p-val (inter)"]<.01)&d["feature name"].ne("intergenic")]
            out.extend((_gene_name(r["feature name"]),_window(r["window"])) for _,r in d.iterrows())
        return out
    d=pd.read_csv(path,sep="\t",header=None,usecols=[0,2])
    return [(_gene_name(r.iloc[0]),_window(r.iloc[1])) for _,r in d.iterrows()]

def _bin(t,pos,geometry):
    tss=t["start"] if t["strand"]=="+" else t["end"]-1
    rel=pos-tss if t["strand"]=="+" else tss-pos
    if geometry=="fixed": return (rel+2000)//100 if -2000<=rel<20000 else None
    length=t["end"]-t["start"]
    if not (-2000<=rel<length+2000): return None
    if rel<0: return (rel+2000)//100
    if rel<length: return 20+min(99,int(rel/length*100))
    return 120+(rel-length)//100

def build(enrichment_dir,models_json,output_dir):
    out=Path(output_dir); out.mkdir(parents=True,exist_ok=True)
    models=json.loads(Path(models_json).read_text()); by_name={}
    for gid,t in sorted(models.items()): by_name.setdefault(t["gene_name"],(gid,t))
    sparse=[]; callsets={}; n_genes=len(models)
    for sample,(context,replicate) in SAMPLES.items():
        for method in ["CA","EBS"]:
            calls=_read_calls(Path(enrichment_dir)/f"{sample}_{method}.tsv",method); windows=set(); genes=set()
            for gene,w in calls:
                if w is None or gene not in by_name: continue
                gid,t=by_name[gene]; chrom,start,end=w
                if chrom!=t["chrom"]: continue
                pos=(start+end)//2; genes.add(gid); windows.add(f"{chrom}:{start}-{end}")
                for geometry in ["fixed","scaled"]:
                    b=_bin(t,pos,geometry)
                    if b is not None: sparse.append({"sample":sample,"context":context,"replicate":replicate,"method":method,"geometry":geometry,"gene_id":gid,"bin":int(b),"peak_count":1})
            callsets[(context,replicate,method)]={"genes":genes,"windows":windows}
    s=pd.DataFrame(sparse).groupby(["sample","context","replicate","method","geometry","gene_id","bin"],as_index=False).peak_count.sum()
    s.to_csv(out/"historical_peak_gene_bins.tsv.gz",sep="\t",index=False,compression="gzip")
    per=s.groupby(["sample","context","replicate","method","geometry","bin"],as_index=False).peak_count.sum(); per["peaks_per_gene"]=per.peak_count/n_genes
    curves=per.groupby(["context","method","geometry","bin"],as_index=False).agg(mean_peaks_per_gene=("peaks_per_gene","mean"),replicate_sd=("peaks_per_gene","std"))
    curves.to_csv(out/"historical_peak_density_curves.tsv",sep="\t",index=False)
    rows=[]
    for context in ["nodox","dox"]:
        for method in ["CA","EBS"]:
            a=callsets[(context,"A",method)]; b=callsets[(context,"B",method)]
            for level in ["genes","windows"]:
                union=a[level]|b[level]; inter=a[level]&b[level]
                rows.append({"context":context,"method":method,"level":level,"n_repA":len(a[level]),"n_repB":len(b[level]),"n_intersection":len(inter),"n_union":len(union),"jaccard":len(inter)/len(union) if union else np.nan})
    pd.DataFrame(rows).to_csv(out/"historical_replicate_overlap.tsv",sep="\t",index=False)
    fig,axes=plt.subplots(2,2,figsize=(13,8)); colors={"CA":"#0072B2","EBS":"#D55E00"}
    for i,context in enumerate(["nodox","dox"]):
        for j,geometry in enumerate(["fixed","scaled"]):
            ax=axes[i,j]
            for method in ["CA","EBS"]:
                z=curves[(curves.context==context)&(curves.method==method)&(curves.geometry==geometry)].sort_values("bin")
                q=z.bin*100-1950 if geometry=="fixed" else z.bin; ax.plot(q,z.mean_peaks_per_gene,color=colors[method],label=method)
            if geometry=="fixed": ax.axvline(0,color="black",ls=":",lw=.7); ax.set_xlabel("position relative to TSS (bp)")
            else:
                ax.axvline(20,color="black",ls=":",lw=.7); ax.axvline(120,color="black",ls=":",lw=.7)
                ax.set_xticks([0,20,70,120,140],["-2 kb","TSS","50%","TES","+2 kb"]); ax.set_xlabel("scaled transcript position")
            ax.set(ylabel="called peaks per annotated gene",title=f"{context} | {geometry}"); ax.legend(frameon=False)
    fig.suptitle("Historical CA/EBS local-hotspot sensitivity (not continuous occupancy)")
    fig.tight_layout(); fig.savefig(out/"historical_peak_density_metagene.png",dpi=250,bbox_inches="tight"); plt.close(fig)

def compare_stability(archived_dir,stability_dir,permutations,output):
    rows=[]
    for sample in SAMPLES:
        old=_read_calls(Path(archived_dir)/f"{sample}_EBS.tsv","EBS")
        new=_read_calls(Path(stability_dir)/f"{sample}.perm{permutations}.tsv","EBS")
        for level,fn in [("genes",lambda x:x[0]),("windows",lambda x:x[1])]:
            a={fn(x) for x in old if fn(x) is not None}; b={fn(x) for x in new if fn(x) is not None}; union=a|b; inter=a&b
            rows.append({"sample":sample,"level":level,"permutations_original":10,"permutations_stability":permutations,"n_original":len(a),"n_stability":len(b),"n_intersection":len(inter),"n_union":len(union),"jaccard":len(inter)/len(union) if union else np.nan})
    Path(output).parent.mkdir(parents=True,exist_ok=True); pd.DataFrame(rows).to_csv(output,sep="\t",index=False)
