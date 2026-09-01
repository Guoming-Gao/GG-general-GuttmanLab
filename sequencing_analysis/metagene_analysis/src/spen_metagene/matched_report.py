from __future__ import annotations
import base64, html
from pathlib import Path
import pandas as pd

def _img(path):
    p=Path(path)
    if not p.exists(): return f"<p class='missing'>Missing: {html.escape(str(p))}</p>"
    return f"<img src='data:image/png;base64,{base64.b64encode(p.read_bytes()).decode()}'>"

def _table(path,max_rows=30):
    p=Path(path)
    if not p.exists(): return f"<p class='missing'>Missing: {html.escape(str(p))}</p>"
    return pd.read_csv(p,sep="\t").head(max_rows).to_html(index=False,border=0,classes="data",float_format=lambda x:f"{x:.4g}")

def _figure(path,caption): return f"<figure>{_img(path)}<figcaption>{caption}</figcaption></figure>"

def build(output_root,output):
    root=Path(output_root); new=root/"figures/expression_matched"; old=root/"figures/feature_distributions"
    order=[("Aux2h_vs_ctrl","total","all"),("DoxAux24h_vs_Dox","total","autosomes"),("Aux2h_vs_ctrl","nascent","all"),("DoxAux24h_vs_Dox","nascent","autosomes"),("Aux1h_vs_ctrl","total","all"),("Aux1h_vs_ctrl","nascent","all")]
    sections=[]; effects=[]
    for c,r,s in order:
        stem=f"{c}__{r}__{s}__expression_matched"
        bal=new/f"{stem}__balance.tsv"; bins=new/f"{stem}__bin_diagnostics.tsv"; eff=new/f"{stem}__feature_effects.tsv"
        if eff.exists():
            q=pd.read_csv(eff,sep="\t"); q.insert(0,"analysis",f"{c} | {r} | {s}"); effects.append(q)
        caption=(f"{c}, {r}, {s}. The gray distribution contains only no-change genes selected without replacement in 0.1-wide log₁₀(baseMean+1) bins against the pooled directional target. "
                 "All directional genes are retained. In bins with too few no-change genes, every available control is used and residual imbalance remains. KDEs are independently normalized; sparse significant classes are exact impulses. Enrichment limits are identical to the unmatched figure.")
        sections.append(f"<h2>{c} — {r} RNA ({s})</h2><h3>Original unmatched distributions</h3>"+_figure(old/f"{c}__{r}__{s}__feature_distributions.png","Approved unmatched reference figure.")+"<h3>Expression-matched no-change reference</h3>"+_figure(new/f"{stem}.png",caption)+"<h3>Expression balance</h3>"+_table(bal)+"<h3>Bin availability and shortages</h3>"+_table(bins,200))
    all_eff=pd.concat(effects,ignore_index=True) if effects else pd.DataFrame()
    summary="<p>No feature-effect summaries were available.</p>"
    if len(all_eff):
        overview=all_eff.groupby(["analysis","direction","trend_status"]).size().rename("n_features").reset_index()
        summary=overview.to_html(index=False,border=0,classes="data")+"<h3>Per-feature median contrasts</h3>"+all_eff.to_html(index=False,border=0,classes="data",float_format=lambda x:f"{x:.4g}")
    xsections=[]
    for r in ["total","nascent"]:
        c="DoxAux24h_vs_Dox"; s="chrX"; stem=f"{c}__{r}__{s}__expression_matched"
        xsections.append(_figure(new/f"{stem}.png",f"Chromosome-X control, {r}. Matching uses the same pooled-direction, without-replacement rule; small classes and bin shortages make this exploratory."))
    css="body{font:15px system-ui;max-width:1280px;margin:auto;padding:30px;color:#222}h1,h2{color:#17365d}img{max-width:100%;border:1px solid #ddd}figure{margin:25px 0}figcaption{line-height:1.5}.note{background:#fff4cc;padding:14px;border-left:5px solid #e69f00}.data{border-collapse:collapse;font-size:12px;display:block;overflow-x:auto}.data th,.data td{border:1px solid #ddd;padding:4px 7px}.data th{background:#eef3f8}.missing{color:#a00}"
    intro="""<p>This companion analysis asks whether SPEN CLAP enrichment patterns persist after reducing the higher-expression bias of the no-change group. It does not redefine DE classes or create new discoveries.</p><div class='note'>All elevated and decreased classes are pooled only to define the target expression distribution. Directional genes are never downsampled. No-change genes are selected without replacement within fixed 0.1-wide log₁₀(baseMean+1) bins using seed 20260812. Where controls are insufficient, all available genes are retained and the mismatch is reported.</div><p>Persistence after matching supports an association not explained solely by measured baseline expression. Weakening or disappearance supports expression dependence. This procedure does not control gene length, intron architecture, or other covariates.</p>"""
    body=intro+"<h2>Cross-contrast trend summary</h2>"+summary+"".join(sections)+"<h2>Chromosome-X controls</h2>"+"".join(xsections)+f"<h2>Reproducibility</h2><p>Selected genes, bin shortages, balance diagnostics, vector figures, and plotted data are under <code>{html.escape(str(new))}</code>.</p>"
    doc=f"<!doctype html><html><head><meta charset='utf-8'><title>SPEN expression-matched enrichment analysis</title><style>{css}</style></head><body><h1>SPEN expression-matched enrichment analysis</h1>{body}</body></html>"
    p=Path(output); p.parent.mkdir(parents=True,exist_ok=True); p.write_text(doc)
