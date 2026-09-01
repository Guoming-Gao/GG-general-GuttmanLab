from __future__ import annotations
import base64, html
from pathlib import Path
import pandas as pd

def _img(path):
    p=Path(path)
    if not p.exists(): return f"<p class='missing'>Missing: {html.escape(str(p))}</p>"
    return f"<img src='data:image/png;base64,{base64.b64encode(p.read_bytes()).decode()}' alt='{html.escape(p.stem)}'>"

def _table(path,max_rows=100):
    p=Path(path)
    if not p.exists(): return f"<p class='missing'>Missing: {html.escape(str(p))}</p>"
    return pd.read_csv(p,sep="\t").head(max_rows).to_html(index=False,border=0,classes="data",float_format=lambda x:f"{x:.4g}")

def _figure(path,caption): return f"<figure>{_img(path)}<figcaption>{caption}</figcaption></figure>"

def _definitions():
    return """<h2>Definitions and interpretation</h2><dl>
<dt>DESeq2 MLE log₂ fold change</dt><dd>The ordinary, unregularized maximum-likelihood Aux-versus-control expression estimate. It defines the five primary gene classes using the effect boundary ±log₂(1.2).</dd>
<dt>apeglm log₂ fold change</dt><dd>An empirical-Bayes estimate based on a heavy-tailed prior. It pulls uncertain effects toward zero while retaining well-supported large effects. It is shown as an effect-size stability diagnostic and never changes the MLE-defined colors.</dd>
<dt>BH FDR and not tested</dt><dd>Benjamini–Hochberg adjusted p-values control the expected false-discovery proportion among called genes. Displayed significance is padj&lt;0.05. Genes lacking padj after independent filtering or insufficient information are retained as not tested and excluded from the five plotted biological classes. Nonsignificant directional classes are estimated directions, not discoveries or proof of no change.</dd>
<dt>SPEN CLAP enrichment over Input</dt><dd>For each matched replicate, enrichment is log₂[(TMM-normalized CLAP CPM + 0.1)/(TMM-normalized Input CPM + 0.1)]; the two replicate-specific log ratios are then averaged. Zero means equal normalized signal, +1 approximately twice as much CLAP as Input, and −1 approximately half as much. The 0.1-CPM pseudocount prevents undefined ratios and limits extreme low-count ratios; this is not a variance-stabilizing transform, shrinkage estimate, or peak call.</dd>
<dt>Equal-gene weighting</dt><dd>Each eligible gene supplies one observation per feature or positional bin. Class summaries average genes rather than reads, preventing long or deeply covered genes from automatically dominating.</dd>
<dt>Bootstrap interval</dt><dd>Scaled-profile 95% bands use 100 seeded resamples of whole genes. The displayed center uses all eligible genes. These intervals are exploratory.</dd></dl>"""

def _three_signal_explanation():
    return """<h2>Supplement: why three feature-level views are retained</h2>
<p><b>1. SPEN occupancy enrichment is the primary view.</b> The first panel shows the paired, pseudocount-adjusted log₂(CLAP/Input) ratio defined above. Input controls for RNA abundance, accessibility, fragmentation, mapping, and other background contributions. A positive ratio therefore asks whether SPEN CLAP is enriched beyond matched Input rather than merely whether reads were observed. The new KDE figures above are the preferred visualization of this signal because they show the full gene-level distribution.</p>
<p><b>2. CLAP and Input coverage density diagnoses the numerator and denominator.</b> A ratio can increase because CLAP rises, because Input falls, or both. The second panel therefore shows the two components separately. Counts are normalized with edgeR TMM effective library sizes, converted to CPM, and divided by the genomic width of each feature in kilobases, producing an FPKM-like coverage density. For compact comparison, log₂(1+density) CLAP is drawn above zero and Input is reflected below zero. Negative Input values are only a plotting convention; biological coverage is never negative. This component view should not be interpreted as another enrichment statistic.</p>
<p><b>3. Within-gene distribution of positive SPEN enrichment asks where binding lies.</b> Total binding magnitude can obscure location: a strongly bound gene can dominate a class even if its binding is not preferentially early. For each gene, negative enrichment values are first set to zero. Positive enrichment across the non-overlapping 5′ UTR, CDS exons, all introns, and 3′ UTR partition is then divided by that gene's summed positive enrichment, making the four fractions sum to one. Each gene consequently receives equal total mass. Genes with no positive enrichment have an undefined location distribution and are excluded only from this panel. First/later exon and intron categories overlap broader categories, so they remain useful in occupancy and coverage panels but are not included in the normalization denominator.</p>
<p>Together, the three views answer different questions: enrichment asks <i>whether CLAP exceeds Input</i>; component density asks <i>which assay drives the ratio</i>; within-gene allocation asks <i>where a gene's positive enrichment is located independent of magnitude</i>.</p>"""

def _shell(body,root):
    css="body{font:15px system-ui;max-width:1280px;margin:auto;padding:30px;color:#222}h1,h2{color:#17365d}figure{margin:30px 0}img{max-width:100%;height:auto;border:1px solid #ddd}figcaption{line-height:1.5;margin-top:9px}.note{background:#fff4cc;padding:14px;border-left:5px solid #e69f00}.data{border-collapse:collapse;font-size:12px;display:block;overflow-x:auto}.data th,.data td{border:1px solid #ddd;padding:4px 7px}.data th{background:#eef3f8}dt{font-weight:700;margin-top:11px}.missing{color:#a00}"
    return f"<!doctype html><html><head><meta charset='utf-8'><title>SPEN CLAP metagene analysis</title><style>{css}</style></head><body><h1>SPEN CLAP metagene analysis</h1>{body}<h2>Reproducibility</h2><p>Exact plotted data, vector figures, matrices, logs, tests, and provenance are under <code>{html.escape(str(root))}</code>.</p></body></html>"

def build_report(output_root,report_path,trimmed_path=None):
    root=Path(output_root); mf=pd.read_csv(root/"metagene/metagene_manifest.tsv",sep="\t"); counts=pd.read_csv(root/"tables/class_counts.tsv",sep="\t")
    standard=[("Aux1h_vs_ctrl","nascent"),("Aux1h_vs_ctrl","total"),("Aux2h_vs_ctrl","nascent"),("Aux2h_vs_ctrl","total"),("DoxAux24h_vs_Dox","nascent"),("DoxAux24h_vs_Dox","total")]
    priority=[("Aux2h_vs_ctrl","total"),("DoxAux24h_vs_Dox","total"),("Aux2h_vs_ctrl","nascent"),("DoxAux24h_vs_Dox","nascent"),("Aux1h_vs_ctrl","total"),("Aux1h_vs_ctrl","nascent")]
    def ns(c,r):
        q=counts[(counts.contrast==c)&(counts.readout==r)&(counts.threshold=="fdr05")]
        return ", ".join(f"{x['class'].replace('_',' ')}={int(x['n']):,}" for _,x in q.iterrows())
    intro="""<p><b>Primary significance:</b> BH FDR &lt;0.05. Early contrasts contain little FDR-supported differential expression; nonsignificant elevation and decrease groups describe MLE direction beyond ±log₂(1.2), not confirmed changes.</p><div class='note'>Only 1 h and 2 h define the early response. Dox+Aux 24 h versus Dox is independent, potentially includes secondary effects, and lacks a matched no-Dox 24 h condition. Its headline uses autosomes; chromosome X is treated separately as the Xist/SPEN control.</div>"""
    volcanos="".join(_figure(root/f"figures/volcano/{c}__{r}__fdr05.png",f"{c}, {r}. Left: DESeq2 MLE LFC; right: apeglm LFC. Both use −log₁₀(BH padj) and identical MLE/FDR05 colors. not-tested genes cannot be positioned on the padj axis. DE-table counts: {ns(c,r)}.") for c,r in standard)
    distributions=""
    for c,r in priority:
        subset="autosomes" if c=="DoxAux24h_vs_Dox" else "all"
        p=root/f"figures/feature_distributions/{c}__{r}__{subset}__feature_distributions.png"
        distributions+=_figure(p,f"{c}, {r}, {subset}. Rows 1–8 show independently normalized KDEs of gene-level SPEN CLAP enrichment over Input for classes with n≥30; smaller classes are exact impulses, not smoothed estimates. The eight panels share pooled 0.15th–99.94th-percentile limits, retaining approximately 99.8% of feature observations. Row 9 shows log₁₀(DESeq2 baseMean+1) for the same CLAP-eligible genes on a separate axis. Feature eligibility causes n to vary by row; exact per-row counts and observations are in the plotted TSV. Colors are MLE/FDR05 classes; directional nonsignificant curves are descriptive, not discoveries.")
    def meta(kind):
        out=""
        for c,r in standard:
            subset="autosomes" if c=="DoxAux24h_vs_Dox" else "all"; q=mf[(mf.contrast==c)&(mf.readout==r)&(mf.fdr=="fdr05")&(mf.subset==subset)&(mf.geometry==kind)]
            if len(q): out+=_figure(q.iloc[0].figure,f"{c}, {r}, {subset}. "+("Transcript-oriented 100-bin TSS-to-TES profile with 2-kb flanks. Genes are equally weighted; bands are exploratory 100-resample gene-bootstrap intervals." if kind=="scaled" else "Supplementary gene-level feature distributions for occupancy enrichment, component coverage density, and within-gene positive-enrichment allocation. See the preceding calculation rationale and exclusions."))
        return out
    agreement="".join(_figure(root/f"metagene/early_agreement__{r}__fdr05__all__scaled.png",f"Early 1 h/2 h MLE directional agreement, {r}; 24 h is excluded. Equal-gene SPEN CLAP enrichment over Input across scaled transcripts with exploratory 100-resample intervals. Directional agreement does not imply FDR significance.") for r in ["nascent","total"])
    xctrl=""
    for r in ["total","nascent"]:
        dp=root/f"figures/feature_distributions/DoxAux24h_vs_Dox__{r}__chrX__feature_distributions.png"
        xctrl+=_figure(dp,f"Dox+Aux 24 h versus Dox, {r}, chromosome X only. KDE and sparse-impulse rules match the autosomal headline; small chromosome-X classes require cautious interpretation.")
        q=mf[(mf.contrast=="DoxAux24h_vs_Dox")&(mf.readout==r)&(mf.fdr=="fdr05")&(mf.subset=="chrX")&(mf.geometry=="scaled")]
        if len(q): xctrl+=_figure(q.iloc[0].figure,f"Dox+Aux 24 h versus Dox, {r}, chromosome-X scaled positional control. This is an internal Xist/SPEN control, not a third early time point.")
    body=intro+_definitions()+"<h2>Expression contrasts and five MLE-defined classes</h2>"+volcanos+"<h2>Headline gene-level SPEN enrichment distributions</h2>"+distributions+"<h2>Scaled positional metagenes</h2>"+meta("scaled")+"<h2>Early 1 h/2 h agreement</h2>"+agreement+"<h2>Chromosome-X 24 h control</h2>"+xctrl+_three_signal_explanation()+meta("features")+"<h2>Execution audit</h2>"+_table(root/"tables/execution_audit.tsv")+"<h2>Input validation</h2>"+_table(root/"inputs_manifest/validation.tsv")+"<h2>Historical CA/EBS sensitivity</h2><p>Archived CA and EBS-10 calls are local-hotspot sensitivity analyses, not inputs to continuous occupancy. The redundant EBS-1000 rerun remains disabled.</p>"+_img(root/"historical/historical_peak_density_metagene.png")+"<h2>Limits</h2><p>KDEs are descriptive and independently normalized within class; curve height does not encode class abundance. Association between expression class and CLAP profile is not causal evidence. Sparse classes, especially chromosome-X and early significant groups, must not be interpreted as estimated population distributions.</p>"
    Path(report_path).parent.mkdir(parents=True,exist_ok=True); Path(report_path).write_text(_shell(body,root))
