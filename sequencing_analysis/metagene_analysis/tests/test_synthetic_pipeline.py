import pandas as pd

from spen_metagene.classes import CLASS_ORDER
from spen_metagene.suite import _plot_coordinates, make_suite


def test_synthetic_five_class_metagene_end_to_end(tmp_path):
    genes=[f"g{i}" for i in range(5)]
    classes=pd.DataFrame({
        "gene_id":genes,"gene_name":[x.upper() for x in genes],"chrom":["chr1"]*5,
        "baseMean":[10,20,30,40,50],"class_fdr05":CLASS_ORDER,"class_fdr10":CLASS_ORDER,
    })
    class_file=tmp_path/"classes.tsv"; classes.to_csv(class_file,sep="\t",index=False)
    common={"gene_length":25000,"intron_count":2,"chrom":"chr1"}
    fixed=[]; scaled=[]; features=[]
    for i,g in enumerate(genes):
        early=5.0 if g=="g4" else float(i)
        for b,v in [(20,early),(69,early),(70,0.0),(219,0.0)]:
            fixed.append({"gene_id":g,"bin":b,"enrichment":v,"clap_cpm":v+1,"input_cpm":1,**common})
        for b,v in [(20,early),(39,early),(40,0.0),(100,0.0),(119,0.0)]:
            scaled.append({"gene_id":g,"bin":b,"enrichment":v,"clap_cpm":v+1,"input_cpm":1,**common})
        for feature,v in [("first_exon",early),("later_exons",0),("first_intron",early),("later_introns",0),("all_introns",early)]:
            features.append({"gene_id":g,"feature":feature,"enrichment":v,"clap_cpm":v+1,"input_cpm":1,**common})
    signals={}
    for context in ["nodox","dox"]:
        for geometry,rows in [("fixed",fixed),("scaled",scaled),("features",features)]:
            p=tmp_path/f"{context}.{geometry}.tsv"; pd.DataFrame(rows).to_csv(p,sep="\t",index=False); signals[f"{context}_{geometry}"]=str(p)
    out=tmp_path/"out"
    make_suite(signals,{"Aux1h_vs_ctrl__nascent":str(class_file)},out,bootstrap=10,seed=9)
    manifest=pd.read_csv(out/"metagene_manifest.tsv",sep="\t")
    # FDR10 remains available in gene-class tables but is intentionally not
    # rendered; one contrast therefore yields fixed/scaled/features at FDR05.
    assert len(manifest)==3
    curve=pd.read_csv(out/"Aux1h_vs_ctrl__nascent__fdr05__all__scaled.tsv",sep="\t")
    elevated=curve[curve["class"].eq("significant_elevation")].set_index("bin")
    assert elevated.loc[20,"mean_enrichment"] > elevated.loc[100,"mean_enrichment"]
    assert _plot_coordinates("fixed",[0,20,219]).tolist()==[-1950.0,50.0,19950.0]
