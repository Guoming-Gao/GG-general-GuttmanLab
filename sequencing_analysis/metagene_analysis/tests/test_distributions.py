import numpy as np
import pandas as pd
from spen_metagene.distributions import build_distribution_data, _density

def test_density_integrates_to_one():
    grid=np.linspace(-5,5,1000); y=_density(np.linspace(-2,2,50),grid)
    assert np.isclose(np.trapezoid(y,grid),1,atol=1e-4)

def test_sparse_encoding_and_basemean(tmp_path):
    features=pd.DataFrame({"gene_id":["a","b","c"],"feature":["first_exon"]*3,"enrichment":[-1,0,1]})
    classes=pd.DataFrame({"gene_id":["a","b","c"],"gene_name":["A","B","C"],"baseMean":[0,9,99],"class_fdr05":["significant_decrease","no_change","no_change"]})
    fp=tmp_path/"features.tsv"; cp=tmp_path/"classes.tsv"; features.to_csv(fp,sep="\t",index=False); classes.to_csv(cp,sep="\t",index=False)
    d=build_distribution_data(fp,cp,min_kde_n=2)
    assert d.loc[d["class"].eq("significant_decrease"),"encoding"].eq("exact_impulse").all()
    base=d[d.value_type.eq("log10_basemean_plus_1")].sort_values("gene_id")
    assert np.allclose(base.value,[0,1,2])
    assert base.panel.eq("Expression").all()
    assert not d.panel.str.contains("_").any()
