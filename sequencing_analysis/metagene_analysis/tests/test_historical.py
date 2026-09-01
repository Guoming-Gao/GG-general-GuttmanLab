import pandas as pd

from spen_metagene.historical import summarize
from spen_metagene.historical_profiles import _bin, _gene_name, _window


def test_historical_ca_thresholds_and_ebs_discovery(tmp_path):
    ca=pd.DataFrame({
        "sample count":[5,4], "enrichment (intra)":[2.1,3],
        "enrichment (inter)":[2.1,3], "p-val (intra)":[.009,.001],
        "p-val (inter)":[.009,.001], "feature name":["NM_1|GeneA_intron","NM_2|GeneB_exon"],
    })
    ca.to_csv(tmp_path/"clap_endospen_1_CA.tsv",sep="\t",index=False)
    pd.DataFrame([["NM_1|GeneA_intron",1],["NM_2|GeneB_exon",2]]).to_csv(
        tmp_path/"clap_endospen_1_EBS.tsv",sep="\t",index=False,header=False)
    out=tmp_path/"summary.tsv"; summarize(str(tmp_path),str(out)); got=pd.read_csv(out,sep="\t")
    assert set(got.method)=={"CA","EBS"}
    assert got.loc[got.method.eq("CA"),"retained_rows"].item()==1
    assert got.loc[got.method.eq("EBS"),"permutations"].item()==10


def test_historical_peak_coordinate_parsing_and_orientation():
    assert _gene_name("NM_1|GeneA_first_intron".replace("_first",""))=="GeneA"
    assert _window("chr2:100-200")==('chr2',100,200)
    plus={"strand":"+","start":1000,"end":11000}
    minus={"strand":"-","start":1000,"end":11000}
    assert _bin(plus,1050,"fixed")==20
    assert _bin(minus,10949,"fixed")==20
