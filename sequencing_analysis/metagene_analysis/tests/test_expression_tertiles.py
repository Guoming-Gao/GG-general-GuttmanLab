import numpy as np
import pandas as pd

from spen_metagene.expression_tertiles import assign_tertiles, baseline_expression, eligible_gene_ids


def test_baseline_uses_only_requested_untreated_columns(tmp_path):
    p=tmp_path/"norm.tsv"
    pd.DataFrame({"gene_id":["a","b","c"],"ctrlA":[0,2,4],"ctrlB":[0,4,6],"Aux":[999,999,999]}).to_csv(p,sep="\t",index=False)
    d=baseline_expression(p,["ctrlA","ctrlB"])
    assert d.gene_id.tolist()==["b","c"]
    assert d.untreated_mean_normalized_count.tolist()==[3,5]
    assert "Aux" not in d


def test_value_tertiles_keep_ties_together():
    d=pd.DataFrame({"gene_id":list("abcdefghi"),"untreated_mean_normalized_count":[1,1,1,2,2,2,3,3,3]})
    d["log10_untreated_expression_plus_1"]=np.log10(d.untreated_mean_normalized_count+1)
    out,cut=assign_tertiles(d)
    assert out.groupby("untreated_mean_normalized_count").expression_tertile.nunique().max()==1
    assert set(out.expression_tertile)=={"Low","Middle","High"}


def test_missing_untreated_column_fails(tmp_path):
    p=tmp_path/"norm.tsv"; pd.DataFrame({"gene_id":["a"],"ctrlA":[1]}).to_csv(p,sep="\t",index=False)
    try: baseline_expression(p,["ctrlA","ctrlB"])
    except ValueError as e: assert "ctrlB" in str(e)
    else: raise AssertionError("missing untreated replicate was accepted")


def test_eligible_gene_union(tmp_path):
    a=tmp_path/"a.tsv"; b=tmp_path/"b.tsv"
    pd.DataFrame({"gene_id":["a","b"]}).to_csv(a,sep="\t",index=False)
    pd.DataFrame({"gene_id":["b","c"]}).to_csv(b,sep="\t",index=False)
    assert eligible_gene_ids([a,b])=={"a","b","c"}
