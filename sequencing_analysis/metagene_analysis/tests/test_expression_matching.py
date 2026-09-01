import pandas as pd
from spen_metagene.matched import select_no_change, balance_table

def _base():
    return pd.DataFrame({
        "gene_id":["n1","n2","n3","n4","d1","d2","e1","e2"],
        "class":["no_change"]*4+["nonsignificant_decrease"]*2+["nonsignificant_elevation"]*2,
        "value":[.11,.12,.25,.35,.13,.14,.16,.26],
    })

def test_matching_is_deterministic_unique_and_bin_limited():
    a,bins,target,controls=select_no_change(_base(),.1,17)
    c,_,_,_=select_no_change(_base(),.1,17)
    assert a.gene_id.tolist()==c.gene_id.tolist()
    assert a.gene_id.is_unique
    assert set(a.gene_id)<=set(controls.gene_id)
    assert (bins.n_no_change_selected<=bins.n_no_change_available).all()
    assert (bins.n_no_change_selected<=bins.n_directional).all()
    assert len(target)==4

def test_shortages_and_balance_are_reported():
    selected,bins,target,controls=select_no_change(_base(),.1,17)
    assert bins.shortage.sum()>0
    out=balance_table(target,controls,selected)
    assert set(out.group)=={"directional_target","no_change_unmatched","no_change_matched"}
    assert out.loc[out.group.eq("no_change_matched"),"ks_vs_directional"].notna().all()
