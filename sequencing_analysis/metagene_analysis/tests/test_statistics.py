import pandas as pd
from spen_metagene.statistics import stratified_tests


def test_stratified_test_is_deterministic():
    metrics=pd.DataFrame({"gene_id":[f"g{i}" for i in range(20)],"first_5kb":[0]*10+[2]*10})
    classes=pd.DataFrame({"gene_id":metrics.gene_id,"class_fdr05":["no_change"]*10+["significant_elevation"]*10,"baseMean":range(20),"gene_length":[1000]*20,"intron_count":[2]*20})
    a=stratified_tests(metrics,classes,"class_fdr05",50,7)
    b=stratified_tests(metrics,classes,"class_fdr05",50,7)
    assert a.equals(b)
    assert a.loc[0,"median_difference"]==2
