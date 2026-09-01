import numpy as np
import pandas as pd

from spen_metagene.classes import classify_frame, early_agreement


def test_five_bins_and_not_tested():
    d=pd.DataFrame({"log2FoldChange_mle":[-.5,-.5,0,.5,.5,.2,np.nan],"padj":[.01,.2,.8,.2,.01,np.nan,.01]})
    assert classify_frame(d,.05,.263).tolist()==[
        "significant_decrease","nonsignificant_decrease","no_change",
        "nonsignificant_elevation","significant_elevation","not_tested","not_tested"]


def test_boundary_is_no_change():
    d=pd.DataFrame({"log2FoldChange_mle":[-.263,.263],"padj":[.001,.001]})
    assert classify_frame(d,.05,.263).tolist()==["no_change","no_change"]


def test_early_directional_agreement_does_not_require_two_significant_times():
    one=pd.DataFrame({"gene_id":["g"],"gene_name":["G"],"chrom":["chr1"],"log2FoldChange_mle":[.5],"padj":[.01],"class_fdr05":["significant_elevation"]})
    two=pd.DataFrame({"gene_id":["g"],"gene_name":["G"],"chrom":["chr1"],"log2FoldChange_mle":[.4],"padj":[.2],"class_fdr05":["nonsignificant_elevation"]})
    out=early_agreement(one,two,.05,.263)
    assert out.loc[0,"agreement"]=="concordant_elevation"
    assert out.loc[0,"significance_support"]=="one_timepoint"


def test_mle_not_apeglm_defines_primary_bin():
    d=pd.DataFrame({"log2FoldChange_mle":[.8],"log2FoldChange_apeglm":[.02],"padj":[.2]})
    assert classify_frame(d,.05,.263).item()=="nonsignificant_elevation"
