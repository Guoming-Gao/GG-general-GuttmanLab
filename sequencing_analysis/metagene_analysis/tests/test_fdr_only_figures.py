from pathlib import Path
import numpy as np
import pandas as pd

from spen_metagene.plotting import paired_volcano
from spen_metagene.suite import _scaled_width_kb


def test_paired_volcano_uses_padj_and_fixed_classes(tmp_path: Path):
    d=pd.DataFrame({
        "gene_id":["a","b"], "gene_name":["A","B"],
        "log2FoldChange_mle":[1.0,-1.0], "log2FoldChange_apeglm":[0.4,-0.2],
        "padj":[0.01,0.2], "pvalue":[1e-20,1e-30],
        "class_fdr05":["significant_elevation","nonsignificant_decrease"],
    })
    source=tmp_path/"classes.tsv"; output=tmp_path/"paired.png"
    d.to_csv(source,sep="\t",index=False)
    paired_volcano(source,output,"class_fdr05","test",np.log2(1.2),.05)
    plotted=pd.read_csv(tmp_path/"paired.tsv",sep="\t")
    assert np.allclose(plotted["plotted_y"],-np.log10(d.padj))
    assert plotted["class_fdr05"].tolist()==d.class_fdr05.tolist()
    assert plotted["class_color_source"].eq("MLE LFC + BH padj").all()


def test_scaled_component_width_normalization():
    d=pd.DataFrame({"bin":[0,20,119,120],"gene_length":[20000]*4})
    assert np.allclose(_scaled_width_kb(d),[.1,.2,.2,.1])
