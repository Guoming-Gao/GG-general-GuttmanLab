from pathlib import Path
import pandas as pd


def write_audit(root, output):
    root=Path(root)
    rows=[
      ("six paired MLE/apeglm FDR05 volcanoes","complete",root/"figures/volcano","Both panels use BH padj and identical MLE/FDR05 class colors"),
      ("FDR10 sensitivity tables","complete",root/"gene_classes","Machine-readable only; not displayed"),
      ("early and late fragment count matrices","reused",root/"intermediate/count_matrices","Validated existing results; BAM recount avoided"),
      ("CLAP positional and normalized matrices","reused",root/"clap_matrices","Validated existing results; BAM recount avoided"),
      ("five MLE-based gene classes","complete",root/"gene_classes","Corrected after prior apeglm classification was invalidated"),
      ("1h/2h early agreement","complete",root/"gene_classes","MLE direction; 24h excluded"),
      ("scaled continuous/component/shape metagenes","complete",root/"metagene","100 seeded bootstrap iterations; exploratory intervals"),
      ("categorical violin-plus-box features","complete",root/"metagene","No feature categories are connected by lines; component density is width normalized"),
      ("24h autosome and chrX panels","complete",root/"metagene","Autosomes headline; chrX control"),
      ("archived CA and EBS-10 sensitivity","reused",root/"historical","Prior caller output reused"),
      ("EBS-1000 stability rerun","trimmed",root/"historical","Disabled: redundant and approximately one hour per sample"),
      ("full-file count/reference checksums","complete",root/"inputs_manifest/input_manifest.tsv","BAMs use identity metadata plus first-MiB digest"),
      ("full and trimmed self-contained reports","complete",root/"report","Trimmed report contains the six paired volcano, scaled-profile, and feature figures"),
    ]
    out=[]
    for item,status,path,note in rows:
        exists=path.exists()
        if status=="complete" and not exists: status="missing"
        out.append({"planned_item":item,"status":status,"artifact":str(path),"artifact_exists":exists,"note":note})
    p=Path(output); p.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(out).to_csv(p,sep="\t",index=False)
