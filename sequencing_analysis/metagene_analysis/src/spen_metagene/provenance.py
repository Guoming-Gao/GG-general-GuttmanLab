from __future__ import annotations
import json, platform, subprocess, sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


def write_versions(output):
    packages={}
    for p in ["spen-metagene","numpy","pandas","matplotlib","seaborn","pysam","pyBigWig","snakemake","pytest"]:
        try: packages[p]=version(p)
        except PackageNotFoundError: packages[p]=None
    def cmd(x):
        try: return subprocess.check_output(x,text=True,stderr=subprocess.STDOUT).strip()
        except Exception as e: return repr(e)
    data={"python":sys.version,"platform":platform.platform(),"packages":packages,
          "samtools":cmd(["samtools","--version"]),"bedtools":cmd(["bedtools","--version"]),
          "java":cmd(["java","-version"]),
          "R_packages":cmd(["/opt/miniconda3/envs/bioinfo/bin/Rscript","-e",'cat(R.version.string,"\\n"); for(x in c("DESeq2","apeglm","edgeR")) cat(x,as.character(packageVersion(x)),"\\n")'])}
    Path(output).parent.mkdir(parents=True,exist_ok=True); Path(output).write_text(json.dumps(data,indent=2))
