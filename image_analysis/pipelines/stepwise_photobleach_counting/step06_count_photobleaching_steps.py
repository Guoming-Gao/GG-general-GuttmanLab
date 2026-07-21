"""Step 6: run quickPBSA SIC plus Bayesian photobleaching-step counting."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from pbsa_shared import atomic_csv, condition_metadata, load_config, matplotlib_setup, output_action, output_dirs, stage_timer
from quickpbsa_compat import result_summary, run_quickpbsa
from step01_inspect_inputs import load_manifest
from step04_extract_background_corrected_traces import trace_paths
from step05_quickpbsa_threshold_pilot import selected_thresholds_path


def count_paths(output_dir:str|Path,fov:str)->dict[str,Path]:
    root=output_dirs(output_dir)["step_counts"]/fov
    return {"canonical":root/f"{fov}__photobleaching_step_counts.csv","plot":root/f"{fov}__photobleaching_step_counting_QC.png","raw":root/"quickpbsa_native"}


def configured_threshold(cfg:dict[str,Any],profile:str)->float:
    value=cfg["quickpbsa"].get("thresholds",{}).get(profile)
    if value is None:
        path=selected_thresholds_path(cfg)
        if path.exists():
            table=pd.read_csv(path); matched=table[table.acquisition_profile==profile]
            if len(matched): value=matched.iloc[0].selected_threshold
    if value is None: raise RuntimeError(f"Missing quickPBSA threshold for acquisition profile '{profile}'. Run the pooled stage-05 pilot or add an explicit config override.")
    return float(value)


def save_fit_plot(input_table:pd.DataFrame,result:pd.DataFrame,canonical:pd.DataFrame,path:Path,dataset:dict[str,Any],maximum:int)->None:
    matplotlib_setup(dataset["output_dir"]); import matplotlib.pyplot as plt
    frame_cols=[column for column in input_table.columns if str(column).isdigit()]; full=result[result.type=="fluors_full"].reset_index(drop=True)
    result_frame_cols=[column for column in full.columns if str(column).isdigit()]
    fig,axes=plt.subplots(2,1,figsize=(11,7),squeeze=False); axes=axes[:,0]
    for roi_class,ax in zip(["point_punctum","extended_condensate"],axes):
        indices=np.flatnonzero(canonical.roi_class.to_numpy()==roi_class)[:maximum]
        for idx in indices:
            ax.plot(input_table.loc[idx,frame_cols].to_numpy(float),alpha=.28,color="gray")
            if idx<len(full): ax.plot(full.loc[idx,result_frame_cols].to_numpy(float),alpha=.8,lw=1)
        ax.set_title(roi_class); ax.set_xlabel("Frame"); ax.set_ylabel("Integrated intensity")
    fig.tight_layout(); path.parent.mkdir(parents=True,exist_ok=True); fig.savefig(path,dpi=200); plt.close(fig)


def process_file(row:dict[str,Any],cfg:dict[str,Any],dataset:dict[str,Any],*,resume:bool,force:bool)->dict[str,Any]:
    paths=count_paths(dataset["output_dir"],row["fov"]); required=[paths["canonical"],paths["plot"]]
    if output_action(required,resume=resume,force=force)=="skip":
        table=pd.read_csv(paths["canonical"]); accepted=table[table.pbsa_flag==1]
        return {"filename":row["filename"],"fov":row["fov"],"status":"complete","rois":len(table),"accepted":len(accepted),
                "accepted_fraction":len(accepted)/len(table) if len(table) else 0,
                "point_puncta_accepted":int((accepted.roi_class=="point_punctum").sum()),
                "extended_condensates_accepted":int((accepted.roi_class=="extended_condensate").sum()),
                "outside_validated_range":int((accepted.validated_range_status=="outside_validated_range").sum()),"plot":str(paths["plot"])}
    infile=trace_paths(dataset["output_dir"],row["fov"])["integrated_difference"]
    if not infile.exists(): raise FileNotFoundError(f"Run trace extraction first: {infile}")
    threshold=configured_threshold(cfg,row["acquisition_profile"]); input_table=pd.read_csv(infile)
    result=run_quickpbsa(infile,paths["raw"],threshold,cfg["quickpbsa"]); summary=result_summary(result)
    metadata=[column for column in ["roi_id","roi_class","area_px","equivalent_diameter_px","centroid_y_px","centroid_x_px","seed_count"] if column in input_table]
    canonical=input_table[metadata].reset_index(drop=True).copy(); canonical["pbsa_flag"]=pd.to_numeric(summary.get("flag"),errors="coerce")
    canonical["photobleaching_step_count"]=pd.to_numeric(summary.get("photobleaching_step_count"),errors="coerce"); canonical["quickpbsa_threshold"]=threshold
    canonical["dataset"]=dataset["name"]; canonical["acquisition_profile"]=row["acquisition_profile"]
    for key,value in condition_metadata(row["filename"]).items(): canonical[key]=value
    canonical["accepted_count"]=canonical.pbsa_flag==1
    maximum=int(cfg["quickpbsa"]["maximum_validated_count"])
    canonical["validated_range_status"] = np.select(
        [
            canonical.photobleaching_step_count.isna(),
            canonical.photobleaching_step_count > maximum,
        ],
        ["not_counted", "outside_validated_range"],
        default="within_validated_range",
    )
    canonical.insert(0,"filename",row["filename"]); canonical.insert(1,"fov",row["fov"]); atomic_csv(canonical,paths["canonical"])
    save_fit_plot(input_table,result,canonical,paths["plot"],dataset,min(12,int(cfg["qc"]["max_trace_panels_per_class"])))
    accepted=canonical[canonical.accepted_count]
    return {"filename":row["filename"],"fov":row["fov"],"status":"complete","rois":len(canonical),"accepted":len(accepted),"accepted_fraction":len(accepted)/len(canonical) if len(canonical) else 0,
            "point_puncta_accepted":int(((accepted.roi_class=="point_punctum")).sum()),"extended_condensates_accepted":int(((accepted.roi_class=="extended_condensate")).sum()),
            "outside_validated_range":int((accepted.validated_range_status=="outside_validated_range").sum()),"plot":str(paths["plot"])}


def save_stage_report(summary:pd.DataFrame,dataset:dict[str,Any])->None:
    matplotlib_setup(dataset["output_dir"]); import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    root=output_dirs(dataset["output_dir"])["step_counts"]; report=root/"photobleaching_step_counting_QC_report.pdf"
    with PdfPages(report) as pdf:
        fig,ax=plt.subplots(figsize=(11,8.5)); ax.axis("off"); ax.text(.02,.98,f"Photobleaching step-counting QC — {dataset['name']}\n\n"+summary.to_string(index=False),va="top",family="monospace",fontsize=7); pdf.savefig(fig); plt.close(fig)
        for path in summary["plot"]:
            fig,ax=plt.subplots(figsize=(11,8.5)); ax.imshow(plt.imread(path)); ax.axis("off"); pdf.savefig(fig); plt.close(fig)


def step_counting_stage(cfg:dict[str,Any],dataset:dict[str,Any],manifest:pd.DataFrame,*,resume:bool=False,force:bool=False,max_files:int|None=None)->pd.DataFrame:
    rows=manifest[manifest.status=="accepted"].to_dict("records")[:max_files]; root=output_dirs(dataset["output_dir"])["step_counts"]; root.mkdir(parents=True,exist_ok=True); records=[]
    missing=[]
    for profile in sorted({row["acquisition_profile"] for row in rows}):
        try: configured_threshold(cfg,profile)
        except RuntimeError: missing.append(profile)
    if missing: raise RuntimeError("Missing quickPBSA thresholds for: "+", ".join(missing))
    with stage_timer(dataset["output_dir"],"06_photobleaching_step_counts",{"files":len(rows)}):
        for index,row in enumerate(rows,1): print(f"[photobleaching step counting {index}/{len(rows)}] {row['filename']}",flush=True); records.append(process_file(row,cfg,dataset,resume=resume,force=force))
        summary=pd.DataFrame.from_records(records); atomic_csv(summary,root/"photobleaching_step_counting_QC_summary.csv"); save_stage_report(summary,dataset)
    return summary


def main(argv:list[str]|None=None)->None:
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--config",required=True); parser.add_argument("--max-files",type=int); policy=parser.add_mutually_exclusive_group(); policy.add_argument("--resume",action="store_true"); policy.add_argument("--force",action="store_true")
    args=parser.parse_args(argv); cfg=load_config(args.config)
    for dataset in cfg["datasets"]: step_counting_stage(cfg,dataset,load_manifest(dataset),resume=args.resume,force=args.force,max_files=args.max_files)


if __name__=="__main__": main()
