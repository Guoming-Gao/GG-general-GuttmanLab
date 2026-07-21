"""Step 4: extract point-ring and mask-band locally background-corrected traces."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import distance_transform_edt
from scipy.spatial.distance import cdist

from pbsa_shared import atomic_csv, load_config, matplotlib_setup, output_action, output_dirs, stage_timer
from step01_inspect_inputs import load_manifest
from step02_drift_correction import drift_correction_paths
from step03_detect_and_grow_rois import roi_paths


TRACE_KINDS = ("peak_mean", "background_mean", "mean_difference", "integrated_difference")


def trace_paths(output_dir: str | Path, fov: str) -> dict[str, Path]:
    root = output_dirs(output_dir)["traces"] / fov
    return {
        "peak_mean": root / f"{fov}__peak_mean.csv",
        "background_mean": root / f"{fov}__background_mean.csv",
        "mean_difference": root / f"{fov}__mean_difference.csv",
        "integrated_difference": root / f"{fov}__integrated_difference.csv",
        "roi_qc": root / f"{fov}__trace_extraction_QC.csv",
        "signal_mask": root / f"{fov}__signal_ROIs.tif",
        "background_mask": root / f"{fov}__background_ROIs.tif",
        "trace_plot": root / f"{fov}__background_corrected_trace_QC.png",
    }


def _point_indices(shape: tuple[int, int], y: float, x: float, radii: dict[str, float], foreground: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    outer = float(radii["point_background_outer_radius_px"]); y0=max(0,int(np.floor(y-outer))); y1=min(shape[0],int(np.ceil(y+outer))+1)
    x0=max(0,int(np.floor(x-outer))); x1=min(shape[1],int(np.ceil(x+outer))+1); yy,xx=np.ogrid[y0:y1,x0:x1]; distance=np.hypot(yy-y,xx-x)
    signal=(distance<=float(radii["point_radius_px"])); background=(distance>=float(radii["point_background_inner_radius_px"]))&(distance<=outer)&(~foreground[y0:y1,x0:x1])
    flat=np.arange(shape[0]*shape[1]).reshape(shape)
    return flat[y0:y1,x0:x1][signal],flat[y0:y1,x0:x1][background]


def build_roi_indices(labels: np.ndarray, seeds: pd.DataFrame, regions: pd.DataFrame, cfg: dict[str, Any]) -> tuple[list[dict[str, Any]], np.ndarray, np.ndarray]:
    radii=cfg["trace_extraction"]; shape=labels.shape; foreground=labels>0; flat_labels=labels.ravel(); entries=[]
    seed_xy=seeds[["y_px","x_px"]].to_numpy(float) if not seeds.empty else np.empty((0,2)); distances=cdist(seed_xy,seed_xy) if len(seed_xy) else np.empty((0,0))
    if distances.size: np.fill_diagonal(distances,np.inf)
    nearest={int(seed_id):float(distances[i].min()) for i,seed_id in enumerate(seeds.seed_id)} if len(seeds) else {}
    signal_display=np.zeros(shape,np.uint16); background_display=np.zeros(shape,np.uint16)
    for region in regions.itertuples():
        status="accepted"; reason=""; signal_idx=np.array([],int); background_idx=np.array([],int)
        group=seeds[seeds.roi_id==region.roi_id]
        if region.roi_class=="point_punctum":
            seed=group.iloc[0] if len(group) else None
            if seed is None: status,reason="rejected","missing_seed"
            elif nearest.get(int(seed.seed_id),np.inf)<float(radii["point_minimum_neighbor_distance_px"]): status,reason="rejected","close_neighbor"
            else: signal_idx,background_idx=_point_indices(shape,float(seed.y_px),float(seed.x_px),radii,foreground)
        else:
            mask=labels==int(region.roi_id); signal_idx=np.flatnonzero(mask); distance=distance_transform_edt(~mask)
            gap=float(radii["mask_background_gap_px"]); width=float(radii["mask_background_width_px"])
            background_idx=np.flatnonzero((distance>gap)&(distance<=gap+width)&(~foreground))
        if status=="accepted" and not len(signal_idx): status,reason="rejected","empty_signal_roi"
        if status=="accepted" and not len(background_idx): status,reason="rejected","empty_background_roi"
        if status=="accepted":
            signal_display.ravel()[signal_idx]=int(region.roi_id); background_display.ravel()[background_idx]=int(region.roi_id)
        entries.append({"roi_id":int(region.roi_id),"roi_class":region.roi_class,"trace_status":status,"trace_rejection_reason":reason,
                        "signal_pixels":len(signal_idx),"background_pixels":len(background_idx),"signal_indices":signal_idx,"background_indices":background_idx})
    return entries,signal_display,background_display


def extract_traces(stack_path: str | Path, entries: list[dict[str, Any]], max_frames: int | None = None) -> dict[str, np.ndarray]:
    accepted=[entry for entry in entries if entry["trace_status"]=="accepted"]
    with tifffile.TiffFile(stack_path) as tif:
        frames=min(len(tif.pages),max_frames) if max_frames else len(tif.pages); n=len(accepted)
        output={name:np.empty((n,frames),np.float32) for name in TRACE_KINDS}
        for frame in range(frames):
            image=tif.pages[frame].asarray().ravel().astype(np.float32,copy=False)
            for i,entry in enumerate(accepted):
                signal=image[entry["signal_indices"]]; background=image[entry["background_indices"]]
                peak_mean=float(np.mean(signal)); bg_mean=float(np.mean(background))
                output["peak_mean"][i,frame]=peak_mean; output["background_mean"][i,frame]=bg_mean
                output["mean_difference"][i,frame]=peak_mean-bg_mean
                output["integrated_difference"][i,frame]=float(np.sum(signal,dtype=np.float64)-len(signal)*bg_mean)
    return output


def trace_base_table(entries: list[dict[str, Any]], regions: pd.DataFrame) -> pd.DataFrame:
    accepted=pd.DataFrame([{key:value for key,value in entry.items() if not key.endswith("_indices")} for entry in entries if entry["trace_status"]=="accepted"])
    columns=[column for column in ["roi_id","roi_class","area_px","equivalent_diameter_px","centroid_y_px","centroid_x_px","seed_count"] if column in regions.columns]
    return accepted.merge(regions[columns],on=["roi_id","roi_class"],how="left") if not accepted.empty else accepted


def write_trace_csv(base: pd.DataFrame, values: np.ndarray, path: Path) -> None:
    frames=pd.DataFrame(values,columns=[str(i) for i in range(values.shape[1])]); atomic_csv(pd.concat([base.reset_index(drop=True),frames],axis=1),path)


def trace_metrics(base: pd.DataFrame, traces: dict[str,np.ndarray]) -> pd.DataFrame:
    records=[]
    for i,row in base.iterrows():
        peak=traces["peak_mean"][i]; bg=traces["background_mean"][i]; diff=traces["integrated_difference"][i]; n=len(diff); late=max(5,n//10)
        x=np.arange(n,dtype=float); bg_slope=float(np.polyfit(x,bg,1)[0]) if n>1 else np.nan; late_slope=float(np.polyfit(x[-late:],diff[-late:],1)[0]) if late>1 else np.nan
        noise=float(1.4826*np.median(np.abs(np.diff(diff[-late:])-np.median(np.diff(diff[-late:]))))) if late>2 else np.nan
        records.append({**row.to_dict(),"initial_peak_mean":float(np.mean(peak[:min(5,n)])),"initial_background_mean":float(np.mean(bg[:min(5,n)])),
                        "background_decay_slope_per_frame":bg_slope,"late_integrated_baseline":float(np.mean(diff[-late:])),
                        "late_integrated_slope_per_frame":late_slope,"late_difference_noise":noise})
    return pd.DataFrame.from_records(records)


def save_trace_plot(qc:pd.DataFrame,traces:dict[str,np.ndarray],path:Path,dataset:dict[str,Any],maximum:int)->None:
    matplotlib_setup(dataset["output_dir"]); import matplotlib.pyplot as plt
    path.parent.mkdir(parents=True,exist_ok=True); classes=[value for value in ["point_punctum","extended_condensate"] if value in set(qc.roi_class)]
    rows=max(1,len(classes)); fig,axes=plt.subplots(rows,3,figsize=(12,3.5*rows),squeeze=False)
    for r,roi_class in enumerate(classes):
        indices=np.flatnonzero(qc.roi_class.to_numpy()==roi_class)[:maximum]
        for idx in indices: axes[r,0].plot(traces["peak_mean"][idx],alpha=.35); axes[r,1].plot(traces["background_mean"][idx],alpha=.35); axes[r,2].plot(traces["integrated_difference"][idx],alpha=.35)
        axes[r,0].set_ylabel(roi_class); axes[r,0].set_title("Peak mean"); axes[r,1].set_title("Local background mean"); axes[r,2].set_title("Integrated difference")
    for ax in axes.ravel(): ax.set_xlabel("Frame")
    fig.tight_layout(); fig.savefig(path,dpi=200); plt.close(fig)


def process_file(row:dict[str,Any],cfg:dict[str,Any],dataset:dict[str,Any],*,resume:bool,force:bool,max_frames:int|None)->dict[str,Any]:
    paths=trace_paths(dataset["output_dir"],row["fov"])
    if output_action(list(paths.values()),resume=resume,force=force)=="skip":
        qc=pd.read_csv(paths["roi_qc"]); return {"filename":row["filename"],"fov":row["fov"],"status":"skipped","accepted_rois":len(qc),"point_puncta":int((qc.roi_class=="point_punctum").sum()),"extended_condensates":int((qc.roi_class=="extended_condensate").sum()),"trace_plot":str(paths["trace_plot"])}
    corrected=drift_correction_paths(dataset["output_dir"],row["fov"])["stack"]; rpaths=roi_paths(dataset["output_dir"],row["fov"])
    if not corrected.exists() or not rpaths["mask"].exists(): raise FileNotFoundError("Run drift correction and ROI detection first")
    labels=tifffile.imread(rpaths["mask"]); seeds=pd.read_csv(rpaths["seeds"]); regions=pd.read_csv(rpaths["regions"])
    entries,signal_mask,bg_mask=build_roi_indices(labels,seeds,regions,cfg); traces=extract_traces(corrected,entries,max_frames=max_frames); base=trace_base_table(entries,regions)
    for kind in TRACE_KINDS: write_trace_csv(base,traces[kind],paths[kind])
    qc=trace_metrics(base,traces); rejected=[{key:value for key,value in entry.items() if not key.endswith("_indices")} for entry in entries if entry["trace_status"]!="accepted"]
    if rejected:
        rejected_df=pd.DataFrame(rejected); rejected_df["accepted_for_trace_extraction"]=False
    else:
        rejected_df=pd.DataFrame()
    qc["accepted_for_trace_extraction"]=True; atomic_csv(pd.concat([qc,rejected_df],ignore_index=True,sort=False),paths["roi_qc"])
    paths["signal_mask"].parent.mkdir(parents=True,exist_ok=True); tifffile.imwrite(paths["signal_mask"],signal_mask,metadata={"axes":"YX"}); tifffile.imwrite(paths["background_mask"],bg_mask,metadata={"axes":"YX"})
    save_trace_plot(qc,traces,paths["trace_plot"],dataset,int(cfg["qc"]["max_trace_panels_per_class"]))
    return {"filename":row["filename"],"fov":row["fov"],"status":"complete","accepted_rois":len(qc),"rejected_rois":len(rejected),
            "point_puncta":int((qc.roi_class=="point_punctum").sum()),"extended_condensates":int((qc.roi_class=="extended_condensate").sum()),
            "median_late_baseline":float(qc.late_integrated_baseline.median()) if len(qc) else np.nan,"median_late_slope":float(qc.late_integrated_slope_per_frame.median()) if len(qc) else np.nan,"trace_plot":str(paths["trace_plot"])}


def save_stage_report(summary:pd.DataFrame,dataset:dict[str,Any])->None:
    matplotlib_setup(dataset["output_dir"]); import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    root=output_dirs(dataset["output_dir"])["traces"]; report=root/"background_corrected_trace_QC_report.pdf"
    with PdfPages(report) as pdf:
        fig,ax=plt.subplots(figsize=(11,8.5)); ax.axis("off"); ax.text(.02,.98,f"Background-corrected trace QC — {dataset['name']}\n\n"+summary.to_string(index=False),va="top",family="monospace",fontsize=7); pdf.savefig(fig); plt.close(fig)
        for path in summary.trace_plot:
            fig,ax=plt.subplots(figsize=(11,8.5)); ax.imshow(plt.imread(path)); ax.axis("off"); pdf.savefig(fig); plt.close(fig)


def trace_extraction_stage(cfg:dict[str,Any],dataset:dict[str,Any],manifest:pd.DataFrame,*,resume:bool=False,force:bool=False,max_files:int|None=None,max_frames:int|None=None)->pd.DataFrame:
    rows=manifest[manifest.status=="accepted"].to_dict("records")[:max_files]; root=output_dirs(dataset["output_dir"])["traces"]; root.mkdir(parents=True,exist_ok=True); records=[]
    with stage_timer(dataset["output_dir"],"04_background_corrected_traces",{"files":len(rows),"max_frames":max_frames}):
        for index,row in enumerate(rows,1): print(f"[background-corrected traces {index}/{len(rows)}] {row['filename']}",flush=True); records.append(process_file(row,cfg,dataset,resume=resume,force=force,max_frames=max_frames))
        summary=pd.DataFrame.from_records(records); atomic_csv(summary,root/"background_corrected_trace_QC_summary.csv"); save_stage_report(summary,dataset)
    return summary


def main(argv:list[str]|None=None)->None:
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--config",required=True); parser.add_argument("--max-files",type=int); parser.add_argument("--max-frames",type=int); policy=parser.add_mutually_exclusive_group(); policy.add_argument("--resume",action="store_true"); policy.add_argument("--force",action="store_true")
    args=parser.parse_args(argv); cfg=load_config(args.config)
    for dataset in cfg["datasets"]: trace_extraction_stage(cfg,dataset,load_manifest(dataset),resume=args.resume,force=args.force,max_files=args.max_files,max_frames=args.max_frames)


if __name__=="__main__": main()
