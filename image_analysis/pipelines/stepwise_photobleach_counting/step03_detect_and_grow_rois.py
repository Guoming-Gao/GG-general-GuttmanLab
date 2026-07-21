"""Step 3: detect Spotiflow seeds and grow local point/condensate ROIs."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import gaussian_filter, label
from skimage.measure import regionprops
from skimage.segmentation import find_boundaries

from pbsa_shared import atomic_csv, load_config, matplotlib_setup, output_action, output_dirs, stage_timer
from step01_inspect_inputs import load_manifest
from step02_drift_correction import drift_correction_paths


def roi_paths(output_dir: str | Path, fov: str) -> dict[str, Path]:
    root = output_dirs(output_dir)["roi_detection"]
    return {
        "image": root / "detection_images" / f"{fov}__first_frames_mean.tif",
        "seeds": root / "tables" / f"{fov}__spotiflow_seeds.csv",
        "regions": root / "tables" / f"{fov}__roi_properties.csv",
        "mask": root / "masks" / f"{fov}__roi_labels.tif",
        "overlay": root / "overlays" / f"{fov}__roi_detection_and_growth.png",
    }


def first_frames_mean(path: str | Path, count: int) -> np.ndarray:
    with tifffile.TiffFile(path) as tif:
        n = min(count, len(tif.pages))
        return np.mean([tif.pages[i].asarray().astype(np.float32) for i in range(n)], axis=0, dtype=np.float32)


def load_spotiflow(cfg: dict[str, Any]) -> SimpleNamespace:
    # Spotiflow imports plotting modules; keep their cache inside the result tree
    # in environments where the user-level Matplotlib directory is read-only.
    matplotlib_setup(cfg["datasets"][0]["output_dir"])
    from spotiflow.model import Spotiflow
    settings = cfg["spotiflow"]
    model = Spotiflow.from_pretrained(settings["pretrained_model"], cache_dir=settings["model_cache_dir"], map_location="cpu", verbose=False)
    return SimpleNamespace(model=model)


def detect_seeds(image: np.ndarray, cfg: dict[str, Any], context: SimpleNamespace) -> pd.DataFrame:
    settings = cfg["spotiflow"]
    points, details = context.model.predict(
        image, prob_thresh=float(settings["probability_threshold"]), min_distance=int(settings["min_distance"]),
        exclude_border=False, normalizer=settings["normalizer"], device=settings["device"], verbose=False,
    )
    def detail(name: str, index: int) -> float:
        value = getattr(details, name, None)
        if value is None or len(value) <= index: return float("nan")
        current = np.asarray(value[index]).squeeze()
        return float(current) if current.size == 1 else float("nan")
    rows = []
    for index, point in enumerate(np.asarray(points)):
        rows.append({"seed_id": index + 1, "y_px": float(point[0]), "x_px": float(point[1]),
                     "spotiflow_probability": detail("prob", index), "spotiflow_intensity": detail("intens", index)})
    return pd.DataFrame.from_records(rows, columns=["seed_id","y_px","x_px","spotiflow_probability","spotiflow_intensity"])


def robust_background(values: np.ndarray) -> tuple[float, float]:
    data = values[np.isfinite(values)].astype(np.float32)
    if not data.size: return float("nan"), float("nan")
    for _ in range(3):
        median = float(np.median(data)); sigma = float(1.4826 * np.median(np.abs(data - median)))
        if sigma <= 0: break
        clipped = data[data <= median + 3 * sigma]
        if len(clipped) == len(data) or not len(clipped): break
        data = clipped
    median = float(np.median(data)); sigma = float(1.4826 * np.median(np.abs(data - median)))
    return median, max(sigma, 1e-6)


def grow_seed_regions(image: np.ndarray, seeds: pd.DataFrame, cfg: dict[str, Any]) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    settings = cfg["roi_growth"]; smooth = gaussian_filter(image.astype(np.float32), float(settings["smooth_sigma_px"]))
    radius = int(settings["maximum_radius_px"]); outer_width = int(settings["outer_background_width_px"])
    union = np.zeros(image.shape, bool); seed_records = []
    h, w = image.shape
    for row in seeds.itertuples():
        yi, xi = int(round(row.y_px)), int(round(row.x_px)); y0=max(0,yi-radius); y1=min(h,yi+radius+1); x0=max(0,xi-radius); x1=min(w,xi+radius+1)
        crop=smooth[y0:y1,x0:x1]; yy,xx=np.ogrid[y0:y1,x0:x1]; dist=np.hypot(yy-yi,xx-xi)
        outer=(dist >= max(1,radius-outer_width)) & (dist <= radius)
        background, noise=robust_background(crop[outer]); threshold=background + float(settings["foreground_sigma"])*noise
        binary=crop > threshold; sy,sx=yi-y0,xi-x0
        if sy<0 or sx<0 or sy>=binary.shape[0] or sx>=binary.shape[1] or not binary[sy,sx]:
            seed_records.append({"seed_id":row.seed_id,"growth_status":"seed_below_local_threshold","local_background":background,"local_noise":noise,"growth_threshold":threshold,"growth_limit_reached":False}); continue
        labels,_=label(binary); component=labels==labels[sy,sx]
        if int(component.sum()) < int(settings["minimum_area_px"]):
            seed_records.append({"seed_id":row.seed_id,"growth_status":"region_too_small","local_background":background,"local_noise":noise,"growth_threshold":threshold,"growth_limit_reached":False}); continue
        touches=bool(component[0].any() or component[-1].any() or component[:,0].any() or component[:,-1].any())
        union[y0:y1,x0:x1] |= component
        seed_records.append({"seed_id":row.seed_id,"growth_status":"grown_with_warning" if touches else "grown","local_background":background,"local_noise":noise,"growth_threshold":threshold,"growth_limit_reached":touches})
    roi_labels,n=label(union)
    seed_qc_columns = [
        "seed_id", "growth_status", "local_background", "local_noise",
        "growth_threshold", "growth_limit_reached",
    ]
    seed_qc=pd.DataFrame.from_records(seed_records, columns=seed_qc_columns)
    seeds_out=seeds.merge(seed_qc,on="seed_id",how="left"); assignments=[]
    for row in seeds_out.itertuples():
        yi,xi=int(round(row.y_px)),int(round(row.x_px)); roi_id=int(roi_labels[yi,xi]) if 0<=yi<h and 0<=xi<w else 0
        assignments.append(roi_id)
    seeds_out["roi_id"]=assignments
    seed_counts=seeds_out[seeds_out.roi_id>0].groupby("roi_id").seed_id.count().to_dict()
    limit_flags=seeds_out[seeds_out.roi_id>0].groupby("roi_id").growth_limit_reached.any().to_dict()
    records=[]
    for prop in regionprops(roi_labels):
        seed_count=int(seed_counts.get(prop.label,0)); eq=float(prop.equivalent_diameter_area)
        roi_class="point_punctum" if seed_count==1 and eq<=float(settings["point_max_equivalent_diameter_px"]) else "extended_condensate"
        records.append({"roi_id":int(prop.label),"roi_class":roi_class,"area_px":int(prop.area),"equivalent_diameter_px":eq,
                        "centroid_y_px":float(prop.centroid[0]),"centroid_x_px":float(prop.centroid[1]),"seed_count":seed_count,
                        "growth_limit_reached":bool(limit_flags.get(prop.label,False)),"bbox_min_y":int(prop.bbox[0]),"bbox_min_x":int(prop.bbox[1]),
                        "bbox_max_y":int(prop.bbox[2]),"bbox_max_x":int(prop.bbox[3])})
    region_columns = [
        "roi_id", "roi_class", "area_px", "equivalent_diameter_px",
        "centroid_y_px", "centroid_x_px", "seed_count",
        "growth_limit_reached", "bbox_min_y", "bbox_min_x",
        "bbox_max_y", "bbox_max_x",
    ]
    return roi_labels.astype(np.uint16), seeds_out, pd.DataFrame.from_records(records, columns=region_columns)


def normalize01(image: np.ndarray) -> np.ndarray:
    low,high=np.percentile(image,(1,99.8)); return np.clip((image-low)/(high-low),0,1) if high>low else np.zeros_like(image)


def save_overlay(image: np.ndarray, labels: np.ndarray, seeds: pd.DataFrame, regions: pd.DataFrame, path: Path, dataset: dict[str,Any], title: str) -> None:
    matplotlib_setup(dataset["output_dir"]); import matplotlib.pyplot as plt
    path.parent.mkdir(parents=True,exist_ok=True); rgb=np.repeat(normalize01(image)[...,None],3,axis=-1)
    class_by_id=dict(zip(regions.roi_id,regions.roi_class)) if not regions.empty else {}
    boundary=find_boundaries(labels,mode="outer"); rgb[boundary]=(1,1,0)
    fig,ax=plt.subplots(figsize=(7,10)); ax.imshow(rgb)
    for roi_id, group in seeds[seeds.roi_id>0].groupby("roi_id"):
        color="cyan" if class_by_id.get(roi_id)=="point_punctum" else "magenta"
        ax.scatter(group.x_px,group.y_px,s=18,marker="+",c=color,linewidths=.7)
    rejected=seeds[seeds.roi_id==0]
    if not rejected.empty: ax.scatter(rejected.x_px,rejected.y_px,s=18,marker="x",c="gray",linewidths=.7)
    ax.set_title(title); ax.axis("off"); fig.tight_layout(); fig.savefig(path,dpi=250,bbox_inches="tight"); plt.close(fig)


def process_file(row:dict[str,Any],cfg:dict[str,Any],dataset:dict[str,Any],context:SimpleNamespace,*,resume:bool,force:bool)->dict[str,Any]:
    corrected=drift_correction_paths(dataset["output_dir"],row["fov"])["stack"]
    if not corrected.exists(): raise FileNotFoundError(f"Run drift correction first: {corrected}")
    paths=roi_paths(dataset["output_dir"],row["fov"])
    if output_action(list(paths.values()),resume=resume,force=force)=="skip":
        regions=pd.read_csv(paths["regions"]); return {"filename":row["filename"],"fov":row["fov"],"status":"skipped","seeds":len(pd.read_csv(paths["seeds"])),"rois":len(regions),"point_puncta":int((regions.roi_class=="point_punctum").sum()),"extended_condensates":int((regions.roi_class=="extended_condensate").sum()),"overlay":str(paths["overlay"])}
    image=first_frames_mean(corrected,int(cfg["spotiflow"]["detection_frames"])); seeds=detect_seeds(image,cfg,context); labels,seeds,regions=grow_seed_regions(image,seeds,cfg)
    pixel_nm=float(row["pixel_size_um"])*1000
    if not seeds.empty: seeds["x_nm"]=seeds.x_px*pixel_nm; seeds["y_nm"]=seeds.y_px*pixel_nm
    paths["image"].parent.mkdir(parents=True,exist_ok=True); tifffile.imwrite(paths["image"],image,metadata={"axes":"YX"})
    paths["mask"].parent.mkdir(parents=True,exist_ok=True); tifffile.imwrite(paths["mask"],labels,metadata={"axes":"YX"})
    for table in (seeds,regions): table.insert(0,"filename",row["filename"]); table.insert(1,"fov",row["fov"])
    atomic_csv(seeds,paths["seeds"]); atomic_csv(regions,paths["regions"])
    save_overlay(image,labels,seeds,regions,paths["overlay"],dataset,f"{row['fov']} — Spotiflow 0.4 — {len(regions)} ROIs")
    return {"filename":row["filename"],"fov":row["fov"],"status":"complete","seeds":len(seeds),"rejected_seeds":int((seeds.roi_id==0).sum()),"rois":len(regions),
            "point_puncta":int((regions.roi_class=="point_punctum").sum()) if not regions.empty else 0,"extended_condensates":int((regions.roi_class=="extended_condensate").sum()) if not regions.empty else 0,
            "growth_limit_warnings":int(regions.growth_limit_reached.sum()) if not regions.empty else 0,"overlay":str(paths["overlay"])}


def save_stage_report(summary:pd.DataFrame,dataset:dict[str,Any])->None:
    matplotlib_setup(dataset["output_dir"]); import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    root=output_dirs(dataset["output_dir"])["roi_detection"]; report=root/"roi_detection_and_growth_QC_report.pdf"
    with PdfPages(report) as pdf:
        fig,ax=plt.subplots(figsize=(11,8.5)); ax.axis("off"); ax.text(.02,.98,f"ROI detection and growth QC — {dataset['name']}\n\n"+summary.to_string(index=False),va="top",family="monospace",fontsize=7); pdf.savefig(fig); plt.close(fig)
        for path in summary.overlay:
            fig,ax=plt.subplots(figsize=(8.5,11)); ax.imshow(plt.imread(path)); ax.axis("off"); pdf.savefig(fig); plt.close(fig)


def roi_detection_stage(cfg:dict[str,Any],dataset:dict[str,Any],manifest:pd.DataFrame,*,resume:bool=False,force:bool=False,max_files:int|None=None)->pd.DataFrame:
    rows=manifest[manifest.status=="accepted"].to_dict("records")[:max_files]; root=output_dirs(dataset["output_dir"])["roi_detection"]; root.mkdir(parents=True,exist_ok=True)
    context=load_spotiflow(cfg); records=[]
    with stage_timer(dataset["output_dir"],"03_roi_detection_and_growth",{"files":len(rows),"spotiflow_threshold":0.4}):
        for index,row in enumerate(rows,1): print(f"[ROI detection and growth {index}/{len(rows)}] {row['filename']}",flush=True); records.append(process_file(row,cfg,dataset,context,resume=resume,force=force))
        summary=pd.DataFrame.from_records(records); atomic_csv(summary,root/"roi_detection_and_growth_QC_summary.csv"); save_stage_report(summary,dataset)
    return summary


def main(argv:list[str]|None=None)->None:
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--config",required=True); parser.add_argument("--max-files",type=int); policy=parser.add_mutually_exclusive_group(); policy.add_argument("--resume",action="store_true"); policy.add_argument("--force",action="store_true")
    args=parser.parse_args(argv); cfg=load_config(args.config)
    for dataset in cfg["datasets"]: roi_detection_stage(cfg,dataset,load_manifest(dataset),resume=args.resume,force=args.force,max_files=args.max_files)


if __name__=="__main__": main()
