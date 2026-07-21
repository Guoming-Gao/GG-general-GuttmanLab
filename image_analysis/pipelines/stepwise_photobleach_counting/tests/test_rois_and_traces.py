from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from skimage.draw import disk

from step03_detect_and_grow_rois import grow_seed_regions
from step04_extract_background_corrected_traces import build_roi_indices, extract_traces


def roi_config():
    return {
        "roi_growth":{"smooth_sigma_px":0.5,"foreground_sigma":3,"maximum_radius_px":15,"outer_background_width_px":5,"point_max_equivalent_diameter_px":5,"minimum_area_px":2},
        "trace_extraction":{"point_radius_px":2,"point_background_inner_radius_px":3.5,"point_background_outer_radius_px":5,
                            "point_minimum_neighbor_distance_px":3.5,"mask_background_gap_px":2,"mask_background_width_px":2},
    }


def test_seed_growth_classifies_point_and_extended_regions():
    rng=np.random.default_rng(2); image=rng.normal(100,1,(64,64)).astype(np.float32)
    image[15,15]=170; rr,cc=disk((42,42),7,shape=image.shape); image[rr,cc]=150
    seeds=pd.DataFrame({"seed_id":[1,2,3],"y_px":[15.,42.,44.],"x_px":[15.,42.,44.],"spotiflow_probability":[.8,.9,.7],"spotiflow_intensity":[170,150,150]})
    labels,seeds_out,regions=grow_seed_regions(image,seeds,roi_config())
    assert set(regions.roi_class)=={"point_punctum","extended_condensate"}
    extended=regions[regions.roi_class=="extended_condensate"].iloc[0]
    assert extended.seed_count==2 and extended.area_px>50
    assert seeds_out.roi_id.nunique()==2 and labels.max()==2


def test_seed_growth_handles_an_image_with_no_spotiflow_seeds():
    image=np.full((32,40),100,np.float32)
    seeds=pd.DataFrame(columns=["seed_id","y_px","x_px","spotiflow_probability","spotiflow_intensity"])
    labels,seeds_out,regions=grow_seed_regions(image,seeds,roi_config())
    assert labels.max()==0 and seeds_out.empty and regions.empty
    assert {"roi_id","roi_class","area_px"}.issubset(regions.columns)


def test_integrated_trace_is_area_independent(tmp_path: Path):
    labels=np.zeros((40,48),np.uint16); labels[10:12,10:12]=1; labels[24:31,28:36]=2
    regions=pd.DataFrame([{"roi_id":1,"roi_class":"point_punctum","area_px":4,"equivalent_diameter_px":2.3,"centroid_y_px":10.5,"centroid_x_px":10.5,"seed_count":1},
                          {"roi_id":2,"roi_class":"extended_condensate","area_px":56,"equivalent_diameter_px":8.4,"centroid_y_px":27.,"centroid_x_px":31.5,"seed_count":1}])
    seeds=pd.DataFrame({"seed_id":[1,2],"y_px":[10.5,27.],"x_px":[10.5,31.5],"roi_id":[1,2]})
    entries,signal_mask,bg_mask=build_roi_indices(labels,seeds,regions,roi_config()); accepted=[entry for entry in entries if entry["trace_status"]=="accepted"]
    assert len(accepted)==2 and np.all(bg_mask[labels>0]==0)
    frames=[]
    for level in [300.,200.,100.,0.]:
        image=np.full(labels.shape,50.,np.float32)
        for entry in accepted: image.ravel()[entry["signal_indices"]]+=level/len(entry["signal_indices"])
        frames.append(image)
    stack=tmp_path/"traces.tif"
    with tifffile.TiffWriter(stack) as writer:
        for frame in frames: writer.write(frame.astype(np.float32),metadata=None)
    traces=extract_traces(stack,entries)
    assert np.allclose(traces["integrated_difference"][0],[300,200,100,0],atol=.01)
    assert np.allclose(traces["integrated_difference"][1],[300,200,100,0],atol=.01)
