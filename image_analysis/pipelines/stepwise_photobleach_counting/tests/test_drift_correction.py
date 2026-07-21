from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
from scipy.ndimage import gaussian_filter, shift as ndi_shift

from conftest import write_stack
from step02_drift_correction import estimate_block_shifts, frame_shifts, residual_shifts, write_corrected_stack


def synthetic_image(shape=(64,72)):
    image=np.zeros(shape,np.float32)
    for y,x,value in [(15,18,1000),(43,51,800),(28,34,600)]: image[y,x]=value
    return gaussian_filter(image,1.2)+100


def settings():
    return {"block_frames":5,"highpass_sigma_px":6,"upsample_factor":10,"max_abs_shift_px":10,
            "interpolation_order":1,"compression":"zlib","compression_level":1}


def test_phase_cross_correlation_recovers_and_corrects_drift(tmp_path: Path):
    base=synthetic_image(); offsets=[(0,0),(1.2,-.8),(2.0,-1.5)]
    frames=[]
    for offset in offsets:
        frames.extend([ndi_shift(base,offset,order=1,mode="nearest") for _ in range(5)])
    stack=np.clip(np.stack(frames),0,65535).astype(np.uint16); source=tmp_path/"source.tif"; write_stack(source,stack)
    table,blocks=estimate_block_shifts(source,settings())
    assert np.allclose(table.interpolated_y_shift_px.to_numpy(),-np.asarray([value[0] for value in offsets]),atol=.2)
    assert np.allclose(table.interpolated_x_shift_px.to_numpy(),-np.asarray([value[1] for value in offsets]),atol=.2)
    ry,rx=residual_shifts(blocks,table,settings()); assert np.percentile(np.hypot(ry,rx),95)<=.2
    destination=tmp_path/"corrected.tif"; margin,error=write_corrected_stack(source,destination,frame_shifts(table,len(stack)),settings())
    assert destination.exists() and margin>=3 and error<2
    corrected=tifffile.imread(destination); assert corrected.shape[0]==len(stack)
