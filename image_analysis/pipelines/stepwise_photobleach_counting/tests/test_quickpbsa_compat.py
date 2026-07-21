from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from quickpbsa_compat import result_summary, run_quickpbsa


def test_quickpbsa_runs_on_numpy2_and_python313(tmp_path: Path):
    rng=np.random.default_rng(4); frames=180; traces=[]
    specifications = [([30,75,125],90.),([20,50,85,120,150],100.),([15,35,60,85,110,135,160],110.)]
    for positions, step_size in specifications:
        level=len(positions)*step_size; values=np.empty(frames); index=0
        for frame in range(frames):
            while index<len(positions) and frame>=positions[index]: level-=step_size; index+=1
            values[frame]=level+rng.normal(0,2)
        traces.append(values)
    infile=tmp_path/"synthetic_difference.csv"; pd.DataFrame(traces).to_csv(infile,index=False)
    settings={"maxiter":20,"num_cores":2,"preliminary":{"crop":False,"bg_frames":20},"filter":{"subtracted":True,"percentile_step":100,"length_laststep":10},"refinement":{}}
    result=run_quickpbsa(infile,tmp_path/"result",40,settings)
    summary=result_summary(result)
    assert len(summary)==3 and "photobleaching_step_count" in summary
    accepted = summary[summary.flag == 1]
    assert set(accepted.photobleaching_step_count.dropna().astype(int))=={3,5}
    assert summary.loc[summary.flag != 1,"photobleaching_step_count"].isna().all()
