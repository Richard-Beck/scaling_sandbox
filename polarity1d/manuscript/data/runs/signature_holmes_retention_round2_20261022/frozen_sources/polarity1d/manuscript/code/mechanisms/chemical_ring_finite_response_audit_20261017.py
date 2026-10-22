"""Retain matched no-cue controls and verify finite-time cue increments."""
from pathlib import Path
import argparse
import csv
import json
import os
import sys
import time
for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[k]="1"
import numpy as np


def main(run):
    assert os.environ.get("SLURM_JOB_ID"),"Native integration requires SLURM"
    sys.path.insert(0,str(run/"frozen_sources/code/mechanisms"))
    import chemical_ring_functional_v2_20261017 as c
    m=c.check_sources(run);r=run/"finite_response_audit";records=[];sources={};worst=0.
    for job in m["jobs"]:
        s=c.read(run/"workers"/f'{job["slot"]:02d}'/"status.json")
        if s["status"]!="complete":continue
        checkpoint=Path(s["checkpoint"]);assert c.sha(checkpoint)==s["checkpoint_sha256"]
        with np.load(checkpoint) as z:original=z["state"].copy()
        model=c.ring.native_model(job["model"]);g=c.ring.RingGrid(job["circumference"],job["cells"])
        times,states,stats=c.integrate(model,g,original,np.zeros(g.cells),10.,np.array([0.,1.,10.]),
            c.ring.pool_amounts(model,g,original),time.monotonic()+3300)
        path=r/f'{job["model"]}_C{job["circumference"]}.npz'
        np.savez_compressed(path,time=times,states=states,x=g.x)
        sources[str(path)]=c.sha(path)
        for response in s["responses"]:
            if response["status"]!="complete":continue
            cuepath=Path(next(p["path"] for p in response["phases"] if p["phase"]=="cue"))
            with np.load(cuepath) as z:
                ts=z["time"] if "states" in z.files else z["snapshot_times"]
                values=z["states"] if "states" in z.files else z["snapshot_states"]
                weight=c.previous.weights_about(g,float(g.x[response["cue_center_cell"]]),-1.5,1.5)
                for t in (1.,10.):
                    cued=values[np.flatnonzero(np.isclose(ts,t,rtol=0,atol=1e-10))[0]]
                    control=states[np.flatnonzero(times==t)[0]]
                    observed=float((cued[0]-control[0])@weight/weight.sum())
                    recorded=response["gain"][f"observed_patch_mean_extra_at_{t:g}"]
                    error=abs(observed-recorded);worst=max(worst,error)
                    assert error<1e-11
                    records.append(dict(model=job["model"],circumference=job["circumference"],fraction=response["fraction"],
                        time=t,observed_patch_mean_extra=observed,worker_recorded_extra=recorded,
                        absolute_error=error,control_path=str(path),control_sha256=sources[str(path)],
                        checkpoint_sha256=s["checkpoint_sha256"],cue_path=str(cuepath)))
    with (r/"finite_response_validation.csv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    c.write(r/"audit.json",dict(generated_at=c.now(),passed=True,matched_no_cue_controls=len(sources),
        finite_time_readouts=len(records),max_absolute_worker_difference=worst,source_sha256=sources,
        frozen_native_sources_verified=True,runner_sha256=c.sha(Path(__file__))))
    print(json.dumps(dict(passed=True,controls=len(sources),readouts=len(records),max_error=worst)))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--run",type=Path,required=True)
    main(parser.parse_args().run.resolve())
