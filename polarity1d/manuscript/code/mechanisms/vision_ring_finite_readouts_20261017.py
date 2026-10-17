"""Read-only finite-window local classifications and duration-control table."""
from pathlib import Path
import json
import os
import sys
for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[k]="1"
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import chemical_ring_functional_v2_20261017 as c
import chemical_ring_functional_collect_20261017 as oldcollect

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"data/runs/chemical_ring_functional_v2_20261017"
OUT=ROOT/"data/runs/vision_ring_duration_20261017/results"


def main():
    matrix=pd.read_csv(SOURCE/"results/matrix.csv");rows=[];basals={};hashes={}
    for _,x in matrix[matrix.status=="complete"].iterrows():
        key=(x.model,int(x.circumference))
        if key not in basals:
            s=c.read(Path(x.status_path));checkpoint=Path(s["checkpoint"])
            assert c.sha(checkpoint)==x.checkpoint_sha256
            with np.load(checkpoint) as z:basal=z["state"].copy()
            grid=c.ring.RingGrid(int(x.circumference),int(x.cells));core=c.previous.original_core(grid,basal[0])
            basals[key]=(basal,grid,core)
        basal,grid,core=basals[key];path=Path(x.recovery_path)
        hashes[str(path)]=c.sha(path)
        with np.load(path) as z:
            times=z["time"];values=z["states"][:,0] if "states" in z.files else z["primary"]
            for t in (3000.,8000.,13000.):
                i=int(np.flatnonzero(times==t)[0]);u=values[i]
                identity=oldcollect.identity(u,core,int(x.cue_center_cell),c.spatial(grid,u,core,int(x.cue_center_cell)))
                if t==13000:assert identity["location_association"]==x.final_location_association
                rows.append(dict(model=x.model,circumference=x.circumference,fraction=x.fraction,
                    cue_duration=1000.,cue_free_time=t,time_from_onset=t+1000,
                    resolved_peak_count=identity["resolved_peak_count"],location_association=identity["location_association"],
                    distinct_old_peak_available=identity["distinct_old_peak_available"],
                    original_core_amount=float(u@core["weights"]),
                    original_core_amount_fraction=float(u@core["weights"]/(basal[0]@core["weights"])),
                    recovery_path=str(path),recovery_sha256=hashes[str(path)],checkpoint_sha256=x.checkpoint_sha256))
    q=pd.DataFrame(rows);q.to_csv(OUT/"all_grid_local_readouts.csv",index=False)
    counts=[];exceptions=[]
    final=q[q.cue_free_time==13000].set_index(["model","circumference","fraction"])
    for t in (3000.,8000.):
        a=q[q.cue_free_time==t].set_index(["model","circumference","fraction"])
        for model in c.KEYS:
            current=a.loc[model];end=final.loc[model]
            changed=current.location_association!=end.location_association
            peak_changed=current.resolved_peak_count!=end.resolved_peak_count
            counts.append(dict(model=model,cue_free_time=t,reference_cue_free_time=13000,conditions=len(current),
                geographic_class_changed=int(changed.sum()),resolved_peak_count_changed=int(peak_changed.sum())))
            for key in current.index[changed]:
                x=current.loc[key];y=end.loc[key]
                exceptions.append(dict(model=model,circumference=key[0],fraction=key[1],cue_free_time=t,
                    early_class=x.location_association,final_class=y.location_association,
                    early_peaks=x.resolved_peak_count,final_peaks=y.resolved_peak_count,recovery_path=x.recovery_path))
    pd.DataFrame(counts).to_csv(OUT/"readout_change_counts.csv",index=False)
    pd.DataFrame(exceptions).to_csv(OUT/"readout_changed_cells.csv",index=False)
    controls=pd.read_csv(OUT/"duration_readouts.csv")
    moments=pd.read_csv(ROOT/"data/runs/vision_ring_moments_20261017/results/readout_sensitivity.csv")
    original=[]
    for model,C,f in controls[["model","circumference","fraction"]].drop_duplicates().itertuples(index=False,name=None):
        local=q[(q.model==model)&(q.circumference==C)&(q.fraction==f)]
        for _,r in local.iterrows():
            mv=moments[(moments.model==model)&(moments.circumference==C)&(moments.fraction==f)&(moments.time_after_withdrawal==r.cue_free_time)].iloc[0]
            original.append(dict(model=model,circumference=C,fraction=f,cue_duration=1000.,recovery=13000.,cue_free_time=r.cue_free_time,
                original_core_amount=r.original_core_amount,original_core_amount_fraction=r.original_core_amount_fraction,
                Px=mv.P_x,Py=mv.P_y,R=mv.R,angle_deg=mv.angle_deg_gate_0p1,
                resolved_peak_count=r.resolved_peak_count,location_association=r.location_association,
                history_path=r.recovery_path,history_sha256=r.recovery_sha256))
    combined=pd.concat([controls,pd.DataFrame(original)],ignore_index=True).sort_values(["model","circumference","fraction","cue_duration","cue_free_time"])
    assert len(combined)==81
    combined.to_csv(OUT/"duration_comparison.csv",index=False)
    endpoint=combined[combined.cue_free_time==13000]
    stable=[]
    for key,group in endpoint.groupby(["model","circumference","fraction"]):
        stable.append(dict(model=key[0],circumference=key[1],fraction=key[2],
            same_class_across_500_1000_2000=group.location_association.nunique()==1,
            same_peak_count=group.resolved_peak_count.nunique()==1,R_range=float(group.R.max()-group.R.min())))
    pd.DataFrame(stable).to_csv(OUT/"duration_class_robustness.csv",index=False)
    c.write(OUT/"finite_readout_audit.json",dict(passed=True,complete=True,all_652_endpoint_classes_reproduce=True,
        original_conditions=652,local_readout_rows=len(q),duration_comparison_rows=81,
        readout_counts=counts,duration_anchor_classes=stable,source_recovery_sha256=hashes,
        script_sha256=c.sha(Path(__file__)),interpretation="Descriptive spatial-association classes at recorded times, not tracked fate, asymptotic stability or a threshold campaign."))
    print(json.dumps(dict(complete=True,counts=counts,all9duration_classes_stable=all(r["same_class_across_500_1000_2000"] for r in stable)),indent=2))


if __name__=="__main__":main()
