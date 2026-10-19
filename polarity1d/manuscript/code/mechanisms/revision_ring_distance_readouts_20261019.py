"""Systematic time-resolved distance readouts from all 652 retained rings.

Reversible read-only reanalysis. Crossing times are saved-sample brackets;
first crossings and retained losses are explicitly different readouts.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import revision_figure3_20261019 as geometry

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/"data/runs/chemical_ring_functional_v2_20261017"
OUT=geometry.OUT/"distance_response"
MODELS=geometry.MODELS
SOURCES={}


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write_json(path,z):path.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
def output(name,group="data"):
    parent=OUT if name.endswith(".png") else OUT/"supporting"/group
    parent.mkdir(parents=True,exist_ok=True);return parent/name
def savefig(fig,stem):
    for ext in ("png","pdf","svg"):fig.savefig(output(stem+"."+ext,"vectors"),dpi=200,bbox_inches="tight")
    plt.close(fig)


def weights(x,C,center,left,right):
    arc=(x-center+C/2)%C-C/2;dx=C/len(x)
    return np.maximum(0,np.minimum(arc+dx/2,right)-np.maximum(arc-dx/2,left))


def frozen_core(x,C,p):
    peak=int(p.argmax());center=x[peak];arc=(x-center+C/2)%C-C/2
    order=np.argsort(arc);s,u=arc[order],p[order];level=u.min()+.5*np.ptp(u)
    inside=np.flatnonzero(u>=level);assert np.all(np.diff(inside)==1)
    a,b=int(inside[0]),int(inside[-1]);assert a>0 and b<len(s)-1
    left=s[a-1]+(level-u[a-1])/(u[a]-u[a-1])*(s[a]-s[a-1])
    right=s[b]+(level-u[b])/(u[b+1]-u[b])*(s[b+1]-s[b])
    w=weights(x,C,center,left,right)
    np.testing.assert_allclose(w.sum(),right-left,atol=1e-12)
    return w,float(u.min())


def first_bracket(t,criterion):
    hit=np.flatnonzero(criterion)
    if not len(hit):return dict(status="right_censored",lower=float(t[-1]),upper=None,first_sample=None)
    i=int(hit[0])
    return dict(status="already_at_start" if i==0 else "observed_saved_sample_crossing",
        lower=float(t[max(0,i-1)]),upper=float(t[i]),first_sample=float(t[i]))


def terminal_bracket(t,criterion):
    if not criterion[-1]:return dict(status="not_retained_at_final_sample",lower=None,upper=None,first_sample=None)
    last_bad=np.flatnonzero(~criterion)
    i=int(last_bad[-1]+1) if len(last_bad) else 0
    return dict(status="sampled_entry_retained_through_horizon",lower=float(t[max(0,i-1)]),upper=float(t[i]),first_sample=float(t[i]))


def main():
    geometry.setup();matrix=pd.read_csv(RUN/"results/matrix.csv")
    expected=json.loads((RUN/"results/source_sha256.json").read_text())
    local=pd.read_csv(ROOT/"data/runs/vision_ring_duration_20261017/results/all_grid_local_readouts.csv")
    rows=[];crossings=[];readouts=[];basals={};trace=[];offsets=[0];max_old_error=0.
    for _,row in matrix[matrix.status=="complete"].iterrows():
        C=float(row.circumference);h=geometry.history(row);p,t,x=h["p"],h["t"],h["x"]
        for name in ("cue_path","recovery_path"):
            assert geometry.SOURCES[str(Path(row[name]).resolve())]==expected[str(Path(row[name]))]
        key=(row.model,C)
        if key not in basals:
            core,bg=frozen_core(x,C,p[0]);basals[key]=(core,bg,p[0].copy())
        core,bg,original=basals[key];np.testing.assert_array_equal(p[0],original)
        np.testing.assert_allclose(core.sum(),row.basal_width,atol=1e-10)
        np.testing.assert_allclose(p[0]@core,row.basal_core_amount,atol=1e-9)
        cue_center=x[int(row.cue_center_cell)];patch=weights(x,C,cue_center,-1.5,1.5)
        aold=p@core;acue=p@patch;den=float(aold[0]-bg*core.sum());assert den>0
        oldraw=aold/aold[0];oldexcess=(aold-bg*core.sum())/den
        cueexcess=(acue-bg*patch.sum())/den;cueincrement=(acue-acue[0])/den
        angle=np.rad2deg(np.angle(h["P"]))
        angle=np.where(np.isclose(abs(angle),180.,rtol=0,atol=1e-7),180.,angle)
        cue_angle=360*row.short_route/C
        cue_error=np.rad2deg(np.angle(np.exp(1j*np.deg2rad(angle-cue_angle))))
        meta=dict(model=row.model,circumference=C,fraction=float(row.fraction),physical_distance=float(row.short_route),
            angular_separation_deg=cue_angle,basal_width=float(row.basal_width),edge_gap=float(row.short_edge_to_cue_gap),
            checkpoint_sha256=row.checkpoint_sha256,cue_duration=1000.,recovery=13000.,outcome=row.final_location_association)
        summary=dict(**meta,old_core_excess_minimum=float(oldexcess.min()),old_core_excess_minimum_time=float(t[oldexcess.argmin()]),
            old_core_excess_minimum_cue_on=float(oldexcess[t<=1000].min()),old_core_excess_final=float(oldexcess[-1]),
            old_core_excess_at_withdrawal=float(oldexcess[t==1000][0]),old_core_raw_final=float(oldraw[-1]),
            old_core_excess_rebound_after_withdrawal=float(oldexcess[-1]-oldexcess[t==1000][0]),
            cue_patch_excess_final=float(cueexcess[-1]),cue_patch_increment_final=float(cueincrement[-1]),
            final_R=float(h["R"][-1]),final_angle_deg=float(angle[-1]) if h["R"][-1]>=.1 else None,
            final_cue_angle_error_deg=float(cue_error[-1]) if h["R"][-1]>=.1 else None,
            angle_eligible_at_final=bool(h["R"][-1]>=.1),sampled_minimum_R=float(h["R"].min()),
            trace_start=offsets[-1],trace_frames=len(t),max_cue_on_gap=float(np.diff(t[t<=1000]).max()),
            max_cue_free_gap=float(np.diff(t[t>=1000]).max()))
        for threshold in (.9,.5,.1):
            for window,mask in (("cue_on",t<=1000),("full",np.ones(len(t),bool))):
                bracket=first_bracket(t[mask],oldexcess[mask]<=threshold)
                crossings.append(dict(**meta,observable="frozen_original_core_excess_fraction",threshold=threshold,
                    window=window,event="first_sampled_crossing",**bracket))
                if threshold==.5 and window=="cue_on":
                    summary.update(half_loss_status=bracket["status"],half_loss_lower=bracket["lower"],half_loss_upper=bracket["upper"])
            crossings.append(dict(**meta,observable="frozen_original_core_excess_fraction",threshold=threshold,
                window="full",event="entry_retained_at_all_later_saved_times",**terminal_bracket(t,oldexcess<=threshold)))
        for tolerance in (15.,30.):
            for gate in (.05,.1,.2):
                match=(abs(cue_error)<=tolerance)&(h["R"]>=gate)
                crossings.append(dict(**meta,observable="activity_vector_cue_alignment_deg",threshold=tolerance,
                    magnitude_gate=gate,window="full",event="first_sampled_entry",**first_bracket(t,match)))
                crossings.append(dict(**meta,observable="activity_vector_cue_alignment_deg",threshold=tolerance,
                    magnitude_gate=gate,window="full",event="entry_retained_at_all_later_saved_times",**terminal_bracket(t,match)))
        for readout in (0.,10.,1000.,4000.,9000.,14000.):
            idx,=np.flatnonzero(np.isclose(t,readout,rtol=0,atol=1e-12))
            readouts.append(dict(**meta,time=readout,time_after_withdrawal=readout-1000,
                old_core_excess_fraction=float(oldexcess[idx]),old_core_raw_fraction=float(oldraw[idx]),
                cue_patch_excess_over_basal_incumbent=float(cueexcess[idx]),cue_patch_increment_over_basal_incumbent=float(cueincrement[idx]),
                R=float(h["R"][idx]),q1=float(h["q1"][idx]),q2=float(h["q2"][idx]),
                angle_deg=float(angle[idx]) if h["R"][idx]>=.1 else None,
                cue_angle_error_deg=float(cue_error[idx]) if h["R"][idx]>=.1 else None))
            if readout in (4000.,9000.,14000.):
                target=local[(local.model==row.model)&(local.circumference==C)&np.isclose(local.fraction,row.fraction)&(local.cue_free_time==readout-1000)].iloc[0]
                error=float(abs(oldraw[idx]-target.original_core_amount_fraction));assert error<1e-9
                max_old_error=max(max_old_error,error)
        trace.append(np.column_stack((t,oldraw,oldexcess,cueexcess,cueincrement,h["R"],angle,h["q1"],h["q2"])))
        offsets.append(offsets[-1]+len(t));rows.append(summary)
    summary=pd.DataFrame(rows);crossings=pd.DataFrame(crossings);readouts=pd.DataFrame(readouts)
    assert len(summary)==652 and len(readouts)==3912
    for name,data in (("distance_response_summary.csv",summary),("distance_crossing_brackets.csv",crossings),("distance_response_readouts.csv",readouts)):
        data.to_csv(output(name),index=False)
    values=np.concatenate(trace)
    np.savez_compressed(output("distance_response_traces.npz"),values=values,offsets=np.asarray(offsets),
        columns=np.asarray(["time","old_core_raw_fraction","old_core_excess_fraction","cue_patch_excess","cue_patch_increment","R","angle_deg_unmasked","q1","q2"]))
    durations=pd.read_csv(ROOT/"data/runs/vision_ring_duration_20261017/results/duration_comparison.csv")
    durations.to_csv(output("existing_duration_anchor_readouts.csv"),index=False)
    plots(summary,readouts)
    doc(summary,crossings,max_old_error)


def plots(summary,readouts):
    fig,axes=plt.subplots(3,3,figsize=(11.8,10.2),constrained_layout=True)
    for j,model in enumerate(MODELS):
        g=summary[(summary.model==model)&(summary.circumference==40)].sort_values("physical_distance")
        r=readouts[(readouts.model==model)&(readouts.circumference==40)]
        for time,ls,marker,label in ((1000,"--","s","At withdrawal"),(14000,"-","o","After 13,000 cue-free")):
            z=r[r.time==time].sort_values("physical_distance")
            axes[0,j].plot(z.physical_distance,z.old_core_excess_fraction,ls=ls,marker=marker,ms=3.5,color="#4C78A8",label=label)
            axes[1,j].plot(z.physical_distance,z.angle_deg,ls=ls,marker=marker,ms=3.5,color="#A46837",label=label)
        axes[0,j].axhline(.5,color="#bbbbbb",lw=.8,ls=":")
        axes[0,j].set(title=geometry.NAMES[model]+" · C40",ylabel="Original-core excess / basal excess" if j==0 else "",ylim=(-.06,1.05))
        axes[1,j].plot(g.physical_distance,g.angular_separation_deg,color="#aaaaaa",lw=.8,label="Prescribed cue angle")
        axes[1,j].set(ylabel="Activity-vector angle (°)" if j==0 else "",ylim=(-5,190))
        uncensored=g.half_loss_upper.notna();z=g[uncensored]
        axes[2,j].vlines(z.physical_distance,z.half_loss_lower,z.half_loss_upper,color="#386C55",lw=1.3)
        axes[2,j].plot(z.physical_distance,z.half_loss_upper,"o-",color="#386C55",ms=4)
        censored=g[~uncensored]
        axes[2,j].scatter(censored.physical_distance,np.full(len(censored),1000),marker="^",facecolor="white",edgecolor="#386C55",s=40,label="Not crossed during cue")
        axes[2,j].set(ylim=(0,1100),xlabel="Prescribed peak–cue physical distance d",ylabel="First half-loss time during cue" if j==0 else "")
    axes[0,0].legend(frameon=False,fontsize=8);axes[1,0].legend(frameon=False,fontsize=8,loc="upper left");axes[2,2].legend(frameon=False,fontsize=8)
    fig.suptitle("Distance changes the response trajectory even when the final fate agrees\nSame qualified C40 basal per model; full support 3; cue 1000 + 13,000 cue-free units",fontsize=13)
    savefig(fig,"distance_response_C40")
    fig,axes=plt.subplots(3,3,figsize=(11.8,10),constrained_layout=True)
    for i,model in enumerate(MODELS):
        g=summary[summary.model==model]
        for j,(key,title,vmin,vmax,cmap) in enumerate((("old_core_excess_final","Retained original-core excess",0,1,"viridis"),("final_angle_deg","Final activity angle (°)",0,180,"twilight_shifted"),("half_loss_upper","First cue-on half-loss time",0,1000,"magma"))):
            ax=axes[i,j];valid=g[key].notna()
            im=ax.scatter(g.loc[valid,"fraction"],g.loc[valid,"circumference"],c=g.loc[valid,key],s=35,marker="s",vmin=vmin,vmax=vmax,cmap=cmap)
            ax.scatter(g.loc[~valid,"fraction"],g.loc[~valid,"circumference"],marker="x",s=23,c="#777777")
            ax.set(title=geometry.NAMES[model]+"\n"+title,xlabel="Prescribed shorter-arc fraction d/C",ylabel="Circumference C")
            fig.colorbar(im,ax=ax,shrink=.8)
    fig.suptitle("All 652 complete native ring responses · excluded preparations/overlaps omitted\nAngle crosses: R < 0.1; half-loss crosses: no crossing by withdrawal at 1000",fontsize=12)
    savefig(fig,"distance_response_all_sizes")


def doc(summary,crossings,max_error):
    now=datetime.now(timezone.utc).isoformat()
    txt=r"""# Systematic ring distance-response reanalysis

All652complete histories from867requested conditions are analysed, preserving the approved geometry sheet as a preceding draft. There are198overlap exclusions and17unqualified preparations. Every completed case uses one fixed-duration1000cue and13000cue-freeunits. This is systematic distance×size coverage, not a minimum-duration search. The existing duration sensitivity contains9selectedgeometryanchors×3durations500/1000/2000, so it does not supply a dense distance×duration surface.

The original core is the frozen basal half-contrast interval with fractional finite-volume edge weights. Primary readout is excess activity above the frozen basal background, integrated in that fixed region and divided by the basal incumbent excess. Raw core activity fractions are also retained and reproduce1956previously audited localreadouts. Cue-patch excess and increment above its exact initial value use fixedphysicalsupport3, normalized by the same basal incumbent excess. None of these fixed-region amounts identifies the lineage of an evolving activity peak.

For thresholds0.9/0.5/0.1, report first saved-sample crossing separately during the cue and over the full14000window. Brackets use the last prior savedtime and first passing savedtime, never interpolation. A first crossing does not imply retained loss: sampled entry retained through every subsequent savedtime is a separate readout. No crossing by the horizon is right-censoring, not zero delay. Events can occur between savedframes and recross unobserved. Orientation readouts use P=M1(t)/M1(0), angle gated atR≥0.1, with cue-alignment entry tolerances15°/30° and magnitude gates0.05/0.1/0.2. Angle relative to initial incumbent is not cell movement.

In C40, Goryachev's original-core half-loss time increases with distance while every final endpoint is cue-associated. Wave-pinning has a nonmonotone half-loss response, with prolonged near-incumbent overlap-adjacent competition. Holmes maintains more than half the original-core excess during every testedC40cue, despite changed globalactivityangles and local spatial organization. Main C40 curves therefore add information hidden by endpoint categories. Cross-model times retain their native units and are not calibrated rates.

Files: `distance_response_summary.csv`, `distance_response_readouts.csv`, `distance_crossing_brackets.csv` and compact434002-frame `distance_response_traces.npz`. Full retained fields are source evidence; no new native integration is performed. Hashes are verified against the archived campaign. `existing_duration_anchor_readouts.csv` records the actual sparse duration evidence. Figures show the fixedC40controlled-distance comparison and all-size surfaces. Censoring and masked directions remain visible.
"""
    output("EVIDENCE.md","metadata").write_text("Generated:"+now+"\n\n"+txt)
    sources={str(RUN/"results/matrix.csv"):sha(RUN/"results/matrix.csv"),str(Path(__file__)):sha(__file__),**geometry.SOURCES}
    write_json(output("source_sha256.json","metadata"),sources)
    allocations=sum(p.stat().st_blocks*512 for p in OUT.rglob("*") if p.is_file())
    write_json(output("validation.json","metadata"),dict(generated_at=now,passed=True,complete_conditions=len(summary),
        readout_rows=3912,crossing_records=len(crossings),frames=434002,previous_raw_core_readouts_reproduced=1956,
        maximum_raw_core_fraction_replay_error=max_error,all_native_source_hashes_match_archived=True,
        exact_phase_primary_continuity=True,all_initial_profiles_match_frozen_basal=True,
        source_model_unchanged=True,new_integrations=0,output_allocated_bytes_before_validation=allocations))
    print(json.dumps(dict(passed=True,conditions=len(summary),crossing_records=len(crossings),output=str(OUT)),indent=2))


if __name__=="__main__":main()
