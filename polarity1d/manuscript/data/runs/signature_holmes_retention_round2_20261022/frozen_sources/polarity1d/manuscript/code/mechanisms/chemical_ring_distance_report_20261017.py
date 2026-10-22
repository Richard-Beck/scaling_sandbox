"""Actual-data report for the first fixed-dose chemical ring distance pilot."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
import sys

for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key]="1"
sys.dont_write_bytecode=True
import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from length_dev_output_layout import dev_output_path,dev_artifact_key
import chemical_ring_distance_pilot_20261017 as pilot

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/"data/runs/chemical_ring_distance_pilot_20261017"
OUT=ROOT/"figures/length_mechanisms_20261008/dev/chemical_ring_distance_20261017"
NAMES={"wave_pinning":"Wave-pinning","goryachev":"Goryachev","holmes_model3":"Holmes 3"}
SOURCES={}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def source(path,expected=None):
    digest=sha(path)
    if expected is not None:assert digest==expected,path
    SOURCES[str(Path(path).relative_to(ROOT))]=digest
    return Path(path)
def read(path):return json.loads(source(path).read_text())
def csvwrite(name,rows):
    path=dev_output_path(OUT,name)
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    return path
def jsonwrite(name,value):
    path=dev_output_path(OUT,name);path.write_text(json.dumps(value,indent=2)+"\n");return path


def load_and_check():
    manifest=read(RUN/"manifest.json");launch=read(RUN/"launch.json")
    assert sha(RUN/"manifest.json")==launch["manifest_sha256"]
    for path,digest in manifest["source_sha256"].items():source(RUN/"frozen_sources"/path,digest)
    grid=pilot.ring.RingGrid(40,800);data={};summary=[];traces=[];audits=[]
    for job in manifest["jobs"]:
        base=manifest["models"][job["model"]]
        with np.load(source(RUN/base["checkpoint"],base["checkpoint_sha256"])) as z:original=z["state"].copy()
        record=read(RUN/"workers"/f'{job["slot"]:02d}'/"status.json")
        assert record["complete"]
        model=pilot.ring.native_model(job["model"]);core=pilot.original_core(grid,original[0])
        cuew=pilot.weights_about(grid,float(grid.x[job["cue_center_cell"]]),-1.5,1.5)
        scales=np.array([base["fixed_field_scales"][name] for name in base["state_names"]])
        initial_pools=pilot.ring.pool_amounts(model,grid,original)
        case=dict(job=job,record=record,original=original,core=core,phases={})
        previous=original.copy();pool_error=0.;minimum=1e10
        for phase in record["phases"]:
            path=source(RUN/phase["path"],phase["sha256"])
            with np.load(path) as z:
                t,states=z["time"].copy(),z["states"].copy()
                np.testing.assert_array_equal(states[0],previous)
                np.testing.assert_array_equal(z["original_core_weights"],core["weights"])
                assert np.isfinite(states).all();minimum=min(minimum,float(states.min()))
                pools=np.column_stack([states[:,list(group),:].sum(axis=(1,2))*grid.dx for group in model.conserved_groups])
                pool_error=max(pool_error,float(np.max(np.abs(pools/initial_pools-1))))
                np.testing.assert_allclose(z["pools"],pools,rtol=2e-14,atol=2e-13)
                np.testing.assert_allclose(z["old_core_amount"],states[:,0]@core["weights"],rtol=2e-14,atol=2e-13)
                np.testing.assert_allclose(z["cue_region_amount"],states[:,0]@cuew,rtol=2e-14,atol=2e-13)
                np.testing.assert_array_equal(z["peak_value"],states[z["peak_frame"],0,z["peak_cell"]])
                np.testing.assert_array_equal(z["peak_left_valley_value"],states[z["peak_frame"],0,z["peak_left_valley_cell"]])
                np.testing.assert_array_equal(z["peak_right_valley_value"],states[z["peak_frame"],0,z["peak_right_valley_cell"]])
                assert np.max(np.diff(t))<=5.00000001
                assert np.max(np.diff(t[t<=10]))<=0.10000001
                assert np.max(np.diff(t[(t>=10)&(t<=100)]))<=1.00000001
                p=phase["phase"]
                if phase["status"]=="complete":assert t[-1]==(1000 if p=="cue" else 13000)
                scalar={key:z[key].copy() for key in ("old_core_mean","old_core_amount","cue_region_mean","cue_region_amount","global_max_cell","harmonic_real","harmonic_imag","harmonic_magnitude","peaks_over_native_prominence_reference","frozen_basal_half_contrast_components")}
                peak_unwrapped=np.unwrap(2*np.pi*(scalar["global_max_cell"]-core["peak_cell"])/grid.cells)*40/(2*np.pi)
                moment_unwrapped=np.unwrap(np.arctan2(scalar["harmonic_imag"],scalar["harmonic_real"]))*40/(2*np.pi)
                for i,time_value in enumerate(t):
                    traces.append(dict(model=job["model"],signed_distance=job["signed_distance"],role=job["role"],phase=p,time_in_phase=float(time_value),
                        old_core_mean=float(scalar["old_core_mean"][i]),old_core_integral=float(scalar["old_core_amount"][i]),
                        cue_region_mean=float(scalar["cue_region_mean"][i]),cue_region_integral=float(scalar["cue_region_amount"][i]),
                        global_maximum_unwrapped_arc=float(peak_unwrapped[i]),first_harmonic_unwrapped_arc=float(moment_unwrapped[i]),
                        first_harmonic_magnitude=float(scalar["harmonic_magnitude"][i]),
                        maxima_over_native_prominence_reference=int(scalar["peaks_over_native_prominence_reference"][i]),
                        frozen_half_contrast_regions=int(scalar["frozen_basal_half_contrast_components"][i])))
                case["phases"][p]=dict(time=t,primary=states[:,0].copy(),scalars=scalar,path=path)
                previous=states[-1].copy()
                if p=="recovery":
                    tail=pilot.tail_diagnostic(t,states,scales)
                    assert tail==record["tail"]
            del states
        assert pool_error<=1e-7 and minimum>=-5e-10
        cue,rec=case["phases"]["cue"],case["phases"].get("recovery")
        baseline=float(cue["scalars"]["old_core_mean"][0])
        summary.append(dict(slot=job["slot"],model=job["model"],signed_distance=job["signed_distance"],role=job["role"],status=record["status"],
            original_core_width=core["width"],initial_old_core_mean=baseline,initial_cue_region_mean=float(cue["scalars"]["cue_region_mean"][0]),
            old_core_mean_at_withdrawal=float(cue["scalars"]["old_core_mean"][-1]),
            cue_mean_at_withdrawal=float(cue["scalars"]["cue_region_mean"][-1]),
            final_old_core_mean=float(rec["scalars"]["old_core_mean"][-1]) if rec else None,
            final_old_core_mean_over_basal=float(rec["scalars"]["old_core_mean"][-1]/baseline) if rec else None,
            final_cue_region_mean=float(rec["scalars"]["cue_region_mean"][-1]) if rec else None,
            final_maxima_over_native_prominence_reference=int(rec["scalars"]["peaks_over_native_prominence_reference"][-1]) if rec else None,
            max_pool_relative_error=pool_error,minimum_concentration=minimum))
        audits.append(dict(slot=job["slot"],checkpoint_identity=True,exact_phase_continuity=True,
            fractional_core_and_cue_integrals_recomputed=True,raw_peak_and_valley_values_match_native_fields=True,
            full_native_tail_metrics_recomputed=True,max_pool_relative_error=pool_error,minimum_concentration=minimum))
        data[job["slot"]]=case
    return manifest,grid,data,summary,traces,audits


def reflection_checks(manifest,grid):
    rows=[]
    for model,base in manifest["models"].items():
        positive=next(job for job in manifest["jobs"] if job["model"]==model and job["signed_distance"]==12)
        negative=next(job for job in manifest["jobs"] if job["model"]==model and job["signed_distance"]==-12)
        reflect=(2*base["incumbent_peak_cell"]-np.arange(grid.cells))%grid.cells
        scales=np.array([base["fixed_field_scales"][name] for name in base["state_names"]])
        for phase in ("cue","recovery"):
            with np.load(RUN/f'workers/{positive["slot"]:02d}/{phase}.npz') as a,np.load(RUN/f'workers/{negative["slot"]:02d}/{phase}.npz') as b:
                np.testing.assert_array_equal(a["time"],b["time"])
                av,bv=a["states"],b["states"]
                maximum=float(np.max(np.abs(av[:,:,reflect]-bv)/scales[None,:,None]))
                final=float(np.max(np.abs(av[-1][:,reflect]-bv[-1])/scales[:,None]))
                rows.append(dict(model=model,phase=phase,control="Reflected cue, same unmodified basal checkpoint",max_full_field_scaled_error=maximum,final_scaled_error=final,
                    initial_scaled_reflection_error=base["initial_full_field_reflection_error"]))
    return rows


def two_site_metrics(grid,values,distance):
    peaks,props=find_peaks(np.tile(values,3),prominence=.05)
    keep=(peaks>=800)&(peaks<1600);peaks=peaks[keep]-800
    if len(peaks)!=2:return None
    origin=grid.x[0];old=int(peaks[np.argmin(grid.distance(grid.x[peaks],origin))]);new=int(peaks[np.argmin(grid.distance(grid.x[peaks],origin+distance))])
    if old==new:return None
    path1=(old+np.arange((new-old)%800+1))%800;path2=(new+np.arange((old-new)%800+1))%800
    valley1=int(path1[np.argmin(values[path1])]);valley2=int(path2[np.argmin(values[path2])])
    floor=max(values[valley1],values[valley2])
    widths=[]
    for peak in (old,new):
        level=floor+.5*(values[peak]-floor);left=right=peak
        while values[(left-1)%800]>=level:left=(left-1)%800
        while values[(right+1)%800]>=level:right=(right+1)%800
        lfrac=(level-values[(left-1)%800])/(values[left]-values[(left-1)%800])
        rfrac=(values[right]-level)/(values[right]-values[(right+1)%800])
        widths.append(float(((right-left)%800+1-lfrac+rfrac)*grid.dx))
    q=(grid.x-origin)%40;lo=distance/2;hi=lo+20
    neww=np.maximum(0,np.minimum(q+grid.dx/2,hi)-np.maximum(q-grid.dx/2,lo));oldw=grid.dx-neww
    oldamount=float(values@oldw);newamount=float(values@neww)
    oldexcess=float((values-values.min())@oldw);newexcess=float((values-values.min())@neww)
    return dict(old_peak_arc=float(grid.arc(grid.x[old],origin)),cue_peak_arc=float(grid.arc(grid.x[new],origin)),
        peak_shortest_separation=float(grid.distance(grid.x[old],grid.x[new])),old_peak_height=float(values[old]),cue_peak_height=float(values[new]),
        old_half_local_prominence_width=widths[0],cue_half_local_prominence_width=widths[1],
        short_arc_valley_arc=float((grid.x[valley1]-origin)%40),short_arc_valley_height=float(values[valley1]),
        other_arc_valley_arc=float((grid.x[valley2]-origin)%40),other_arc_valley_height=float(values[valley2]),
        fixed_old_basin_active_amount=oldamount,fixed_cue_basin_active_amount=newamount,
        signed_cue_minus_old_active_amount=newamount-oldamount,
        signed_raw_amount_imbalance=(newamount-oldamount)/(newamount+oldamount),
        signed_above_minimum_amount_imbalance=(newexcess-oldexcess)/(newexcess+oldexcess))


def holmes_details(grid,data):
    rows=[];tails=[]
    for case in data.values():
        job=case["job"]
        if job["model"]!="holmes_model3" or job["signed_distance"] not in (16,20):continue
        for phase in ("cue","recovery"):
            z=case["phases"][phase]
            indices=[len(z["time"])-1] if phase=="cue" else np.flatnonzero(z["time"]>=12000)
            for i in indices:
                metric=two_site_metrics(grid,z["primary"][i],job["signed_distance"])
                if metric is not None:rows.append(dict(distance=job["signed_distance"],phase=phase,time_in_phase=float(z["time"][i]),**metric))
        tail=[row for row in rows if row["distance"]==job["signed_distance"] and row["phase"]=="recovery"]
        t=np.array([row["time_in_phase"] for row in tail]);y=np.array([row["signed_raw_amount_imbalance"] for row in tail])
        slope=float(np.dot(t-t.mean(),y-y.mean())/np.dot(t-t.mean(),t-t.mean()))
        raw=np.array([row["signed_cue_minus_old_active_amount"] for row in tail]);rawslope=float(np.dot(t-t.mean(),raw-raw.mean())/np.dot(t-t.mean(),t-t.mean()))
        tails.append(dict(distance=job["signed_distance"],start_time=float(t[0]),end_time=float(t[-1]),
            signed_raw_imbalance_at_start=float(y[0]),signed_raw_imbalance_at_end=float(y[-1]),
            signed_imbalance_change=float(y[-1]-y[0]),OLS_imbalance_drift_per_native_time=slope,
            signed_cue_minus_old_active_amount_at_start=float(raw[0]),signed_cue_minus_old_active_amount_at_end=float(raw[-1]),
            signed_cue_minus_old_active_amount_change=float(raw[-1]-raw[0]),OLS_signed_active_amount_drift_per_native_time=rawslope,
            observation="Finite two-site tail with slight amplitude/amount asymmetry; not a stability eigenvalue or proof against slower competition."))
    return rows,tails


def heatmaps(manifest,grid,data,early):
    fig,axes=plt.subplots(3,4,figsize=(8.3,6.8),sharex=True,sharey=True)
    fig.subplots_adjust(left=.10,right=.91,bottom=.085,top=.93,hspace=.32,wspace=.12)
    end=2000 if early else 14000
    for row,(model,label) in enumerate(NAMES.items()):
        cases=[case for case in data.values() if case["job"]["model"]==model and case["job"]["role"]=="distance case"]
        vmax=max(max(case["phases"]["cue"]["primary"].max(),case["phases"]["recovery"]["primary"].max()) for case in cases)
        norm=Normalize(0,vmax)
        for col,distance in enumerate((8,12,16,20)):
            ax=axes[row,col]
            case=next((case for case in cases if case["job"]["signed_distance"]==distance),None)
            if case is None:
                ax.text(.5,.5,"Excluded:\ncue overlaps\nincumbent core",transform=ax.transAxes,ha="center",va="center",fontsize=8,color="#555555")
                ax.set_facecolor("#f4f4f4")
            else:
                cue,rec=case["phases"]["cue"],case["phases"]["recovery"]
                t=np.r_[cue["time"],1000+rec["time"][1:]];u=np.concatenate([cue["primary"],rec["primary"][1:]])
                keep=t<=end;origin=case["core"]["origin"]+distance/2
                x=grid.arc(grid.x,origin);order=np.argsort(x)
                image=ax.pcolormesh(x[order],t[keep],u[keep][:,order],shading="auto",cmap="viridis",norm=norm,rasterized=True)
                ax.axvline(-distance/2,color="white",ls=":",lw=.65,alpha=.7)
                ax.axvline(distance/2,color="#F4C66A",ls=":",lw=.65,alpha=.9)
                ax.axhline(1000,color="white",ls="--",lw=.75)
            if row==0:ax.set_title(f"d = {distance}",fontsize=10)
            if col==0:ax.set_ylabel(label+"\nTime after cue onset")
            if row==2:ax.set_xlabel("Ring arc from midpoint")
            ax.set(xlim=(-20,20),ylim=(0,end),xticks=[-20,0,20],yticks=[0,1000,2000] if early else [0,1000,7000,14000])
        box=axes[row,-1].get_position();cax=fig.add_axes([.93,box.y0,.016,box.height])
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap="viridis"),cax=cax);cb.ax.tick_params(labelsize=7)
        cb.set_label("Native primary",fontsize=7)
    stem="distance_spatial_early" if early else "distance_spatial_full"
    paths=[]
    for ext in ("png","pdf","svg"):
        path=dev_output_path(OUT,stem+"."+ext);fig.savefig(path,dpi=240,facecolor="white");paths.append(path)
    plt.close(fig);return paths


def site_organization(grid,values):
    """Native local cap measurements, with valley-bounded spatial catchments."""
    tiled,props=find_peaks(np.tile(values,3),prominence=.05)
    peaks=np.sort(tiled[(tiled>=grid.cells)&(tiled<2*grid.cells)]-grid.cells)
    output=[];floor_global=float(values.min())
    for k,peak in enumerate(peaks):
        if len(peaks)==1:
            left=right=int(np.argmin(values));weights=np.full(grid.cells,grid.dx)
        else:
            previous=int(peaks[(k-1)%len(peaks)]);following=int(peaks[(k+1)%len(peaks)])
            a=(previous+np.arange((int(peak)-previous)%grid.cells+1))%grid.cells
            b=(int(peak)+np.arange((following-int(peak))%grid.cells+1))%grid.cells
            left=int(a[np.argmin(values[a])]);right=int(b[np.argmin(values[b])])
            catchment=(left+np.arange((right-left)%grid.cells+1))%grid.cells
            weights=np.zeros(grid.cells);weights[catchment]=grid.dx;weights[left]*=.5;weights[right]*=.5
        floor_local=float(max(values[left],values[right]));level=floor_local+.5*(values[peak]-floor_local)
        a=b=int(peak)
        while values[(a-1)%grid.cells]>=level:a=(a-1)%grid.cells
        while values[(b+1)%grid.cells]>=level:b=(b+1)%grid.cells
        lfrac=(level-values[(a-1)%grid.cells])/(values[a]-values[(a-1)%grid.cells])
        rfrac=(values[b]-level)/(values[b]-values[(b+1)%grid.cells])
        width=((b-a)%grid.cells+1-lfrac+rfrac)*grid.dx
        output.append(dict(peak_cell=int(peak),peak_arc=float(grid.arc(grid.x[peak],grid.x[0])),peak_height=float(values[peak]),
            local_base_height=floor_local,half_local_contrast_width=float(width),valley_basin_length=float(weights.sum()),
            valley_basin_raw_active_amount=float(values@weights),valley_basin_active_amount_above_profile_minimum=float((values-floor_global)@weights),
            valley_basin_active_amount_above_local_base=float(np.maximum(values-floor_local,0)@weights),
            whole_ring_raw_active_amount=float(values.sum()*grid.dx),whole_ring_minimum=floor_global))
    return output


def holmes_organization(grid,data):
    rows=[];activation=[]
    for distance in (12,16,20):
        case=next(case for case in data.values() if case["job"]["model"]=="holmes_model3" and case["job"]["signed_distance"]==distance)
        for phase,values in (("basal",case["original"][0]),("cue withdrawal",case["phases"]["cue"]["primary"][-1]),("13,000 s cue-free",case["phases"]["recovery"]["primary"][-1])):
            for metric in site_organization(grid,values):rows.append(dict(distance=distance,phase=phase,**metric))
        for time in (0,1,10,100,1000):
            z=case["phases"]["cue"];i=int(np.flatnonzero(z["time"]==time)[0])
            activation.append(dict(distance=distance,time_after_cue_onset=time,cue_region_mean=float(z["scalars"]["cue_region_mean"][i]),old_core_mean=float(z["scalars"]["old_core_mean"][i])))
    return rows,activation


def holmes_candidate(grid,data):
    """Four-panel candidate; full spatial evidence without a front gate."""
    fig,axes=plt.subplots(2,2,figsize=(7.1,5.5),sharex=True)
    fig.subplots_adjust(left=.10,right=.87,bottom=.095,top=.88,hspace=.34,wspace=.24)
    cases=[next(case for case in data.values() if case["job"]["model"]=="holmes_model3" and case["job"]["signed_distance"]==d) for d in (12,16)]
    vmax=max(case["phases"]["cue"]["primary"].max() for case in cases)
    norm=Normalize(0,vmax);rows=[]
    for col,case in enumerate(cases):
        distance=case["job"]["signed_distance"];cue=case["phases"]["cue"];rec=case["phases"]["recovery"]
        origin=case["core"]["origin"];x=(grid.x-origin+10)%40-10;order=np.argsort(x)
        t=np.r_[cue["time"],1000+rec["time"][1:]];u=np.concatenate([cue["primary"],rec["primary"][1:]])
        keep=t<=2000;ax=axes[0,col]
        ax.pcolormesh(x[order],t[keep],u[keep][:,order],shading="auto",cmap="viridis",norm=norm,rasterized=True)
        ax.axhline(1000,color="white",lw=.8,ls="--")
        ax.set_title(f'{chr(65+col)}   Holmes: d = {distance} µm',loc="left",fontsize=10,pad=10)
        ax.set(ylim=(0,2000),yticks=[0,1000,2000])
        if col==0:ax.set_ylabel("Time after cue onset (s)")
        ax=axes[1,col]
        profiles=[("Basal",case["original"][0],"#6e6e6e","--"),("Cue off",cue["primary"][-1],"#c17b0a","-"),("450 s after cue off",rec["primary"][np.flatnonzero(rec["time"]==450)[0]],"#9a65aa","-."),("13,000 s after cue off",rec["primary"][-1],"#087f8c","-")]
        for label,values,color,style in profiles:
            ax.plot(x[order],values[order],color=color,ls=style,lw=1.4,label=label)
            for j in order:rows.append(dict(distance=distance,phase=label,arc_in_common_old_at_zero_frame=float(x[j]),arc_from_original_incumbent=float(grid.arc(grid.x[j],case["core"]["origin"])),native_active_rac=float(values[j])))
        ax.set_title(f'{chr(67+col)}   '+("Valley fills into one broad cap" if distance==12 else "Two narrower caps at endpoint"),loc="left",fontsize=9,pad=9)
        ax.set(ylim=(0,2.5),yticks=[0,1,2]);ax.set_xlabel("Ring arc from original incumbent (µm)")
        if col==0:ax.set_ylabel("Active Rac (µM)")
        for row in range(2):
            ax=axes[row,col]
            ax.axvline(0,color="#e0e0e0" if row==0 else "#777777",ls=":",lw=.8)
            ax.axvline(distance,color="#F4C66A" if row==0 else "#b28015",ls=":",lw=.8)
            ax.set(xlim=(-10,30),xticks=[-10,0,10,20,30])
        axes[0,col].text(0,2030,"old",ha="center",va="bottom",fontsize=7,color="#666666")
        axes[0,col].text(distance,2030,"cue",ha="center",va="bottom",fontsize=7,color="#9b6a00")
    cax=fig.add_axes([.90,.548,.018,.332]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap="viridis"),cax=cax)
    cb.set_label("Active Rac (µM)",fontsize=8);cb.ax.tick_params(labelsize=7)
    handles,labels=axes[1,0].get_legend_handles_labels();fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.48,.99),ncol=4,frameon=False,fontsize=7.5)
    paths=[]
    for ext in ("png","pdf","svg"):
        path=dev_output_path(OUT,"holmes_distance_candidate."+ext);fig.savefig(path,dpi=260,facecolor="white");paths.append(path)
    plt.close(fig);paths.append(csvwrite("holmes_candidate_profiles.csv",rows));return paths


def main():
    source(Path(__file__));source(Path(__file__).with_name("chemical_ring_distance_pilot_20261017.py"));source(Path(__file__).with_name("length_dev_output_layout.py"))
    manifest,grid,data,summary,traces,audits=load_and_check()
    native_units=pilot.ring.native_model("holmes_model3").metadata()["units"]
    assert native_units["concentration"]=="μM" and native_units["time"]=="s" and native_units["length"]=="μm"
    controls=reflection_checks(manifest,grid);sites,tail=holmes_details(grid,data)
    organization,activation=holmes_organization(grid,data)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,"pdf.fonttype":42,"svg.fonttype":"none"})
    files=[csvwrite("case_summary.csv",summary),csvwrite("raw_region_and_supporting_position_traces.csv",traces),
        csvwrite("reflection_controls.csv",controls),csvwrite("holmes_two_site_metrics.csv",sites),csvwrite("holmes_two_site_tail_drift.csv",tail),
        csvwrite("holmes_site_organization.csv",organization),csvwrite("holmes_early_activation.csv",activation)]
    files+=heatmaps(manifest,grid,data,True);files+=heatmaps(manifest,grid,data,False)
    files+=holmes_candidate(grid,data)
    audit=dict(generated_at=datetime.now(timezone.utc).isoformat(),passed=True,cases=audits,
        holmes_native_units=native_units,reflected_controls=controls,model_scope=list(NAMES),distance_cases=11,numerical_controls=3,
        no_retuning=True,baseline_checkpoints_untouched=True,
        peak_tracks="Unwrapped global-maximum and first-harmonic positions are supporting summaries, not identities of moving sites. Maximum switches and weak harmonics can jump; unwrapping is performed separately within each phase.")
    files.append(jsonwrite("verification.json",audit))
    files.append(jsonwrite("readout_definitions.json",dict(
        old_core="Original periodic half-contrast interval with fractional finite-volume edge cells; fixed throughout.",
        site_amounts="Holmes two-site diagnostics use fixed equal-length Voronoi basins of the original incumbent and cue centres (20 units each). Raw active amounts and above-current-minimum amounts are exported; sign is cue minus original.",
        site_organization="Native cap height, half-local-contrast width and raw/above-minimum active amounts in valley-bounded catchments. Active Rac is not separately conserved. Valley boundaries can drift, so fixed equal-length basins remain the tail-imbalance diagnostic.",
        site_width="Each site's width at half its own peak height above the higher of the two adjacent valley minima; linear interpolation between native cell centres.",
        tail_drift="OLS slope and endpoint change of the fixed-basin signed active-amount difference and its normalized imbalance over recovery 12,000–13,000; not an eigenvalue.",
        heatmap="Full raw primary profiles; horizontal coordinate centred halfway between original incumbent and cue so both sites lie inside the display. White dotted line original incumbent, gold dotted line cue centre, dashed horizontal line withdrawal at 1000. Periodic edges are connected. Shared native color range within each model, separately scaled between models.",
        thresholds="Frozen basal half-contrast segmentation and 0.05 native prominence reference are descriptive; no front-establishment threshold is introduced.")))
    table="\n".join(f'| {NAMES[r["model"]]} | {r["signed_distance"]:+g} | {r["role"]} | {r["status"]} | {r["final_old_core_mean_over_basal"]:.4f} | {r["final_maxima_over_native_prominence_reference"]} |' for r in summary)
    tailtext="\n".join(f'- Holmes d={r["distance"]}: normalized signed basin imbalance at 13,000 is {r["signed_raw_imbalance_at_end"]:.8g}; change over the last 1000 is {r["signed_imbalance_change"]:.4g}; fitted drift {r["OLS_imbalance_drift_per_native_time"]:.4g} per native time.' for r in tail)
    md=f'''# First chemical-ring distance pilot: response modes before a distance law

Generated {datetime.now(timezone.utc).isoformat()}; producer finished. This is a reviewable result draft, not a selected main-text figure. No new duration optimization or parameter fitting was performed.

## Proposed message

At fixed circumference and material pools, moving an equal local cue changes how activity reorganizes. The observed distinction is between bridging into one region, loss at the original site with a remaining displaced region, and a two-site pattern that remains through the recorded cue-free window. These trajectories do not support a universal “farther is slower” summary.

![Holmes separation candidate](../../holmes_distance_candidate.png)

This four-panel candidate compares Holmes d = 12 and 16. Both columns use the same fixed 40-µm display frame, with the original incumbent at 0. The two routes between incumbent and cue are 12/28 µm in the left column and 16/24 µm in the right. Top panels show the first 1000 seconds with cue and first 1000 seconds after withdrawal; bottom panels show the common basal profile, cue withdrawal, 450 seconds after withdrawal and the endpoint after 13,000 cue-free seconds. The short view resolves the bridge that develops at d = 12. The endpoint profiles distinguish one broad region from two separated regions without replacing the fields by a binary threshold. The full retained history remains below. Dotted lines identify the original incumbent and cue centre; the dashed horizontal line marks withdrawal. All concentrations are native, unnormalized Rac in µM; time is seconds, as declared by the Holmes adapter. No biological size calibration is added. The two-site endpoint is an observation over this window, not a claim of stability.

The main candidate isolates a within-model comparison: only cue location changes. The following survey shows which parts of this contrast transfer to the other chemical models.

![Early native spatial histories](../../distance_spatial_early.png)

All eleven distance cases are shown; Holmes d=8 is excluded because the cue support intersects the original incumbent half-contrast core. Color is raw active concentration (Rac for Holmes), with one native color range per model. Time is model-native (seconds for Holmes), not biologically calibrated across models. White dotted vertical lines mark the original incumbent; gold marks the cue centre; the horizontal dashed line marks cue withdrawal at 1000. The display origin is halfway between the original incumbent and cue; periodic edges connect. This first figure retains the first 1000 cue-on units and first 1000 after withdrawal, so emergence, merger and movement remain visible.

## Why this design conveys the result

The circumference remains 40 and every position restarts its model's exact same untouched N800 qualified basal array. Only cue location changes within a model. This holds available material, the established incumbent and native kinetics fixed while exposing how two regions interact. Showing the full fields prevents original-core loss from being mistaken for takeover by a distinct challenger. Fixed-core and cue-region measurements remain supporting readouts, not an inferred causal partition.

The symmetric raised-cosine cue has total physical support 3 and fixed native peaks Wave-pinning 0.1, Goryachev 6.4 and Holmes 1. It stays on for 1000 native units, followed by 13,000 cue-free units. Sampling is 0.1 through 10, 1 through 100 and 5 thereafter in each phase. Native reactions, diffusivities and conserved groups are unchanged. Cue rates and time scales are not matched between models. Different cue locations can also encounter different basal tails/substrate concentrations; this experiment does not isolate diffusion travel time from those local chemical interactions.

## Observed modes

- **Wave-pinning:** d=8 rapidly develops one connected moving/merging activity region. Larger offsets show two local maxima/regions followed by loss near the original location and a remaining displaced region; this differs from merger. Fixed original-core activity remains 0.595 of basal at d=8 but only 0.165 at d=20; neither scalar alone identifies which pattern moved or replaced another.
- **Goryachev:** a separate cue-local maximum develops and original-core activity declines at all four distances. The native spatial fields, rather than the frozen half-contrast component count, reveal the initially weaker second maximum. The latter count can stay one while two unequal maxima coexist.
- **Holmes 3:** d=12 fills the intervening valley after cue withdrawal and ends as one region between the original and challenged locations. At d=16 and 20, two activity regions remain through 13,000 cue-free units and pass the two late full-field stationarity diagnostics. These are **finite observed two-site tails**, not established stable equilibria. Very slow competition remains possible, especially near symmetric configurations.

## Full retained spatial histories

![Full cue-on and 13,000-unit cue-free histories](../../distance_spatial_full.png)

The full view includes every retained time sample through 14,000 units after cue onset. The early view above resolves the transient that this long time axis compresses. Full native fields, every raw local maximum with its prominence and adjacent valley bases, and original-core/cue-region integrals are in the worker archives linked below.

| Model | Signed distance | Role | Execution/tail status | Final fixed-core mean / basal | Final maxima above descriptive 0.05 prominence |
|---|---:|---|---|---:|---:|
{table}

The three negative-distance cases are reflected-cue **numerical symmetry controls**, not biological replicates. Positive/negative d=12 full native fields are compared after the corresponding spatial reflection in [control checks](../data/reflection_controls.csv). Initial reflection asymmetry of the unmodified basal arrays is retained and reported.

## Holmes two-site asymmetry and slow-tail limitation

Heights, positions, local half-contrast widths and both valley minima at withdrawal and late recovery are exported in [site metrics](../data/holmes_two_site_metrics.csv). Signed site amounts use two fixed equal-length spatial basins, so movement of a nearly flat valley does not create an artificial imbalance drift. Both raw and above-profile-minimum amounts are retained.

The [basal and final site organization table](../data/holmes_site_organization.csv) compares the original single cap with all final caps. Its amounts integrate raw active Rac, and active Rac above the instantaneous profile minimum, over each valley-bounded spatial catchment. Active Rac is not conserved separately; the native Rac active-plus-inactive pool is conserved. The original single cap has half-local-contrast width 16.5379. At d = 16, final original/cue caps have widths 8.5226 and 7.9210 while peak heights stay near 1.4847. This is cap narrowing with similar peak concentration in this model and state, not a claim of preserved old polarity or a causal material partition.

The [early activation table](../data/holmes_early_activation.csv) shows that d = 12 and 16 have almost matched cue-region activation: means at 1 second are 0.6594 and 0.6577, and at 10 seconds are 1.6938 and 1.6888. Both challenged regions activate before their later spatial fates differ.

At d=16 the final maxima are near −0.10 and +16.00 native length units, hence not antipodal. Heights are approximately 1.484687 and 1.484464, so the pair is not exactly amplitude-symmetric either. Their small tail changes do not exclude a slow competitive mode.

{tailtext}

[Tail-drift values](../data/holmes_two_site_tail_drift.csv) describe the finite recorded interval; they are not an eigenvalue or a stability test. A future bounded antisymmetric conserved-pool perturbation could test this concern, but none was run here.

## What is supported and what remains open

The intervention changes the spatial response mode at the tested locations under one fixed cue protocol. It shows why a finite two-site pattern, bridging and original-site extinction cannot be reduced to a common activation threshold or original-core-loss time.

It does not establish a monotone distance law, minimum successful cue duration, front-establishment gate, global stability, migration persistence or a biological WGD calibration. A final displaced cap alone cannot distinguish continuous motion from competition followed by loss of the old pattern. The late-window gates can miss slower competition. Basal profiles were checked at N400/N800 previously; challenge mesh convergence is untested. No failed or unsettled execution is silently removed: every worker and the explicit geometric exclusion appears above.

## Raw sources and checks

- [Frozen pilot manifest]({os.path.relpath(RUN/'manifest.json',OUT/'supporting/metadata')}) and [launch record]({os.path.relpath(RUN/'launch.json',OUT/'supporting/metadata')}); SLURM array 21083222, 14 one-CPU/3-GB workers, no array throttle.
- [Worker archives]({os.path.relpath(RUN/'workers',OUT/'supporting/metadata')}) contain full `states[time,species,cell]`, native time, stimulus, fixed core weights, all raw peak/valley arrays and pool/readout histories. All cases use copied checkpoint bytes with the original qualifying hashes.
- [Case summary](../data/case_summary.csv), [raw regional and supporting position traces](../data/raw_region_and_supporting_position_traces.csv), [verification](verification.json), [readout definitions](readout_definitions.json).
- Unwrapped global-maximum and circular first-harmonic positions are supporting summaries only. They can jump when peaks exchange rank or the moment becomes weak; they are not tracked front identities. Unwrapping is phase-local, and the two phases are not joined into an inferred trajectory.
- Source/output hashes and reproduction command are recorded in [manifest](manifest.json). Conserved pools, exact starts/phase continuity, fractional-region integrals, native peak/valley values and late-field diagnostics were independently reconstructed from saved arrays. No new integration occurs in this renderer.
'''
    caption="Holmes ring distance candidate. On the same C = 40 µm ring and exact common N = 800 basal state, an equal support-3-µm cue (native peak 1 µM/s) is placed 12 or 16 µm from the incumbent, with complementary routes 28 or 24 µm. A–B: raw active Rac for the first 1000 s cue-on and 1000 s cue-free; the dashed horizontal line marks withdrawal. C–D: common basal, withdrawal, 450 s after withdrawal and endpoint after 13,000 cue-free seconds. Both columns use the same old-at-zero frame; dotted lines mark original and cue sites. The nearer pair fills its valley, whereas two narrower caps remain at the farther endpoint. The native active-plus-inactive Rac pool is conserved; active Rac alone is not. These finite observations do not establish branch stability, challenge mesh convergence or migration persistence. The full 14,000-s histories and other chemical models are supporting results."
    caption_path=dev_output_path(OUT,"caption.txt");caption_path.write_text(caption+"\n");files.append(caption_path)
    path=dev_output_path(OUT,"RESULTS.md");path.write_text(md);files.append(path)
    manifest_out=dict(generated_at=datetime.now(timezone.utc).isoformat(),producer="finished",artifact_root="../..",
        status="First scientific draft awaiting fresh-context critique",source_sha256=SOURCES,
        sha256={dev_artifact_key(OUT,path):sha(path) for path in files},renderer_sha256=sha(Path(__file__)),
        reproduction="PYTHONDONTWRITEBYTECODE=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python "+str(Path(__file__)))
    jsonwrite("manifest.json",manifest_out)
    print(json.dumps(dict(report=str(path),distance_cases=11,controls=3,source_hashes=len(SOURCES),holmes_tail=tail),indent=2))


if __name__=="__main__":main()
