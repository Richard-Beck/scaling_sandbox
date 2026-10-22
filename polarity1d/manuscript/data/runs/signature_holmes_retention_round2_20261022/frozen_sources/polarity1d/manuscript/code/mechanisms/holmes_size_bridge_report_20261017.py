"""Read-only three-size spatial comparison; no integration or tuning."""
from pathlib import Path
import csv
import json
import os
import sys
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths
import holmes_size_bridge_20261017 as pilot

ROOT, RUN = pilot.ROOT, pilot.RUN
OUT = ROOT / "figures/length_mechanisms_20261008/dev/holmes_size_bridge_20261017"
FIELDS = ["Active Rac", "Inactive Rac", "Active Rho", "Inactive Rho"]


def csv_write(path, rows):
    with path.open("w", newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def save(fig, name):
    for ext in ("png","pdf","svg"):
        fig.savefig(OUT/f"{name}.{ext}",dpi=180,bbox_inches="tight")
    plt.close(fig)


def site_metrics(grid, state, origin, circumference, phase):
    u=state[0];tiled=np.tile(u,3)
    loc,prop=find_peaks(tiled,prominence=.05)
    keep=(loc>=grid.cells)&(loc<2*grid.cells);loc=loc[keep]
    widths,levels,left,right=peak_widths(tiled,loc,rel_height=.5)
    rows=[]
    for i,p in enumerate(loc):
        rows.append(dict(circumference=circumference,phase=phase,peak_cell=int(p%grid.cells),
            arc_from_original_peak=float(grid.arc(grid.x[p%grid.cells],origin)),
            height=float(tiled[p]),prominence=float(prop["prominences"][keep][i]),
            half_prominence_width=float(widths[i]*grid.dx),
            left_valley=float(tiled[prop["left_bases"][keep][i]]),
            right_valley=float(tiled[prop["right_bases"][keep][i]])))
    return rows


def load_phase(path, grid, start, end, cue_center, cue_on):
    with np.load(path) as z:
        times=z["time"].copy();states=z["states"].copy();stimulus=z["stimulus"].copy()
    assert states.shape==(len(times),4,grid.cells)
    assert np.isfinite(states).all() and states.min()>=-5e-10
    np.testing.assert_array_equal(states[0],start)
    np.testing.assert_array_equal(times,pilot.distance.sampling(end))
    expected=pilot.ring.ring_cue(grid,cue_center,3,1) if cue_on else np.zeros(grid.cells)
    np.testing.assert_array_equal(stimulus,expected)
    pool=np.column_stack([states[:,group,:].sum(axis=(1,2))*grid.dx for group in ([0,1],[2,3])])
    np.testing.assert_allclose(pool[0]/grid.circumference,[7.5,3.1],rtol=0,atol=1e-12)
    error=float(np.max(np.abs(pool/pool[0]-1)));assert error<=1e-7
    return dict(time=times,states=states,pools=pool,pool_error=error,
        minimum=float(states.min()),cue_integral=float(stimulus.sum()*grid.dx))


def main():
    for name in ("supporting/data","supporting/metadata"):(OUT/name).mkdir(parents=True,exist_ok=True)
    manifest=pilot.read(RUN/"manifest.json")
    assert pilot.sha(RUN/"manifest.json")==pilot.read(RUN/"launch.json")["manifest_sha256"]
    for path,h in manifest["source_sha256"].items():assert pilot.sha(RUN/"frozen_sources"/path)==h,path
    for path,h in manifest["reference_input_sha256"].items():assert pilot.sha(RUN/path)==h,path
    cases=[];summary=[];sites=[];snapshots=[];traces=[];checks=[]
    for c,n,d,slot in ((30,1200,12,0),(40,1600,16,None),(50,2000,20,1)):
        grid=pilot.ring.RingGrid(c,n)
        directory=RUN/("inputs/reference_C40" if slot is None else f"workers/{slot:02d}")
        status=pilot.read(directory/"status.json")
        if status["status"]!="complete":
            summary.append(dict(circumference=c,cells=n,distance=d,execution_status=status["status"],
                basal_core_width="unavailable",rac_total="unavailable",rho_total="unavailable",
                terminal_resolved_peaks="unavailable",terminal_peak_positions="unavailable",
                terminal_half_prominence_widths="unavailable",late_windows_steady="unavailable"))
            continue
        basalpath=directory/("basal.npz" if slot is None else "qualified_basal.npz")
        with np.load(basalpath) as z:basal=z["state"].copy()
        if slot is not None:
            assert pilot.sha(basalpath)==status["checkpoint_sha256"]
            audit=pilot.read(RUN/f"preparation_C{c}/results/independent_audit.json");assert audit["passed"]
        core=pilot.distance.original_core(grid,basal[0]);center=(core["peak_cell"]+int(d/grid.dx))%n
        assert not np.any((core["weights"]>0)&(pilot.ring.ring_cue(grid,center,3,1)>0))
        original=basal.copy();phases={}
        for name,end,on in (("cue",1000,True),("recovery",13000,False)):
            path=directory/f"{name}.npz"
            if slot is not None:
                info=next(p for p in status["phases"] if p["phase"]==name)
                assert pilot.sha(path)==info["sha256"]
            phase=load_phase(path,grid,original,end,center,on);phases[name]=phase
            original=phase["states"][-1].copy()
            checks.append(dict(circumference=c,phase=name,maximum_pool_relative_error=phase["pool_error"],
                minimum_concentration=phase["minimum"],cue_integral=phase["cue_integral"],
                exact_initial_and_sampling=True,source_sha256=pilot.sha(path)))
            primary=phase["states"][:,0]
            old=primary@core["weights"];cuew=pilot.distance.weights_about(grid,grid.x[center],-1.5,1.5)
            cueregion=primary@cuew
            for i,t in enumerate(phase["time"]):
                traces.append(dict(circumference=c,phase=name,phase_time=t,
                    original_fixed_core_active_Rac=old[i],cue_fixed_region_active_Rac=cueregion[i],
                    whole_ring_active_Rac=float(primary[i].sum()*grid.dx),
                    Rac_pool=phase["pools"][i,0],Rho_pool=phase["pools"][i,1]))
        state_snapshots={"basal":basal,"withdrawal":phases["cue"]["states"][-1],"terminal":phases["recovery"]["states"][-1]}
        for name,state in state_snapshots.items():
            sites+=site_metrics(grid,state,core["origin"],c,name)
            arc=grid.arc(grid.x,core["origin"]+d/2)+d/2
            for cell in range(n):
                row=dict(circumference=c,cells=n,distance=d,phase=name,arc_from_original_peak=arc[cell])
                row.update({field:state[k,cell] for k,field in enumerate(FIELDS)});snapshots.append(row)
        endpoints=site_metrics(grid,state_snapshots["terminal"],core["origin"],c,"terminal")
        tail=pilot.distance.tail_diagnostic(phases["recovery"]["time"],phases["recovery"]["states"],pilot.SCALES)
        amount=pilot.ring.pool_amounts(pilot.ring.native_model(pilot.KEY),grid,basal)
        summary.append(dict(circumference=c,cells=n,distance=d,execution_status="complete",
            basal_core_width=core["width"],rac_total=amount[0],rho_total=amount[1],
            terminal_resolved_peaks=len(endpoints),
            terminal_peak_positions="; ".join(f"{p['arc_from_original_peak']:.3f}" for p in endpoints),
            terminal_half_prominence_widths="; ".join(f"{p['half_prominence_width']:.5f}" for p in endpoints),
            late_windows_steady=tail["both_late_windows_steady"]))
        cases.append(dict(c=c,n=n,d=d,grid=grid,core=core,phases=phases,snapshots=state_snapshots))
    assert len(summary)==3
    csv_write(OUT/"supporting/data/size_summary.csv",summary)
    csv_write(OUT/"supporting/data/peak_valley_widths.csv",sites)
    csv_write(OUT/"supporting/data/full_field_snapshots.csv",snapshots)
    csv_write(OUT/"supporting/data/spatial_region_traces.csv",traces)
    plt.rcParams.update({"font.size":10,"axes.titlesize":11})
    fig,axes=plt.subplots(2,3,figsize=(11.5,6.2),constrained_layout=True)
    fullfig,fullaxes=plt.subplots(1,3,figsize=(11.5,3.3),constrained_layout=True)
    for col,case in enumerate(cases):
        g=case["grid"];origin=case["core"]["origin"];d=case["d"]
        arc=g.arc(g.x,origin+d/2)+d/2;order=np.argsort(arc)
        cue=case["phases"]["cue"];rec=case["phases"]["recovery"]
        t=np.r_[cue["time"],1000+rec["time"][1:]]
        activity=np.concatenate([cue["states"][:,0],rec["states"][1:,0]],axis=0)
        early=t<=4000
        for ax,keep in ((axes[0,col],early),(fullaxes[col],np.ones(len(t),dtype=bool))):
            mesh=ax.pcolormesh(arc[order],t[keep],activity[keep][:,order],shading="auto",vmin=.35,vmax=2.35,cmap="viridis",rasterized=True)
            if ax is axes[0,col]:mainmesh=mesh
            else:fullmesh=mesh
            ax.axhline(1000,color="white",ls="--",lw=.8)
            ax.axvline(0,color="white",ls=":",lw=.7);ax.axvline(d,color="#f3c26a",ls=":",lw=.7)
            ax.set(title=f"C = {case['c']} µm, d = {d} µm",xlabel="Arc from incumbent (µm)",xlim=(-15,35))
        axes[0,col].set_ylabel("Time from cue onset (s)" if col==0 else "")
        fullaxes[col].set_ylabel("Time from cue onset (s)" if col==0 else "")
        for phase,color,ls,label in (("basal","#777777","--","Basal"),("withdrawal","#d8793c","-","Withdrawal"),("terminal","#244e78","-","13,000 s off")):
            axes[1,col].plot(arc[order],case["snapshots"][phase][0,order],color=color,ls=ls,lw=1.4,label=label)
        axes[1,col].axvline(0,color=".75",ls=":",lw=.7);axes[1,col].axvline(d,color=".75",ls=":",lw=.7)
        axes[1,col].set(xlim=(-15,35),ylim=(.3,2.4),xlabel="Arc from incumbent (µm)",ylabel="Active Rac (µM)" if col==0 else "")
        axes[1,col].spines[["top","right"]].set_visible(False)
    axes[1,0].legend(frameon=False,fontsize=8,loc="upper left")
    fig.colorbar(mainmesh,ax=list(axes[0,:]),label="Active Rac (µM)",shrink=.9)
    fullfig.colorbar(fullmesh,ax=list(fullaxes),label="Active Rac (µM)",shrink=.9)
    save(fig,"size_bridge_spatial");save(fullfig,"size_bridge_full_history")
    fieldsfig,fieldsaxes=plt.subplots(3,4,figsize=(12,7.5),constrained_layout=True)
    for row,case in enumerate(cases):
        g=case["grid"];arc=g.arc(g.x,case["core"]["origin"]+case["d"]/2)+case["d"]/2;order=np.argsort(arc)
        for k,field in enumerate(FIELDS):
            ax=fieldsaxes[row,k]
            for name,color,ls in (("basal","#777777","--"),("terminal","#244e78","-")):
                ax.plot(arc[order],case["snapshots"][name][k,order],color=color,ls=ls,lw=1.3,label=name)
            ax.set(title=f"C{case['c']}: {field}",xlim=(-15,35),xlabel="Arc (µm)",ylabel="µM")
            ax.spines[["top","right"]].set_visible(False)
    fieldsaxes[0,0].legend(frameon=False,fontsize=8)
    save(fieldsfig,"size_bridge_all_fields")
    maximum_pool=max(p["maximum_pool_relative_error"] for p in checks)
    pilot.write(OUT/"supporting/metadata/verification.json",dict(generated_at=pilot.now(),passed=True,
        frozen_sources_verified=len(manifest["source_sha256"]),frozen_reference_inputs_verified=len(manifest["reference_input_sha256"]),
        all_cases_retained=True,full_field_exact_starts_cues_samples_and_pools_verified=True,
        maximum_sampled_pool_relative_error=maximum_pool,phases=checks))
    table="\n".join(f"| {r['circumference']} | {r['distance']} | {r['cells']} | {r['basal_core_width']:.5f} | {r['rac_total']:.0f} / {r['rho_total']:.0f} | {r['terminal_resolved_peaks']} | {r['terminal_peak_positions']} | {r['terminal_half_prominence_widths']} |" for r in summary)
    text=f"""# A joint size/separation change alters the recorded Holmes spatial mode

Generated {pilot.now()}; producer finished. The two new preparations and two challenges are complete. The C40 reference is the existing fine-mesh result, reused without integration.

At the fixed fractional challenger separation d/C = 0.4, C30 loses the original incumbent region and ends with one Rac cap near the cue, at 11.95 µm, whereas C40 and C50 retain two separated caps through 13,000 cue-free seconds. The response therefore differs along this specified joint size/separation change. The C30 cue-aligned endpoint differs from the between-sites merger in the separate fixed-C40, d12 experiment. This pilot does not isolate a size mechanism: increasing circumference also increases total material, changes the freshly established incumbent width and moves the challenger farther away.

![Size bridge: spatial histories and raw profiles](../../size_bridge_spatial.png)

The top row retains the full first 1000 s with cue and first 3000 s after withdrawal, including the loss of the original C30 region; the dashed horizontal line marks withdrawal. Dotted lines mark the old incumbent and cue. Color is raw active Rac with one common scale. The lower row shows raw basal, withdrawal and terminal profiles. Common physical horizontal limits expose the different circumferences; blank space outside each unwrapped ring segment is outside that ring's displayed coordinate interval. Its left and right edges connect periodically. The bottom profiles distinguish the single cue-aligned broad cap at C30 from the two caps at C40/C50 without treating an original-core integral as survival.

| C (µm) | d (µm) | Cells | Basal half-contrast width (µm) | Total Rac / Rho (µM·µm) | Terminal peaks | Positions from incumbent (µm) | Half-prominence widths (µm) |
|---:|---:|---:|---:|---:|---:|---|---|
{table}

All three circumferences use the same physical cell width 0.025 µm, native pool densities 7.5 µM Rac and 3.1 µM Rho, unchanged native reaction/diffusion parameters, and the same local cue of total support 3 µm and peak 1 µM/s. Preparation starts afresh from the native homogeneous state: 10 s of the preparation cue at cell 0, then the unchanged full-field basal qualification and additional 10,000-s holdout. Both new sizes qualify as one-site basals. Each challenge starts its own exact qualified full-state checkpoint; no state is interpolated or reflected. The challenger is 12, 16 or 20 µm from the measured incumbent at C30, C40 or C50, respectively. Its support remains 3 µm, its peak 1 µM/s, and its duration 1000 s, followed by 13,000 s of native basal dynamics.

The new fixed-density pool amounts are 225/93 at C30 and 375/155 µM·µm at C50, compared with 300/124 at C40. Basal widths and the physical distances change together. Both periodic routes also grow: 12/18, 16/24 and 20/30 µm. The unchanged 3-µm cue occupies 10%, 7.5% and 6% of the circumference. This pilot is a direct observation along one size/separation path, not a causal allocation among material, incumbent width, both travel routes and local basal chemistry. The separate [fixed-circumference experiment](../../../holmes_distance_robustness_20261017/supporting/metadata/MESH_RESULTS.md) already establishes that cue location can alter the spatial mode at C40 while its incumbent and material remain fixed.

![Complete 14,000-s histories](../../size_bridge_full_history.png)

The complete histories retain every saved sample, exposing later changes compressed by the early view. Both new cases and the reference pass the two late full-field stationarity diagnostics, but that finite observation does not establish asymptotic stability. No failed preparation is converted into a biological outcome, and all three specified sizes appear above.

![All four native fields](../../size_bridge_all_fields.png)

This supporting comparison shows basal and terminal concentrations for active/inactive Rac and Rho. Conserved pools refer to active plus inactive forms; active Rac alone is not conserved. The maximum sampled conserved-pool relative error is {maximum_pool:.3g}. The integrated cue rate is 1.5 µM·µm/s at every size. BDF rtol 1e−9, atol 1e−11 and maximum step 10 s are unchanged; each phase samples every 0.1 s through 10, every 1 s through 100, then every 5 s. Field values are not clipped or mass-renormalized.

The result supports the limited statement that this prescribed increase in circumference and separation changes the observed one-cap/two-cap response. It does not identify a unique chemical cause, a transition circumference, a monotone distance law, a stable two-cap branch, or a migration-persistence mechanism. WGD-associated enlargement motivates the size question, but these periodic Rac patterns are distinct from observed migration summaries, conditional OU persistence estimates and interval reversal readouts.

## Reproduction and raw evidence

- [Frozen run manifest](../../../../../../data/runs/holmes_size_bridge_20261017/manifest.json), [launch](../../../../../../data/runs/holmes_size_bridge_20261017/launch.json), and [full worker archives](../../../../../../data/runs/holmes_size_bridge_20261017/workers). Array 21110502 has exactly two QoS-small, one-CPU/3-GB workers, no concurrency throttle. Each worker performs one preparation and one challenge; C40 is reused.
- [C30 preparation audit](../../../../../../data/runs/holmes_size_bridge_20261017/preparation_C30/results/independent_audit.json) and [C50 preparation audit](../../../../../../data/runs/holmes_size_bridge_20261017/preparation_C50/results/independent_audit.json) reconstruct the saved native homogeneous starts, full-field qualification, holdout, pool conservation and exact checkpoint continuity.
- [Size summary](../data/size_summary.csv), [all peak/valley/width values](../data/peak_valley_widths.csv), [four-field basal/withdrawal/end snapshots](../data/full_field_snapshots.csv), [fixed-region and native pool traces](../data/spatial_region_traces.csv), and [verification](verification.json). Exact full fields and every raw local maximum/prominence/valley are retained in the worker NPZ archives.
- Reproduce exports only with `{pilot.PYTHON} {Path(__file__)}`. Frozen integration sources and executable commands are retained in the run; original/native files and preceding result packages are unchanged.
"""
    if (OUT/"supporting/metadata/REVIEW_RESPONSE.md").exists():
        text += "\n[Fresh review and response](REVIEW_RESPONSE.md).\n"
    (OUT/"supporting/metadata/RESULTS.md").write_text(text)
    files=[p for p in OUT.rglob("*") if p.is_file() and p.name!="manifest.json"]
    pilot.write(OUT/"supporting/metadata/manifest.json",dict(generated_at=pilot.now(),producer="finished",
        run_manifest_sha256=pilot.sha(RUN/"manifest.json"),source_sha256={str(Path(__file__).relative_to(ROOT)):pilot.sha(Path(__file__))},
        artifacts={str(p.relative_to(OUT)):pilot.sha(p) for p in files}))
    print(json.dumps(dict(output=str(OUT),summary=summary,maximum_pool_relative_error=maximum_pool),indent=2))


if __name__=="__main__":main()
