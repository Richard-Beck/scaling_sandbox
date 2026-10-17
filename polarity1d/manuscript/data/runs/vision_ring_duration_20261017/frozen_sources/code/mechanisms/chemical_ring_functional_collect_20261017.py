"""Reconstruct all reported ring functional readouts from retained native arrays."""
from pathlib import Path
import argparse
import csv
import json
import os
import sys
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key]="1"
sys.dont_write_bytecode=True
import numpy as np
import chemical_ring_functional_v2_20261017 as campaign
import ring_distance_prototype_20261017 as ring
import chemical_ring_distance_pilot_20261017 as old

ROOT=Path(__file__).resolve().parents[2]


def csvwrite(path, rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=keys);writer.writeheader();writer.writerows(rows)


def native_state(path, time):
    with np.load(path) as z:
        if "states" in z.files:
            index=np.flatnonzero(np.isclose(z["time"],time,rtol=0,atol=1e-10))[0]
            return z["states"][index].copy()
        index=np.flatnonzero(np.isclose(z["snapshot_times"],time,rtol=0,atol=1e-10))[0]
        return z["snapshot_states"][index].copy()


def identity(primary,core,cue_center,stored):
    """Geographic association does not establish survival of a prior site."""
    result=dict(stored)
    components=ring.circular_components(primary>=core["level"])
    spanning=any(core["peak_cell"] in c and cue_center in c for c in components)
    result["single_region_spans_old_and_cue"] = bool(spanning and stored["resolved_peak_count"]==1)
    if result["single_region_spans_old_and_cue"]:
        result["location_association"]="one connected region spanning old and cue (merged or bridged)"
        result["distinct_old_peak_available"]=False
    elif stored["resolved_peak_count"]==1 and stored.get("old_core_peak_cells"):
        result["location_association"]="one old-core-associated peak; incumbent identity unresolved"
    return result


def main(run):
    m=campaign.check_sources(run)
    assert campaign.sha(run/"manifest.json")==campaign.read(run/"launch_expanded.json")["manifest_sha256"]
    for relative,digest in m["source_sha256"].items():
        if relative.startswith(("code/core/models/goryachev", "code/core/models/wave_pinning", "code/core/models/holmes_model3", "code/core/src/polarity1d/holmes", "code/mechanisms/ring_distance_prototype")):
            assert campaign.sha(ROOT/relative)==digest,"Live audit operator differs from frozen native implementation"
    matrix=[];basal=[];stages=[];sources={};probes=0
    maxerr={"pool":0.,"rhs_identity":0.,"probe_damage":0.,"probe_rate":0.,"gain_identity":0.}
    def source(path,digest=None):
        path=Path(path)
        if str(path) not in sources:
            measured=campaign.sha(path)
            if digest is not None:assert measured==digest,str(path)
            sources[str(path)]=measured
        elif digest is not None:
            assert sources[str(path)]==digest

    def check_probe(path,rows,model,grid,core,snapshot,stimulus):
        nonlocal probes
        source(path)
        with np.load(path) as z:
            np.testing.assert_array_equal(z["snapshot"],snapshot)
            np.testing.assert_array_equal(z["core_weights"],core["weights"])
            sham=z["sham_return"];shamrate=ring.native_rhs(model,grid,sham,stimulus)[0]
            sham_reaction=model.reaction(0.,grid.x,sham,stimulus)[0]
            sham_diffusion=model.diffusivities[0]*ring.periodic_laplacian(sham[0],grid.dx)
            primary=next(r for r in rows if r["epsilon"]==.005)
            tight_d=float((z["tight_sham_return"][0]-z["tight_probe_return"][0])@core["weights"])
            tight_n=float((ring.native_rhs(model,grid,z["tight_probe_return"],stimulus)[0]-ring.native_rhs(model,grid,z["tight_sham_return"],stimulus)[0])@core["weights"])
            np.testing.assert_allclose(abs(tight_d-primary["damage"]),primary["tighter_reference_damage_error"],rtol=1e-5,atol=1e-13)
            np.testing.assert_allclose(abs(tight_n-primary["net_rate_delta"]),primary["tighter_reference_net_rate_error"],rtol=1e-5,atol=1e-12)
            for row in rows:
                epsilon=row["epsilon"];returned=z[f"returned_{epsilon:g}"]
                row["sham_absolute_reaction_rate"]=float(sham_reaction@core["weights"])
                row["sham_absolute_diffusion_rate"]=float(sham_diffusion@core["weights"])
                row["probe_absolute_reaction_rate"]=float(model.reaction(0.,grid.x,returned,stimulus)[0]@core["weights"])
                row["probe_absolute_diffusion_rate"]=float(model.diffusivities[0]*ring.periodic_laplacian(returned[0],grid.dx)@core["weights"])
                poolerr=float(np.max(abs(ring.pool_amounts(model,grid,returned)/ring.pool_amounts(model,grid,snapshot)-1)))
                assert poolerr<1e-7;maxerr["pool"]=max(maxerr["pool"],poolerr)
                damage=float((sham[0]-returned[0])@core["weights"])
                net=float((ring.native_rhs(model,grid,returned,stimulus)[0]-shamrate)@core["weights"])
                maxerr["probe_damage"]=max(maxerr["probe_damage"],abs(damage-row["damage"]))
                maxerr["probe_rate"]=max(maxerr["probe_rate"],abs(net-row["net_rate_delta"]))
                np.testing.assert_allclose(damage,row["damage"],rtol=1e-8,atol=1e-12)
                np.testing.assert_allclose(net,row["net_rate_delta"],rtol=1e-7,atol=1e-11)
                if row["normalized_gain_available"]:
                    np.testing.assert_allclose(net/damage,row["net_repair_per_damage"],rtol=1e-7,atol=1e-9)
                probes+=1

    for job in m["jobs"]:
        statuspath=run/"workers"/f'{job["slot"]:02d}'/"status.json"
        source(statuspath);s=campaign.read(statuspath)
        assert s["complete"],str(statuspath)
        common=dict(model=job["model"],circumference=job["circumference"],cells=job["cells"],dx=job["dx"],
            worker_slot=job["slot"],worker_status=s["status"],status_path=str(statuspath))
        if s["status"]!="complete":
            for fraction in m["fractions"]:
                matrix.append(dict(**common,fraction=fraction,status=s["status"],preparation_status=s.get("preparation_status"),
                    requested_short_route=job["circumference"]*fraction))
            continue
        checkpoint=Path(s["checkpoint"]);source(checkpoint,s["checkpoint_sha256"])
        with np.load(checkpoint) as z:state=z["state"].copy()
        grid=ring.RingGrid(job["circumference"],job["cells"]);model=ring.native_model(job["model"])
        core=old.original_core(grid,state[0]);zero=np.zeros(grid.cells)
        common.update(checkpoint_sha256=s["checkpoint_sha256"],basal_width=core["width"],basal_peak=float(state[0].max()),
            basal_background=float(state[0].min()),basal_core_amount=float(state[0]@core["weights"]),
            basal_core_excess_amount=float((state[0]-state[0].min())@core["weights"]))
        check_probe(statuspath.parent/"basal_probe.npz",s["basal_probe"],model,grid,core,state,zero)
        for row in s["basal_probe"]:
            basal.append({**common,**row})
        for r in s["responses"]:
            rec=dict(**common,fraction=r["fraction"],status=r["status"],short_route=r["short_route"],long_route=r["long_route"],
                requested_short_route=r.get("requested_short_route",job["circumference"]*r["fraction"]),
                realized_fraction=r.get("realized_fraction",r["short_route"]/job["circumference"]),
                short_edge_to_cue_gap=r["short_edge_to_cue_gap"],long_edge_to_cue_gap=r["long_edge_to_cue_gap"],
                cue_center_cell=r["cue_center_cell"],overlaps_original_core=r["overlaps_original_core"])
            if r["status"]!="complete":matrix.append(rec);continue
            phases={p["phase"]:p for p in r["phases"]}
            for phase in phases.values():source(phase["path"],phase["sha256"])
            np.testing.assert_array_equal(native_state(phases["cue"]["path"],0.),state)
            np.testing.assert_array_equal(native_state(phases["cue"]["path"],1000.),native_state(phases["recovery"]["path"],0.))
            cue=ring.ring_cue(grid,r["cue_center_cell"],3.,campaign.previous.PEAKS[job["model"]])
            peak=campaign.previous.PEAKS[job["model"]];q=cue/peak
            derivative=ring.native_rhs(model,grid,state,cue)-ring.native_rhs(model,grid,state,zero)
            coefficient=state[1]/(7.5 if job["model"]=="holmes_model3" else 1.)
            err=float(np.max(abs(derivative[0]-cue*coefficient)));maxerr["rhs_identity"]=max(maxerr["rhs_identity"],err)
            assert err<1e-10
            gain=float(derivative[0].sum()/(peak*q.sum()))
            maxerr["gain_identity"]=max(maxerr["gain_identity"],abs(gain-r["gain"]["initial_cue_weighted_gain"]))
            np.testing.assert_allclose(gain,r["gain"]["initial_cue_weighted_gain"],rtol=1e-12,atol=1e-12)
            rec.update({k:v for k,v in r["gain"].items() if not isinstance(v,dict)})
            patch=old.weights_about(grid,float(grid.x[r["cue_center_cell"]]),-1.5,1.5)
            basal_reaction=model.reaction(0.,grid.x,state,zero)[0]
            basal_diffusion=model.diffusivities[0]*ring.periodic_laplacian(state[0],grid.dx)
            rec.update(basal_patch_mean_reaction=float(basal_reaction@patch/patch.sum()),
                basal_patch_mean_diffusion=float(basal_diffusion@patch/patch.sum()),
                basal_patch_mean_net=float((basal_reaction+basal_diffusion)@patch/patch.sum()))
            final_state=native_state(phases["recovery"]["path"],13000.)
            final_identity=identity(final_state[0],core,r["cue_center_cell"],r["final_spatial"])
            rec.update(cue_path=phases["cue"]["path"],recovery_path=phases["recovery"]["path"],
                reused_functional_status=r.get("reused_functional_status"),final_peak_count=r["final_spatial"]["resolved_peak_count"],
                final_location_association=final_identity["location_association"],
                final_single_region_spans_old_and_cue=final_identity["single_region_spans_old_and_cue"],
                final_old_core_amount_fraction=r["final_old_core_amount_fraction"])
            for name,field in zip(model.state_names,state):rec["cue_weighted_basal_"+name]=float((q*field).sum()/q.sum())
            for stage in r["stage_probes"]:
                t=stage["time_from_onset"];stimulus=cue if t<1000 else zero
                snapshot=native_state(phases["cue" if t<1000 else "recovery"]["path"],t if t<1000 else t-1000)
                source(stage["path"],stage["sha256"])
                check_probe(stage["path"],stage["rows"],model,grid,core,snapshot,stimulus)
                ident=identity(snapshot[0],core,r["cue_center_cell"],stage["spatial"])
                for row in stage["rows"]:
                    stages.append(dict(model=job["model"],circumference=job["circumference"],fraction=r["fraction"],
                        short_route=r["short_route"],long_route=r["long_route"],time_from_onset=t,
                        return_time_from_onset=t+.1,cue_on=stage["cue_on"],
                        resolved_peak_count=ident["resolved_peak_count"],location_association=ident["location_association"],
                        distinct_old_peak_available=ident["distinct_old_peak_available"],
                        single_region_spans_old_and_cue=ident["single_region_spans_old_and_cue"],
                        snapshot_absolute_net_rate=stage["snapshot_absolute_net_rate"],probe_path=stage["path"],**row))
                    stages[-1].update(snapshot_absolute_reaction_rate=float(model.reaction(0.,grid.x,snapshot,stimulus)[0]@core["weights"]),
                        snapshot_absolute_diffusion_rate=float(model.diffusivities[0]*ring.periodic_laplacian(snapshot[0],grid.dx)@core["weights"]))
                primary=next(row for row in stage["rows"] if row["epsilon"]==.005)
                for key in ("net_repair_per_damage","raw_ratio_estimate","normalized_gain_available","sham_absolute_net_rate","damage","numerical_gain_error_bound"):
                    rec[f"stage{t:g}_"+key]=primary[key]
                rec[f"stage{t:g}_distinct_old_peak_available"]=ident["distinct_old_peak_available"]
                rec[f"stage{t:g}_location_association"]=ident["location_association"]
            matrix.append(rec)
    assert len(matrix)==len(m["jobs"])*len(m["fractions"])==867
    assert len({(r["model"],r["circumference"],r["fraction"]) for r in matrix})==867
    anchor=next(r for r in matrix if r["model"]=="holmes_model3" and r["circumference"]==40 and r["fraction"]==.3)
    assert anchor["final_peak_count"]==1 and anchor["final_single_region_spans_old_and_cue"]
    out=run/"results"
    for name,rows in (("matrix",matrix),("basal",basal),("stage_probes",stages)):csvwrite(out/f"{name}.csv",rows)
    summary=dict(generated_at=campaign.now(),passed=True,planned_matrix_cells=867,model_size_workers=51,
        counts={name:sum(r["status"]==name for r in matrix) for name in sorted({r["status"] for r in matrix})},
        baseline_rows=len(basal),stage_rows=len(stages),probe_amplitudes_reconstructed=probes,
        exact_checkpoint_and_phase_continuity=True,all_reported_probes_reconstructed=True,
        full_867_cell_schema_coverage=True,max_absolute_errors=maxerr,
        scientific_limits="Algebraic/source/array audit; not continuum convergence, global stability, migration persistence or causal size decomposition.")
    campaign.write(out/"audit.json",summary)
    campaign.write(out/"source_sha256.json",sources)
    campaign.write(run/"status.json",dict(generated_at=campaign.now(),complete=True,producer="finished",audit=summary))
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--run",type=Path,default=campaign.DEFAULT)
    main(parser.parse_args().run.resolve())
