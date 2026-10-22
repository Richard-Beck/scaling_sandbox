"""Prospective three-model size/location comparison with conditional repair probes.

All integration runs in frozen SLURM workers. Existing C40 trajectories are read
by exact source hash, never remeshed. Native models and reviewed figures unchanged.
"""
from pathlib import Path
import argparse
import copy
import csv
import json
import os
import shlex
import shutil
import signal
import subprocess
import sys
import time
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np
from scipy.integrate import BDF, solve_ivp
from scipy.signal import find_peaks
import ring_distance_prototype_20261017 as ring
import ring_preparation_pilot_20261017 as prep
import chemical_ring_distance_pilot_20261017 as previous

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "data/runs/chemical_ring_functional_20261017"
TOPIC = ROOT / "figures/length_mechanisms_20261008/dev/chemical_ring_functional_20261017"
PYTHON = prep.PYTHON
KEYS = ring.KEYS
SIZES = (20, 30, 40, 50)
FRACTIONS = (.3, .4, .5)
EPSILONS = (.001, .005, .01)
STAGES = (10., 1000., 2000.)
read, write, sha, now = prep.read, prep.write, prep.sha, prep.now
STOP = False


def freeze(run):
    paths = [Path(__file__)] + [Path(__file__).with_name(name) for name in (
        "ring_distance_prototype_20261017.py", "ring_preparation_pilot_20261017.py",
        "chemical_ring_distance_pilot_20261017.py", "length_dev_output_layout.py")]
    for folder in (ROOT / "code/core/models", ROOT / "code/core/src/polarity1d"):
        paths += sorted(folder.glob("*.py"))
    hashes = {}
    for source in paths:
        rel = source.relative_to(ROOT); dest = run / "frozen_sources" / rel
        dest.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(source, dest)
        hashes[str(rel)] = sha(dest)
    return hashes


def initialize(run):
    if run.exists():
        raise ValueError("Use a fresh run directory")
    for name in ("inputs", "workers", "logs", "results", "preparations"):
        (run / name).mkdir(parents=True)
    sources = freeze(run)
    oldroot = ROOT / "data/runs/chemical_ring_distance_pilot_20261017"
    oldmanifest = read(oldroot / "manifest.json")
    ptemplate = read(ROOT / "data/runs/ring_preparation_pilot_N800_20261017/manifest.json")
    assert read(ROOT / "data/runs/ring_preparation_pilot_N800_20261017/results/independent_audit.json")["passed"]
    jobs, external = [], {}
    for key in KEYS:
        for circumference in SIZES:
            cells = int(circumference / .05)
            entry = dict(slot=len(jobs), model=key, circumference=circumference, cells=cells,
                dx=.05, fixed_scales=prep.SCALES[key], fractions=list(FRACTIONS))
            if circumference == 40:
                base = oldmanifest["models"][key]
                original = oldroot / base["checkpoint"]
                assert sha(original) == base["checkpoint_sha256"]
                dest = run / "inputs" / f"{key}_C40.npz"
                shutil.copy2(original, dest)
                entry.update(checkpoint=str(dest.relative_to(run)), checkpoint_sha256=sha(dest),
                    preparation="existing_exact_N800", original_checkpoint=str(original))
                external[str(original)] = sha(original)
                reuse = {}
                for fraction in FRACTIONS:
                    d = 40 * fraction
                    oldjob = next(j for j in oldmanifest["jobs"] if j["model"] == key and j["signed_distance"] == d)
                    directory = oldroot / "workers" / f'{oldjob["slot"]:02d}'
                    status = read(directory / "status.json")
                    assert status["complete"] and status["status"].startswith("complete_")
                    paths = {}
                    for phase in ("cue", "recovery"):
                        path = directory / f"{phase}.npz"
                        expected = next(p["sha256"] for p in status["phases"] if p["phase"] == phase)
                        assert sha(path) == expected
                        paths[phase] = str(path); external[str(path)] = expected
                    reuse[str(fraction)] = paths
                entry["reuse"] = reuse
            else:
                directory = run / "preparations" / f"{key}_C{circumference}"
                for name in ("inputs", "workers", "logs", "results"):
                    (directory / name).mkdir(parents=True)
                (directory / "frozen_sources").symlink_to(run / "frozen_sources", target_is_directory=True)
                grid = ring.RingGrid(circumference, cells)
                state = ring.native_rest(key, grid)
                path = directory / "inputs/homogeneous.npz"
                np.savez_compressed(path, state=state, x=grid.x)
                base = copy.deepcopy(next(j for j in ptemplate["jobs"] if j["model"] == key))
                base.update(slot=0, input=str(path.relative_to(directory)), input_sha256=sha(path),
                    native_pool_means=(ring.pool_amounts(ring.native_model(key), grid, state) / circumference).tolist())
                manifest = copy.deepcopy(ptemplate)
                manifest.update(created_at=now(), jobs=[base], circumference=circumference, cells=cells,
                    preparation_durations=[10.], source_sha256=sources,
                    integration_scope="One fixed 10-unit preparation, no alternate duration or parameter tuning")
                write(directory / "manifest.json", manifest)
                entry["preparation"] = str(directory.relative_to(run))
            jobs.append(entry)
    manifest = dict(created_at=now(), experiment="Three-model size/location native-ring functional characterization",
        jobs=jobs, sizes=list(SIZES), fractions=list(FRACTIONS), h=.05,
        cue=dict(width=3., peaks=previous.PEAKS, on=1000., off=13000.),
        probe=dict(epsilons=list(EPSILONS), primary_epsilon=.005, hold=.1,
            stages=list(STAGES), rtol=1e-11, atol=1e-13, max_step=.01,
            tight_reference="Every basal and challenged-stage primary epsilon=.005 pair; tolerances tenfold tighter",
            damage_gate="D>0 and D>100 times absolute primary-core difference from tighter return-state reference, also D>1e-11 times basal core amount",
            numerator_gate="Net-rate difference magnitude >100 times tighter-reference error for resolved signed normalized gain; raw estimate and error retained otherwise",
            zero_removal="Independent epsilon-zero solve must reproduce sham exactly; zero-duration algebraic withdrawal/return identity checked"),
        solver=dict(method="BDF", rtol=1e-9, atol=1e-11, max_step=10.),
        source_sha256=sources, external_sha256=external,
        frozen_plan_sha256=sha(TOPIC / "PLAN.md"),
        output_sampling="0.1 through10, 1 through100, 5 through1000, 50 after withdrawal; exact stages retained",
        interpretation="Conditional native inactive-withdrawal/return response, never generic stability; basal is measured once per size and has no future-distance dependence.")
    write(run / "manifest.json", manifest)
    write(run / "status.json", dict(generated_at=now(), complete=False, producer="not_started"))
    return manifest


def submit(run, slots, label):
    manifest = read(run / "manifest.json")
    runner = run / "frozen_sources" / Path(__file__).relative_to(ROOT)
    command = "exec " + shlex.join(["env", "PYTHONDONTWRITEBYTECODE=1", "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", PYTHON, str(runner), "worker", "--run", str(run)])
    command += ' --slot "$SLURM_ARRAY_TASK_ID"'
    args = ["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=3G",
        "--time=01:00:00", "--signal=B:USR1@60", "--array=" + slots, "--job-name=ring-functional",
        "--output=" + str(run / "logs/%A_%a.out"), "--error=" + str(run / "logs/%A_%a.err"), "--wrap", command]
    job = subprocess.check_output(args, text=True).strip().split(";")[0]
    record = dict(generated_at=now(), job_id=job, slots=slots, command=args,
        manifest_sha256=sha(run / "manifest.json"))
    write(run / f"launch_{label}.json", record)
    for base in manifest["jobs"]:
        if base["preparation"] != "existing_exact_N800":
            p = run / base["preparation"]
            write(p / "launch.json", dict(generated_at=now(), job_id=job,
                manifest_sha256=sha(p / "manifest.json"), outer_record=f"launch_{label}.json"))
    print(json.dumps(record, indent=2))


def check_sources(run):
    manifest = read(run / "manifest.json")
    for path, digest in manifest["source_sha256"].items():
        assert sha(run / "frozen_sources" / path) == digest, path
    return manifest


def integrate(model, grid, initial, cue, end, samples, pools, deadline):
    species = len(model.state_names)
    solver = BDF(lambda t,y:ring.native_rhs(model,grid,y.reshape(species,grid.cells),cue).ravel(),
        0., initial.ravel(), end, rtol=1e-9, atol=1e-11, max_step=10.,
        jac_sparsity=ring.periodic_sparsity(species,grid.cells))
    times, states, index = [0.], [initial.copy()], 1
    minimum, worst_pool = float(initial.min()), 0.
    while solver.status == "running":
        if STOP or time.monotonic() >= deadline:
            raise TimeoutError("Worker soft deadline/SLURM signal")
        message = solver.step()
        if solver.status == "failed":
            raise RuntimeError(message)
        accepted = solver.y.reshape(species, grid.cells)
        minimum = min(minimum, float(accepted.min()))
        worst_pool = max(worst_pool, float(np.max(np.abs(ring.pool_amounts(model,grid,accepted)/pools-1))))
        if minimum < -5e-10 or worst_pool > 1e-7:
            raise FloatingPointError("Concentration or conserved pool error")
        dense = solver.dense_output()
        while index < len(samples) and samples[index] <= solver.t:
            times.append(float(samples[index])); states.append(dense(samples[index]).reshape(species,grid.cells).copy()); index += 1
    return np.asarray(times), np.asarray(states), dict(minimum=minimum, max_pool_error=worst_pool,
        nfev=solver.nfev, njev=solver.njev, nlu=solver.nlu)


def short(model, grid, state, cue, duration=.1, tight=False):
    solution = solve_ivp(lambda t,y:ring.native_rhs(model,grid,y.reshape(state.shape),cue).ravel(),
        (0.,duration),state.ravel(),method="BDF",rtol=1e-12 if tight else 1e-11,
        atol=1e-14 if tight else 1e-13,max_step=.01,
        jac_sparsity=ring.periodic_sparsity(len(model.state_names),grid.cells),t_eval=[duration])
    if not solution.success:
        raise RuntimeError(solution.message)
    result = solution.y[:,-1].reshape(state.shape)
    assert result.min() >= -5e-10 and np.isfinite(result).all()
    np.testing.assert_allclose(ring.pool_amounts(model,grid,result),ring.pool_amounts(model,grid,state),rtol=1e-8,atol=1e-11)
    return result


def spatial(grid, primary, core, cue_center):
    peaks, prop = find_peaks(np.tile(primary,3),prominence=0.)
    keep = (peaks >= grid.cells) & (peaks < 2*grid.cells)
    peaks = peaks[keep] % grid.cells; prominence = prop["prominences"][keep]
    resolved = peaks[prominence >= .05]
    cue_origin = float(grid.x[cue_center])
    records = [dict(cell=int(p),position=float(grid.x[p]),value=float(primary[p]),
        prominence=float(v),distance_old=float(grid.distance(grid.x[p],core["origin"])),
        distance_cue=float(grid.distance(grid.x[p],cue_origin))) for p,v in zip(peaks,prominence)]
    old = [int(p) for p in resolved if core["weights"][p]>0 and
        grid.distance(grid.x[p],core["origin"]) < grid.distance(grid.x[p],cue_origin)]
    challenged = [int(p) for p in resolved if grid.distance(grid.x[p],cue_origin)<=1.5 and
        grid.distance(grid.x[p],cue_origin) < grid.distance(grid.x[p],core["origin"])]
    if old and challenged:
        identity = "distinct peaks associated with old core and cue patch"
    elif len(resolved)==1 and old:
        identity = "one peak associated with old core"
    elif len(resolved)==1 and challenged:
        identity = "one peak associated with cue patch"
    elif len(resolved)==1:
        identity = "one shifted or merged peak; original identity unresolved"
    elif len(resolved)==0:
        identity = "no peak exceeds native prominence reference"
    else:
        identity = "multiple peaks with ambiguous location associations"
    return dict(resolved_peak_count=len(resolved), all_peaks=records,
        location_association=identity,old_core_peak_cells=old,cue_patch_peak_cells=challenged,
        distinct_old_peak_available=bool(old and len(resolved)>1),
        frozen_half_contrast_components=len(ring.circular_components(primary >= core["level"])),
        minimum=float(primary.min()), maximum=float(primary.max()))


def probe(model, grid, snapshot, cue, core, reference, do_tight=False):
    w = core["weights"]
    sham = short(model,grid,snapshot,cue)
    sham_reaction = model.reaction(0.,grid.x,sham,cue)[0]
    sham_diffusion = model.diffusivities[0]*ring.periodic_laplacian(sham[0],grid.dx)
    base_amount = float(reference[0] @ w)
    rows, arrays = [], {"snapshot":snapshot,"sham_return":sham,"core_weights":w}
    zero_removal = short(model,grid,snapshot.copy(),cue)
    np.testing.assert_array_equal(zero_removal,sham)
    tight_sham = short(model,grid,snapshot,cue,tight=True)
    initial_tight = snapshot.copy(); initial_tight[1] *= .995
    tight_return = short(model,grid,initial_tight,cue,tight=True)
    tight_return[1] += .005*snapshot[1]
    arrays["tight_sham_return"] = tight_sham; arrays["tight_probe_return"] = tight_return
    for epsilon in EPSILONS:
        removed = epsilon*snapshot[1]
        initial = snapshot.copy(); initial[1] -= removed
        held = short(model,grid,initial,cue)
        returned = held.copy(); returned[1] += removed
        np.testing.assert_allclose(ring.pool_amounts(model,grid,returned),ring.pool_amounts(model,grid,snapshot),rtol=1e-8,atol=1e-11)
        reaction = model.reaction(0.,grid.x,returned,cue)[0]
        diffusion = model.diffusivities[0]*ring.periodic_laplacian(returned[0],grid.dx)
        damage = float((sham[0]-returned[0])@w)
        reaction_delta = float((reaction-sham_reaction)@w)
        diffusion_delta = float((diffusion-sham_diffusion)@w)
        net_delta = reaction_delta+diffusion_delta
        numerical_core_error, numerical_rate_error = None, None
        if epsilon == .005:
            numerical_core_error = float(abs(((tight_sham[0]-tight_return[0])-(sham[0]-returned[0]))@w))
            tight_delta = float((ring.native_rhs(model,grid,tight_return,cue)[0]-ring.native_rhs(model,grid,tight_sham,cue)[0])@w)
            numerical_rate_error = abs(tight_delta-net_delta)
        floor = max(1e-11*base_amount, 0. if numerical_core_error is None else 100*numerical_core_error)
        usable = bool(damage > floor)
        rows.append(dict(epsilon=epsilon,damage=damage,damage_fraction_basal_core=damage/base_amount,
            damage_fraction_current_sham=damage/float(sham[0]@w),basal_core_amount=base_amount,
            current_sham_core_amount=float(sham[0]@w),probe_core_amount=float(returned[0]@w),
            reaction_rate_delta=reaction_delta,diffusion_rate_delta=diffusion_delta,net_rate_delta=net_delta,
            sham_absolute_net_rate=float((sham_reaction+sham_diffusion)@w),
            probe_absolute_net_rate=float((reaction+diffusion)@w),
            net_repair_per_damage=net_delta/damage if usable else None,
            normalized_gain_available=usable,damage_floor=floor,
            tighter_reference_damage_error=numerical_core_error,tighter_reference_net_rate_error=numerical_rate_error,
            unavailable_reason=None if usable else "Nonpositive damage or damage below numerical gate",
            removed_amount=float(removed.sum()*grid.dx),
            returned_pool_error=float(np.max(np.abs(ring.pool_amounts(model,grid,returned)/ring.pool_amounts(model,grid,snapshot)-1)))))
        arrays[f"returned_{epsilon:g}"] = returned
    primary = next(row for row in rows if row["epsilon"]==.005)
    for row in rows:
        err_d = primary["tighter_reference_damage_error"]
        err_n = primary["tighter_reference_net_rate_error"]
        floor = max(1e-11*base_amount,100*err_d)
        usable_d = row["damage"]>floor
        numerator_resolved = abs(row["net_rate_delta"])>100*err_n
        raw = row["net_rate_delta"]/row["damage"] if row["damage"] != 0 else None
        row.update(damage_floor=floor,normalized_gain_available=bool(usable_d and numerator_resolved),
            numerator_sign_resolved=bool(numerator_resolved),raw_ratio_estimate=raw,
            net_repair_per_damage=raw if usable_d and numerator_resolved else None,
            numerical_gain_error_bound=(err_n+abs(raw)*err_d)/(row["damage"]-err_d) if usable_d else None,
            unavailable_reason=None if usable_d and numerator_resolved else "Damage or net-rate numerator below numerical reliability gate",
            zero_removal_identity=True)
    return rows, arrays


def cue_gain(model, grid, state, cue, core, center, peak, observed, no_cue):
    zero = np.zeros(grid.cells); q = cue/peak
    delta = ring.native_rhs(model,grid,state,cue)-ring.native_rhs(model,grid,state,zero)
    normalization = float(q.sum()*grid.dx)
    substrate = state[1]/(model.rac_total if model.name == "holmes_model3" else 1.)
    np.testing.assert_allclose(delta[0],cue*substrate,rtol=1e-11,atol=1e-12)
    np.testing.assert_allclose(delta[1],-delta[0],rtol=1e-11,atol=1e-12)
    if state.shape[0] > 2:
        np.testing.assert_allclose(delta[2:],0.,rtol=0,atol=1e-12)
    w = previous.weights_about(grid,float(grid.x[center]),-1.5,1.5)
    result = dict(initial_cue_weighted_gain=float(delta[0].sum()*grid.dx/(peak*normalization)),
        cue_weighted_inactive=float((q*state[1]).sum()*grid.dx/normalization),
        cue_total_rate=float(cue.sum()*grid.dx),initial_extra_active_amount_rate=float(delta[0].sum()*grid.dx),
        initial_extra_patch_mean_rate=float(delta[0]@w/w.sum()),
        cue_weighted_basal_active=float((q*state[0]).sum()*grid.dx/normalization),
        rhs_identity_max_error=float(np.max(np.abs(delta[0]-cue*substrate))))
    for t in (1.,10.):
        result[f"observed_patch_mean_extra_at_{t:g}"] = float((observed[t][0]-no_cue[t][0])@w/w.sum())
        result[f"observed_patch_mean_extra_per_peak_at_{t:g}"] = result[f"observed_patch_mean_extra_at_{t:g}"]/peak
    return result


def worker(run, slot):
    assert os.environ.get("SLURM_JOB_ID"), "All integration requires SLURM"
    signal.signal(signal.SIGUSR1,lambda s,f:globals().__setitem__("STOP",True))
    signal.signal(signal.SIGTERM,lambda s,f:globals().__setitem__("STOP",True))
    manifest = check_sources(run); job = manifest["jobs"][slot]
    out = run / "workers" / f"{slot:02d}"; out.mkdir()
    status = dict(generated_at=now(),producer="running",complete=False,job=job,status="preparation",responses=[])
    write(out / "status.json",status)
    began = time.monotonic()
    try:
        if job["preparation"] == "existing_exact_N800":
            checkpoint = run / job["checkpoint"]
            assert sha(checkpoint) == job["checkpoint_sha256"]
        else:
            p = run / job["preparation"]
            prep.worker(p,0)
            pstatus = read(p / "workers/00/status.json")
            status["preparation_status"] = pstatus["status"]
            if pstatus["status"] != "qualifying_single_site":
                status.update(complete=True,producer="finished",status="unqualified_preparation")
                write(out / "status.json",status); return
            assert prep.audit(p)["passed"]
            checkpoint = p / pstatus["checkpoint"]
        with np.load(checkpoint) as z:
            state = z["state"].copy()
        grid, model = ring.RingGrid(job["circumference"],job["cells"]),ring.native_model(job["model"])
        core = previous.original_core(grid,state[0]); zero = np.zeros(grid.cells)
        pools = ring.pool_amounts(model,grid,state)
        status.update(checkpoint=str(checkpoint),checkpoint_sha256=sha(checkpoint),
            core={k:v for k,v in core.items() if k != "weights"},status="basal_probe")
        baseline, arrays = probe(model,grid,state,zero,core,state,do_tight=job["circumference"] in (20,50))
        np.savez_compressed(out / "basal_probe.npz",**arrays)
        status["basal_probe"] = baseline
        t0, y0, _ = integrate(model,grid,state,zero,10.,np.array([0.,1.,10.]),pools,began+3300)
        no_cue = {float(t):y for t,y in zip(t0,y0)}
        for fraction in FRACTIONS:
            response = dict(fraction=fraction,short_route=job["circumference"]*fraction,
                long_route=job["circumference"]*(1-fraction))
            status["responses"].append(response)
            center = (core["peak_cell"]+int(round(fraction*grid.cells)))%grid.cells
            peak = previous.PEAKS[job["model"]]
            cue = ring.ring_cue(grid,center,3.,peak)
            overlap = bool(np.any((core["weights"]>0)&(cue>0)))
            d = response["short_route"]
            response.update(cue_center_cell=int(center),overlaps_original_core=overlap,
                short_edge_to_cue_gap=float(d-1.5-core["right"]),
                long_edge_to_cue_gap=float(job["circumference"]-d-1.5+core["left"]))
            if overlap:
                response.update(status="geometric_overlap_exclusion"); write(out / "status.json",status); continue
            response_dir = out / f"f{fraction:g}"; response_dir.mkdir()
            phase_data, phases = {}, []
            if "reuse" in job:
                for name,path_text in job["reuse"][str(fraction)].items():
                    path = Path(path_text)
                    assert sha(path) == manifest["external_sha256"][path_text]
                    with np.load(path) as z:
                        phase_data[name] = (z["time"].copy(),z["states"].copy())
                    phases.append(dict(phase=name,path=path_text,sha256=sha(path),origin="existing_exact_source"))
                np.testing.assert_array_equal(phase_data["cue"][1][0],state)
            else:
                initial = state.copy()
                for name,end,stimulus,samples in (
                    ("cue",1000.,cue,np.r_[np.arange(0,10.00001,.1),np.arange(11,101),np.arange(105,1001,5)]),
                    ("recovery",13000.,zero,np.r_[0.,np.arange(50,13001,50)])):
                    times, values, stats = integrate(model,grid,initial,stimulus,end,samples,pools,began+3300)
                    path = response_dir / f"{name}.npz"
                    np.savez_compressed(path,time=times,states=values,x=grid.x,stimulus=stimulus)
                    phases.append(dict(phase=name,path=str(path),sha256=sha(path),origin="new_native_integration",**stats))
                    phase_data[name] = (times,values); initial = values[-1].copy()
            response["phases"] = phases
            tc,yc = phase_data["cue"]; tr,yr = phase_data["recovery"]
            observed = {t:yc[np.flatnonzero(np.isclose(tc,t,rtol=0,atol=1e-10))[0]] for t in (1.,10.)}
            response["gain"] = cue_gain(model,grid,state,cue,core,center,peak,observed,no_cue)
            response["stage_probes"] = []
            for stage in STAGES:
                if stage < 1000:
                    snapshot = yc[np.flatnonzero(np.isclose(tc,stage,rtol=0,atol=1e-10))[0]]; stimulus = cue
                else:
                    snapshot = yr[np.flatnonzero(np.isclose(tr,stage-1000,rtol=0,atol=1e-10))[0]]; stimulus = zero
                rows, arrays = probe(model,grid,snapshot,stimulus,core,state,
                    do_tight=(fraction == .4 and stage == 1000 and job["circumference"] in (20,50)))
                path = response_dir / f"probe_t{stage:g}.npz"
                np.savez_compressed(path,**arrays)
                response["stage_probes"].append(dict(time_from_onset=stage,
                    cue_on=bool(stage < 1000),spatial=spatial(grid,snapshot[0],core,center),
                    snapshot_absolute_net_rate=float(ring.native_rhs(model,grid,snapshot,stimulus)[0]@core["weights"]),
                    rows=rows,path=str(path),sha256=sha(path)))
            response["final_spatial"] = spatial(grid,yr[-1,0],core,center)
            response["final_old_core_amount_fraction"] = float(yr[-1,0]@core["weights"]/(state[0]@core["weights"]))
            response["status"] = "complete"
            write(out / "status.json",status)
        status.update(status="complete",producer="finished",complete=True,elapsed_seconds=time.monotonic()-began)
    except Exception as error:
        status.update(status="error",producer="finished",complete=True,error=f"{type(error).__name__}: {error}",elapsed_seconds=time.monotonic()-began)
        write(out / "status.json",status); raise
    status["generated_at"] = now(); write(out / "status.json",status)
    print(json.dumps(dict(slot=slot,status=status["status"],elapsed=status.get("elapsed_seconds"))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action",choices=["initialize","submit","worker"])
    parser.add_argument("--run",type=Path,default=DEFAULT)
    parser.add_argument("--slot",type=int)
    parser.add_argument("--slots",default="0-11")
    parser.add_argument("--label",default="matrix")
    args = parser.parse_args(); run = args.run.resolve()
    if args.action == "initialize":
        manifest = initialize(run); print(json.dumps(dict(jobs=len(manifest["jobs"]),run=str(run))))
    elif args.action == "submit":
        submit(run,args.slots,args.label)
    else:
        worker(run,args.slot)


if __name__ == "__main__":
    main()
