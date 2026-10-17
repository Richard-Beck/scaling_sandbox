"""Common active-deactivation recovery on exact retained WP/G/H ring states.

Sources and inputs are frozen before SLURM integration.  Pilot and full-grid
launches are distinct actions; a passing pilot never launches the grid itself.
No retained model, cue history, checkpoint, or legacy assay is changed.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import time

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np
from scipy.integrate import solve_ivp
import ring_distance_prototype_20261017 as ring
import chemical_ring_distance_pilot_20261017 as old

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "data/runs/vision_ring_recovery_20261017"
SOURCE = ROOT / "data/runs/chemical_ring_functional_v2_20261017"
TOPIC = ROOT / "figures/length_mechanisms_20261008/dev/vision_ring_20261017"
MODELS = ring.KEYS
TAUS = dict(wave_pinning=.1778279410038923, goryachev=.005623413251903491,
            holmes_model3=.1778279410038923)
EPSILONS = (0., .001, .005, .01)
PYTHON = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, obj):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False,
        default=lambda x: x.item() if isinstance(x, np.generic) else str(x)) + "\n")
    temporary.replace(path)


def csvwrite(path, rows):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns); w.writeheader(); w.writerows(rows)


def geometry(grid, basal, center=None):
    """Frozen basal half-contrast interval and independent left/right taper."""
    core = old.original_core(grid, basal[0])
    if center is not None:
        # An explicit center is used only by the exact reflection fixture.
        arc = grid.arc(grid.x, center); order = np.argsort(arc)
        x, u = arc[order], basal[0, order]
        level = float(u.min()+.5*np.ptp(u)); inside = np.flatnonzero(u >= level)
        a, b = int(inside[0]), int(inside[-1])
        left = x[a-1]+(level-u[a-1])/(u[a]-u[a-1])*(x[a]-x[a-1])
        right = x[b]+(level-u[b])/(u[b+1]-u[b])*(x[b+1]-x[b])
        core = dict(origin=float(center), left=float(left), right=float(right),
                    width=float(right-left), level=level,
                    peak_cell=int(np.argmax(basal[0])),
                    weights=old.weights_about(grid, center, left, right))
    arc = grid.arc(grid.x, core["origin"])
    radius = np.where(arc < 0, -core["left"], core["right"])
    ratio = np.abs(arc)/radius
    taper = np.where(ratio < 1, .5*(1+np.cos(np.pi*np.minimum(ratio, 1))), 0.)
    background = float(basal[0].min())
    direction = taper*(basal[0]-background)
    assert np.isfinite(direction).all() and np.min(direction) >= 0
    np.testing.assert_allclose(core["weights"].sum(), core["width"], atol=2e-13)
    return dict(core, background=background, direction=direction, taper=taper)


def kick(state, direction, epsilon):
    removed = epsilon*direction
    if np.any(state[0] < removed):
        raise ValueError("Unavailable frozen basal direction: negative active field; no clipping")
    result = state.copy(); result[0] -= removed; result[1] += removed
    np.testing.assert_allclose(result[0]+result[1], state[0]+state[1], rtol=0, atol=3e-14)
    np.testing.assert_array_equal(result[2:], state[2:])
    assert np.isfinite(result).all() and result.min() >= 0
    if epsilon == 0: np.testing.assert_array_equal(result, state)
    return result


def algebraic_preflight():
    """Exact even-reflection, FV, local-kinetics and periodic-mode checks."""
    import all_five_incumbent_recovery_20261017 as linear
    rows = list(csv.DictReader((ROOT / "data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv").open()))
    results = []
    for model in MODELS:
        row = next(r for r in rows if r["model"] == model and float(r["domain_length"]) == 15.)
        path = ROOT/row["checkpoint"]; assert sha(path) == row["checkpoint_sha256"]
        state, fields, x = linear.load_state(path, model)
        case = dict(model=model, size_native=15., core_width_native=float(row["old_pole_width"]), background_native=float(row["background_active"]))
        lg = linear.initial_geometry(case, state, fields, x)
        mirrored = np.concatenate([state, state[:, ::-1]], axis=1)
        grid = ring.RingGrid(30., 2*state.shape[1]); native = ring.native_model(model)
        geom = geometry(grid, mirrored, center=0.)
        mirror = lambda a: np.concatenate([a, a[..., ::-1]], axis=-1)
        np.testing.assert_allclose(geom["direction"], mirror(lg["direction"]), atol=2e-13, rtol=2e-13)
        np.testing.assert_allclose(geom["weights"], mirror(lg["core_weights"]), atol=2e-13, rtol=2e-13)
        linear_rhs = linear.continuous_rhs(case, state)
        ring_rhs = ring.native_rhs(native, grid, mirrored, np.zeros(grid.cells))
        operator_error = float(np.max(np.abs(ring_rhs-mirror(linear_rhs))))
        np.testing.assert_allclose(ring_rhs, mirror(linear_rhs), atol=2e-12, rtol=2e-12)
        for epsilon in EPSILONS:
            np.testing.assert_allclose(kick(mirrored, geom["direction"], epsilon),
                mirror(linear.kick(state, model, lg["direction"], epsilon)), atol=2e-13, rtol=2e-13)
        # Ring total support6 is the reflection of linear end support3.
        lin_cue = .5*(1+np.cos(np.pi*np.minimum(x/3., 1))) * (x < 3)
        seam_distance = grid.distance(grid.x, 0.)
        mirrored_cue = np.where(seam_distance < 3., .5*(1+np.cos(np.pi*seam_distance/3.)), 0.)
        np.testing.assert_allclose(mirrored_cue, mirror(lin_cue), atol=2e-14)
        # Constant mode has zero diffusion; a nonconstant Fourier mode has
        # the exact discrete eigenvalue, including seam-face orientation.
        np.testing.assert_array_equal(ring.periodic_laplacian(np.ones(grid.cells), grid.dx), np.zeros(grid.cells))
        mode = np.cos(2*np.pi*np.arange(grid.cells)/grid.cells)
        eigenvalue = -4*np.sin(np.pi/grid.cells)**2/grid.dx**2
        np.testing.assert_allclose(ring.periodic_laplacian(mode, grid.dx), eigenvalue*mode, atol=2e-12, rtol=2e-10)
        for cue in (np.zeros(grid.cells), mirrored_cue):
            diff = np.asarray(native.diffusivities)[:, None]*ring.periodic_laplacian(mirrored, grid.dx)
            np.testing.assert_allclose(ring.native_rhs(native, grid, mirrored, cue)-diff,
                                       native.reaction(0., grid.x, mirrored, cue), atol=2e-13)
        results.append(dict(model=model, linear_cells=state.shape[1], ring_cells=grid.cells,
            linear_checkpoint=str(path), sha256=sha(path), mirrored_rhs_maxabs=operator_error,
            core_weights_and_direction_exact_reflection=True, cue_half_width3_ring_total_width6=True,
            conservative_kick_exact_reflection=True, constant_and_fourier_modes_passed=True))
    result = dict(generated_at=now(), passed=True, integrations=0, cases=results,
                  interpretation="Exact matched-mesh operator reflection; separately prepared ring basals are not remapped.")
    write(TOPIC/"supporting/metadata/RECOVERY_ALGEBRAIC_PREFLIGHT.json", result)
    return result


def initialize(run):
    if run.exists(): raise ValueError("Preserve existing run; choose a fresh directory")
    algebraic = algebraic_preflight()
    matrix = list(csv.DictReader((SOURCE/"results/matrix.csv").open()))
    counts = {s: sum(r["status"] == s for r in matrix) for s in sorted({r["status"] for r in matrix})}
    assert counts == dict(complete=652, geometric_overlap_exclusion=198, unqualified_preparation=17)
    for name in ("inputs", "frozen_sources", "workers", "logs", "results"):
        (run/name).mkdir(parents=True)
    sources = [Path(__file__), Path(__file__).with_name("ring_distance_prototype_20261017.py"),
               Path(__file__).with_name("chemical_ring_distance_pilot_20261017.py"),
               Path(__file__).with_name("length_dev_output_layout.py"),
               Path(__file__).with_name("all_five_incumbent_recovery_20261017.py")]
    sources += sorted((ROOT/"code/core").glob("*.py"))+sorted((ROOT/"code/core").glob("*.json"))
    sources += sorted((ROOT/"code/core/models").glob("*.py"))+sorted((ROOT/"code/core/src/polarity1d").glob("*.py"))
    frozen = {}
    for path in sources:
        rel = path.relative_to(ROOT); target = run/"frozen_sources"/rel
        target.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(path, target); frozen[str(rel)] = sha(target)
    jobs, cases, external = [], [], {}
    source_manifest = read(SOURCE/"manifest.json")
    for base in source_manifest["jobs"]:
        status_path = SOURCE/"workers"/f'{base["slot"]:02d}'/"status.json"
        status = read(status_path); external[str(status_path)] = sha(status_path)
        if status["status"] == "unqualified_preparation":
            assert base["model"] == "holmes_model3" and base["circumference"] == 80
            continue
        assert status["complete"] and status["status"] == "complete"
        checkpoint = Path(status["checkpoint"]); assert sha(checkpoint) == status["checkpoint_sha256"]
        job = dict(slot=len(jobs), model=base["model"], circumference=base["circumference"], cells=base["cells"], cases=[])
        basal = run/"inputs"/f'basal_{job["slot"]:02d}.npz'; shutil.copy2(checkpoint, basal)
        external[str(checkpoint)] = sha(checkpoint)
        common = dict(model=job["model"], circumference=job["circumference"], cells=job["cells"],
            size_native=job["circumference"], basal_input=str(basal.relative_to(run)), basal_sha256=sha(basal),
            source_status=str(status_path), source_status_sha256=sha(status_path), tau=TAUS[job["model"]],
            worker_slot=job["slot"])
        basal_case = dict(common, id=len(cases), stage="basal", fraction=None, short_route=None,
            time_from_onset=0., cue_on=False, input=str(basal.relative_to(run)), input_key="state",
            input_sha256=sha(basal), source_input=str(checkpoint), source_input_sha256=sha(checkpoint))
        job["cases"].append(basal_case["id"]); cases.append(basal_case)
        for response in status["responses"]:
            if response["status"] != "complete": continue
            probe = next(p for p in response["stage_probes"] if p["time_from_onset"] == 10.)
            source = Path(probe["path"]); assert sha(source) == probe["sha256"]
            target = run/"inputs"/f'stage_{len(cases):04d}.npz'; shutil.copy2(source, target)
            external[str(source)] = sha(source)
            case = dict(common, id=len(cases), stage="challenged_t10", fraction=response["fraction"],
                short_route=response["short_route"], long_route=response["long_route"], time_from_onset=10.,
                cue_on=True, cue_center_cell=response["cue_center_cell"],
                stage_identity=probe["spatial"]["location_association"],
                distinct_old_peak_available=probe["spatial"]["distinct_old_peak_available"],
                input=str(target.relative_to(run)), input_key="snapshot", input_sha256=sha(target),
                source_input=str(source), source_input_sha256=sha(source))
            job["cases"].append(case["id"]); cases.append(case)
        jobs.append(job)
    assert len(jobs) == 50 and len(cases) == 702
    pilot = []
    for model in MODELS:
        group = [c for c in cases if c["model"] == model]
        largest = max(c["circumference"] for c in group)
        pilot += [next(c["id"] for c in group if c["stage"] == "basal" and c["circumference"] == size) for size in (10, 30, largest)]
        challenged = [c for c in group if c["stage"] == "challenged_t10"]
        pilot.append(min((c for c in challenged if c["circumference"] == 10), key=lambda c:c["fraction"])["id"])
        pilot.append(next(c["id"] for c in challenged if c["circumference"] == largest and c["fraction"] == .5))
    manifest = dict(generated_at=now(), experiment="Common incumbent active-deactivation recovery: ring basal and exact challenged t10",
        geometry="Native periodic WP/G/H, dx0.05; frozen original basal core and direction",
        jobs=jobs, cases=cases, planned_cells=matrix, planned_counts=counts, pilot_case_ids=pilot,
        taus=TAUS, epsilons=list(EPSILONS), primary_epsilon=.005,
        continuous_tight=dict(rtol=1e-10, atol=1e-12, max_step=.01),
        continuous_tighter=dict(rtol=1e-12, atol=1e-14, max_step=.005),
        gate=dict(error_max=1e-5, linearity_relative=.05, linearity_error_allowance=3., saturation_max=.2,
            roundoff_factor=64., prefix=True, minimum_effect_required=False,
            E="max(2|F_tight-F_tighter|,2|D_zero|/D0,64epsmach(|As|+|Ap|+|As0|+|Ap0|)/D0)"),
        outputs="Dense scalar core/whole/rate traces for every trial; complete initial/final all-field states; five primary/sham all-field snapshots.",
        source_sha256=frozen, external_sha256=external, algebraic_preflight=algebraic,
        source_matrix_sha256=sha(SOURCE/"results/matrix.csv"),
        interpretation="Perturbation-specific finite-window recovery, not generic stability; unchanged cue continues during challenged assay.")
    write(run/"manifest.json", manifest); write(run/"provenance.json", dict(manifest_sha256=sha(run/"manifest.json")))
    write(run/"results/status.json", dict(generated_at=now(), complete=False, producer="initialized", next="Await reviewed design before pilot launch"))
    return dict(cases=len(cases), workers=len(jobs), pilot=[cases[i] for i in pilot])


def verify(run):
    manifest = read(run/"manifest.json")
    assert sha(run/"manifest.json") == read(run/"provenance.json")["manifest_sha256"]
    for path, digest in manifest["source_sha256"].items(): assert sha(run/"frozen_sources"/path) == digest, path
    return manifest


def load_case(run, case):
    for key in ("input", "basal_input"):
        assert sha(run/case[key]) == case["input_sha256" if key == "input" else "basal_sha256"]
    with np.load(run/case["input"]) as z: state = z[case["input_key"]].copy()
    with np.load(run/case["basal_input"]) as z: basal = z["state"].copy()
    grid = ring.RingGrid(case["circumference"], case["cells"]); native = ring.native_model(case["model"])
    geom = geometry(grid, basal)
    cue = ring.ring_cue(grid, case["cue_center_cell"], 3., old.PEAKS[case["model"]]) if case["cue_on"] else np.zeros(grid.cells)
    for epsilon in EPSILONS: kick(state, geom["direction"], epsilon)
    return state, basal, grid, native, geom, cue


def samples(tau):
    return np.unique(np.r_[0., np.logspace(-7, np.log10(tau), 121),
                           np.logspace(-6, 0, 25)[np.logspace(-6, 0, 25) < tau], tau])


def integrate(initial, grid, native, geom, cue, times, settings):
    shape = initial.shape
    rhs = lambda t,y:ring.native_rhs(native, grid, y.reshape(shape), cue).ravel()
    solution = solve_ivp(rhs, (0., float(times[-1])), initial.ravel(), method="BDF", t_eval=times,
                         jac_sparsity=ring.periodic_sparsity(shape[0], grid.cells), **settings)
    if not solution.success: raise RuntimeError(solution.message)
    states = solution.y.T.reshape((len(times),)+shape)
    if not np.isfinite(states).all() or states.min() < -max(50*settings["atol"], 1e-7):
        raise FloatingPointError("Nonfinite/materially negative state; no clipping")
    core = states[:, 0]@geom["weights"]; whole = states[:, 0].sum(axis=1)*grid.dx
    core_rate, whole_rate, reaction_rate, diffusion_rate = [], [], [], []
    for state in states:
        reaction = native.reaction(0., grid.x, state, cue)[0]
        diffusion = native.diffusivities[0]*ring.periodic_laplacian(state[0], grid.dx)
        reaction_rate.append(float(reaction@geom["weights"])); diffusion_rate.append(float(diffusion@geom["weights"]))
        core_rate.append(reaction_rate[-1]+diffusion_rate[-1]); whole_rate.append(float((reaction+diffusion).sum()*grid.dx))
    pools = np.array([ring.pool_amounts(native, grid, state) for state in states])
    drift = float(np.max(np.abs(pools/pools[0]-1)))
    if drift > 1e-8: raise ArithmeticError(f"Conserved pool drift {drift}")
    return dict(states=states, core_amount=core, whole_amount=whole,
                core_rate=np.array(core_rate), whole_rate=np.array(whole_rate),
                core_reaction_rate=np.array(reaction_rate), core_diffusion_rate=np.array(diffusion_rate),
                pools=pools, numerical=dict(nfev=solution.nfev, njev=solution.njev, nlu=solution.nlu,
                    minimum=float(states.min()), pool_max_relative_drift=drift))


def run_case(run, case, manifest):
    folder = run/"workers"/f'{case["worker_slot"]:02d}'/f'case_{case["id"]:04d}'
    if (folder/"status.json").exists():
        status = read(folder/"status.json")
        if status.get("complete") and status.get("status") == "complete":
            assert sha(folder/"curves.npz") == status["curves_sha256"]; return status
    folder.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    write(folder/"status.json", dict(generated_at=now(), complete=False, producer="running", case=case))
    try:
        state, basal, grid, native, geom, cue = load_case(run, case)
    except ValueError as e:
        status = dict(generated_at=now(), complete=True, producer="finished", status="unavailable_kick", case=case, reason=str(e))
        write(folder/"status.json", status); return status
    times = samples(case["tau"])
    arrays = dict(time=times, x=grid.x, direction=geom["direction"], taper=geom["taper"],
                  core_weights=geom["weights"], basal_state=basal, source_state=state, stimulus=cue)
    numerical = {}; zero_identity = {}; maximum_field_reference_difference = {}
    for precision in ("tight", "tighter"):
        settings = manifest["continuous_"+precision]
        sham_states = None
        for label, epsilon in [("sham", None)]+[(f"eps{e:g}", e) for e in EPSILONS]:
            initial = state.copy() if epsilon is None else kick(state, geom["direction"], epsilon)
            trial = integrate(initial, grid, native, geom, cue, times, settings)
            name = precision+"_"+label
            for key in ("core_amount", "whole_amount", "core_rate", "whole_rate", "core_reaction_rate", "core_diffusion_rate", "pools"):
                arrays[name+"_"+key] = trial[key]
            arrays[name+"_initial_state"] = trial["states"][0]
            arrays[name+"_final_state"] = trial["states"][-1]
            if label in ("sham", "eps0.005"):
                indices = np.unique([0, len(times)//4, len(times)//2, 3*len(times)//4, len(times)-1])
                arrays[name+"_snapshot_times"] = times[indices]
                arrays[name+"_snapshot_states"] = trial["states"][indices]
            if epsilon is None: sham_states = trial["states"].copy()
            if epsilon == 0:
                np.testing.assert_array_equal(sham_states, trial["states"]); zero_identity[precision] = True
            numerical[name] = trial["numerical"]
            del trial
    for label in ["sham"]+[f"eps{e:g}" for e in EPSILONS]:
        maximum_field_reference_difference[label] = float(np.max(np.abs(arrays["tight_"+label+"_final_state"]-arrays["tighter_"+label+"_final_state"])))
    temporary = folder/"curves.tmp.npz"; np.savez_compressed(temporary, **arrays); temporary.replace(folder/"curves.npz")
    result = dict(generated_at=now(), complete=True, producer="finished", status="complete", case=case,
        curves_sha256=sha(folder/"curves.npz"), wall_seconds=time.monotonic()-start,
        geometry={k:v for k,v in geom.items() if np.isscalar(v)}, numerical=numerical,
        zero_kick_full_field_identity=zero_identity, final_field_tight_tighter_maxabs=maximum_field_reference_difference)
    write(folder/"status.json", result); return result


def evaluate(case, archive, manifest):
    times = archive["time"]; stats = {}; rows = []
    for epsilon in EPSILONS[1:]:
        label = f"eps{epsilon:g}"
        As, Ap = archive["tight_sham_core_amount"], archive["tight_"+label+"_core_amount"]
        Ars, Arp = archive["tighter_sham_core_amount"], archive["tighter_"+label+"_core_amount"]
        D, Dr = As-Ap, Ars-Arp
        assert D[0] > 0 and Dr[0] > 0
        F, Fr = 1-D/D[0], 1-Dr/Dr[0]
        Dzero = As-archive["tight_eps0_core_amount"]
        roundoff = 64*np.finfo(float).eps*(abs(As)+abs(Ap)+abs(As[0])+abs(Ap[0]))/D[0]
        E = np.maximum.reduce([2*abs(F-Fr), 2*abs(Dzero)/D[0], roundoff])
        stats[epsilon] = dict(F=F, E=E)
        whole = archive["tight_sham_whole_amount"]-archive["tight_"+label+"_whole_amount"]
        row = {k:case.get(k) for k in ("id", "model", "circumference", "cells", "stage", "fraction", "short_route", "long_route", "time_from_onset", "cue_on", "stage_identity", "distinct_old_peak_available", "source_input", "source_input_sha256", "worker_slot")}
        row.update(epsilon=epsilon, tau=case["tau"], F=float(F[-1]), K=float(F[-1]/case["tau"]),
            D0=float(D[0]), Dtau=float(D[-1]), E=float(E[-1]),
            whole_active_deficit=float(whole[-1]), outside_core_deficit=float(whole[-1]-D[-1]),
            sham_core_change=float(As[-1]-As[0]), probe_core_change=float(Ap[-1]-Ap[0]),
            sham_core_amount_initial=float(As[0]), probe_core_amount_initial=float(Ap[0]),
            sham_core_amount_final=float(As[-1]), probe_core_amount_final=float(Ap[-1]),
            sham_net_rate_initial=float(archive["tight_sham_core_rate"][0]),
            sham_net_rate_final=float(archive["tight_sham_core_rate"][-1]),
            probe_net_rate_initial=float(archive["tight_"+label+"_core_rate"][0]),
            probe_net_rate_final=float(archive["tight_"+label+"_core_rate"][-1]),
            conditional_net_rate_delta_initial=float(archive["tight_"+label+"_core_rate"][0]-archive["tight_sham_core_rate"][0]),
            conditional_net_rate_delta_final=float(archive["tight_"+label+"_core_rate"][-1]-archive["tight_sham_core_rate"][-1]),
            resolved_above_3E=bool(abs(F[-1]) > 3*E[-1]),
            core_origin=float(archive["x"][np.argmax(archive["basal_state"][0])]),
            core_width=float(archive["core_weights"].sum()),
            direction_whole_integral=float(archive["direction"].sum()*case["circumference"]/case["cells"]))
        rows.append(row)
    prefix = times > 0
    error = max(float(s["E"][prefix].max()) for s in stats.values())
    saturation = max(float((abs(s["F"])+s["E"])[prefix].max()) for s in stats.values())
    excess = -np.inf
    for i, ea in enumerate(EPSILONS[1:]):
        for eb in EPSILONS[i+2:]:
            a,b = stats[ea],stats[eb]
            excess = max(excess, float((abs(a["F"]-b["F"])-.05*np.maximum(abs(a["F"]),abs(b["F"]))-3*(a["E"]+b["E"]))[prefix].max()))
    gate = dict(prefix_error_max=error, prefix_saturation_max=saturation, prefix_linearity_max_excess=excess,
        error_pass=bool(error <= 1e-5), saturation_pass=bool(saturation <= .2), linearity_pass=bool(excess <= 0),
        passes=bool(error <= 1e-5 and saturation <= .2 and excess <= 0))
    for row in rows: row.update(gate)
    return rows, gate


def worker(run, slot, pilot=False):
    assert os.environ.get("SLURM_JOB_ID"), "Integration uses SLURM"
    manifest = verify(run)
    if pilot:
        case = manifest["cases"][manifest["pilot_case_ids"][slot]]
        run_case(run, case, manifest)
    else:
        job = manifest["jobs"][slot]
        for cid in job["cases"]: run_case(run, manifest["cases"][cid], manifest)


def collect(run, pilot=False):
    manifest = verify(run); ids = manifest["pilot_case_ids"] if pilot else [c["id"] for c in manifest["cases"]]
    rows, statuses, unavailable = [], [], []
    for cid in ids:
        case = manifest["cases"][cid]; folder = run/"workers"/f'{case["worker_slot"]:02d}'/f'case_{cid:04d}'
        status = read(folder/"status.json"); assert status["complete"]
        if status["status"] != "complete": unavailable.append(status); continue
        assert sha(folder/"curves.npz") == status["curves_sha256"]
        with np.load(folder/"curves.npz") as z: result, gate = evaluate(case, z, manifest)
        rows += result; statuses.append(dict(id=cid, model=case["model"], stage=case["stage"], circumference=case["circumference"], **gate,
            curves_path=str(folder/"curves.npz"), curves_sha256=status["curves_sha256"]))
    prefix = "pilot" if pilot else "recovery"
    csvwrite(run/"results"/(prefix+".csv"), rows)
    if not pilot:
        csvwrite(run/"results/basal.csv", [r for r in rows if r["stage"] == "basal"])
        csvwrite(run/"results/challenged.csv", [r for r in rows if r["stage"] == "challenged_t10"])
    result = dict(generated_at=now(), complete=True, producer="finished", scope="15-state bounded pilot" if pilot else "50 basal and652 challenged exact states",
        passed=bool(not unavailable and all(s["passes"] for s in statuses)), completed=len(statuses), unavailable=unavailable,
        cases=statuses, csv_sha256=sha(run/"results"/(prefix+".csv")),
        instruction="Review pilot clocks before full launch" if pilot else "Retain failures; do not use per-cell clocks")
    write(run/"results"/(prefix+"_audit.json"), result)
    return result


def launch(run, pilot=False):
    manifest = verify(run)
    if not pilot:
        assert read(run/"results/pilot_audit.json")["passed"], "Joint clock review required after any pilot failure"
    label = "pilot" if pilot else "full"
    if (run/("launch_"+label+".json")).exists(): raise ValueError("Already launched; inspect existing job handles")
    script = run/"frozen_sources/code/mechanisms"/Path(__file__).name
    command = shlex.join([PYTHON, str(script), "pilot-worker" if pilot else "worker", "--run", str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    count = 15 if pilot else len(manifest["jobs"])
    args = ["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=3G", "--time=01:00:00",
        "--array=0-"+str(count-1), "--job-name=ring-active-recovery-"+label,
        "--output="+str(run/"logs"/(label+"_%A_%a.out")), "--error="+str(run/"logs"/(label+"_%A_%a.err")), "--wrap", command]
    job = subprocess.check_output(args, text=True).strip().split(";")[0]
    write(run/("launch_"+label+".json"), dict(generated_at=now(), job_id=job, workers=count, command=args))
    return job


def main():
    p = argparse.ArgumentParser(); p.add_argument("action", choices=("algebraic-preflight", "initialize", "launch-pilot", "launch-full", "pilot-worker", "worker", "collect-pilot", "collect"))
    p.add_argument("--run", type=Path, default=DEFAULT); p.add_argument("--slot", type=int)
    a = p.parse_args(); run = a.run.resolve()
    if a.action == "algebraic-preflight": result = algebraic_preflight()
    elif a.action == "initialize": result = initialize(run)
    elif a.action.startswith("launch"): result = launch(run, a.action == "launch-pilot")
    elif a.action.endswith("worker"): result = worker(run, a.slot, a.action == "pilot-worker")
    else: result = collect(run, a.action == "collect-pilot")
    if result is not None: print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
