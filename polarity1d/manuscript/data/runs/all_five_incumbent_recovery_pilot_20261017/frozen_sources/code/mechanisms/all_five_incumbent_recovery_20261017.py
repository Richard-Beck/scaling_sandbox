"""Bounded 15-case native incumbent-deactivation pilot; never launches a full grid.

All scientific integration runs from frozen sources under SLURM. Preflight is
algebraic only. Native model files and retained checkpoints are never edited.
"""
from pathlib import Path
from datetime import datetime, timezone
from dataclasses import replace
import argparse
import copy
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
from scipy.linalg import solve as dense_solve, solve_banded as scipy_solve_banded

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code/core"))
import native as n
from polarity1d.solvers import laplacian_no_flux, _jacobian_sparsity
from polarity1d.protocols import NoStimulus

DEFAULT = ROOT / "data/runs/all_five_incumbent_recovery_pilot_20261017"
TOPIC = ROOT / "figures/length_mechanisms_20261008/dev/all_five_linear_function_20261017"
METRICS = ROOT / "data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv"
MODELS = ("wave_pinning", "goryachev", "debelly", "holmes_model3", "spring")
EPSILONS = (0., .001, .005, .01)
CANDIDATES = np.logspace(-6., 0., 25)
DB_CANDIDATES = np.array([.5, 1., 2., 5., 10.])
ROUND_FACTOR = 64.


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def csvwrite(path, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def read(path):
    return json.loads(Path(path).read_text())


def band_to_dense(ab, lower, upper):
    """Independent SciPy storage definition: ab[upper+i-j,j] = A[i,j]."""
    ab = np.array(ab, copy=True)
    size = ab.shape[1]
    result = np.zeros((size, size), dtype=ab.dtype)
    for j in range(size):
        for i in range(max(0, j-upper), min(size, j+lower+1)):
            result[i, j] = ab[upper+i-j, j]
    return result


class AuditedBanded:
    def __init__(self, dense=False):
        self.dense = dense; self.calls = 0; self.maximum_scaled_residual = 0.

    def __call__(self, bandwidth, ab, b, **kwargs):
        lower, upper = bandwidth
        # Snapshot before the native overwrite_ab=True solve can mutate storage.
        ab0 = np.array(ab, dtype=float, copy=True)
        b0 = np.array(b, dtype=float, copy=True)
        matrix = band_to_dense(ab0, lower, upper)
        if self.dense:
            result = dense_solve(matrix.copy(), b0.copy(), assume_a="gen", check_finite=False)
        else:
            result = scipy_solve_banded(bandwidth, ab, b, **kwargs)
        residual = np.max(np.abs(matrix @ result - b0))
        denominator = np.linalg.norm(matrix, ord=np.inf) * np.max(np.abs(result)) + np.max(np.abs(b0))
        scaled = float(residual / max(denominator, np.finfo(float).tiny))
        self.maximum_scaled_residual = max(self.maximum_scaled_residual, scaled)
        self.calls += 1
        if scaled > 1e-12:
            raise ArithmeticError(f"Discrete linear solve scaled residual {scaled} exceeds 1e-12")
        return result


def dense_wrapper_fixture():
    # Non-symmetric and non-constant entries expose the upper/lower indexing.
    expected = np.array([[4., 2., 0., 0.], [3., 8., 5., 0.],
                         [0., 7., 12., 11.], [0., 0., 13., 18.]])
    ab = np.array([[999., 2., 5., 11.], [4., 8., 12., 18.], [3., 7., 13., 999.]])
    np.testing.assert_array_equal(band_to_dense(ab, 1, 1), expected)
    b = expected @ np.array([1., -2., 3., -.5])
    checks = []
    for use_dense in (False, True):
        audit = AuditedBanded(use_dense)
        result = audit((1, 1), ab.copy(), b.copy(), overwrite_ab=True, check_finite=False)
        np.testing.assert_allclose(result, [1., -2., 3., -.5], rtol=1e-13, atol=1e-13)
        checks.append(dict(dense=use_dense, residual=audit.maximum_scaled_residual))
    return checks


def core_material_weights(nodes, edge):
    """Exact integral of piecewise-linear nodal basis on [0, frozen edge]."""
    q = np.linspace(0., 1., nodes)
    weights = np.zeros(nodes)
    for i in range(nodes - 1):
        width = q[i+1] - q[i]
        covered = min(max(edge - q[i], 0.), width)
        weights[i] += covered - covered**2 / (2*width)
        weights[i+1] += covered**2 / (2*width)
    np.testing.assert_allclose(weights.sum(), edge, rtol=1e-13, atol=1e-14)
    return weights


def load_state(path, model):
    with np.load(path) as z:
        fields = {key[6:]: z[key][-1].copy() for key in z.files if key.startswith("field_")}
        state = {name: z["state_"+name].copy() for name in n.DB_NAMES} if model == "debelly" else z["state"].copy()
        x = z["x"].copy()
    return state, fields, x


def active(state, model):
    return state["rac"] if model == "debelly" else state[:-2:2] if model == "spring" else state[0]


def initial_geometry(case, state, fields, x):
    model = case["model"]; size = case["size_native"]
    u = active(state, model)
    physical = float(state[-1]-state[-2]) if model == "spring" else size
    q = np.linspace(0, 1, len(u)) if model == "spring" else x/size
    b0 = float(u.min()); level = b0 + .5*np.ptp(u)
    crossing = np.flatnonzero(u < level)
    if not len(crossing) or crossing[0] == 0:
        raise ValueError("Expected a qualifying boundary-connected incumbent core")
    j = int(crossing[0])
    qedge = float(q[j-1] + (level-u[j-1])/(u[j]-u[j-1])*(q[j]-q[j-1]))
    width = qedge * physical
    np.testing.assert_allclose(width, case["core_width_native"], rtol=2e-11, atol=2e-12)
    np.testing.assert_allclose(b0, case["background_native"], rtol=2e-11, atol=2e-12)
    taper = np.where(q < qedge, .5*(1+np.cos(np.pi*np.minimum(q/qedge, 1.))), 0.)
    direction = taper * (u-b0)  # Do not clip/replace the retained raw profile.
    weights = core_material_weights(len(u), qedge) if model == "spring" else np.clip(width-np.arange(len(u))*size/len(u), 0, size/len(u))
    return dict(q=q, core_edge_material=qedge, width_native=width, background=b0,
                level=level, direction=direction, core_weights=weights,
                core_sample_centers=int(np.sum(q < qedge)), physical_length=physical)


def kick(state, model, direction, epsilon):
    result = copy.deepcopy(state)
    removed = float(epsilon)*direction
    if np.any(active(state, model) < removed):
        raise ValueError("Positive active-site removal is unavailable: no clipping permitted")
    if model == "debelly":
        result["rac"] -= removed
        for name in n.DB_NAMES:
            if name != "rac":
                np.testing.assert_array_equal(result[name], state[name])
    elif model == "spring":
        result[:-2:2] -= removed; result[1:-2:2] += removed
        np.testing.assert_array_equal(result[-2:], state[-2:])
        np.testing.assert_allclose(result[:-2:2]+result[1:-2:2], state[:-2:2]+state[1:-2:2], rtol=0, atol=2e-15)
    else:
        result[0] -= removed; result[1] += removed
        np.testing.assert_allclose(result[0]+result[1], state[0]+state[1], rtol=0, atol=3e-14)
        if result.shape[0] > 2:
            np.testing.assert_array_equal(result[2:], state[2:])
    chemical = [result[name] for name in ("rac", "rho", "mca_bound", "mca_unbound")] if model == "debelly" else [result[:-2]] if model == "spring" else [result]
    if any(np.min(value) < 0 for value in chemical):
        raise ValueError("Kick produced a negative chemical field")
    if epsilon == 0:
        if isinstance(result, dict):
            for name in result: np.testing.assert_array_equal(result[name], state[name])
        else:
            np.testing.assert_array_equal(result, state)
    return result


def select_cases():
    prov = read(METRICS.parent / "report_provenance.json")
    assert sha(METRICS) == prov["summary_csv_sha256"]
    rows = list(csv.DictReader(METRICS.open())); assert len(rows) == 97
    cases, complete = [], []
    for model in MODELS:
        group = sorted([r for r in rows if r["model"] == model], key=lambda r: float(r["domain_length"]))
        selected = [group[i] for i in (0, (len(group)-1)//2, len(group)-1)]
        complete.append(dict(model=model, eligible_sizes=[float(r["domain_length"]) for r in group],
                             selected_sizes=[float(r["domain_length"]) for r in selected]))
        for label, row in zip(("short", "middle", "long"), selected):
            path = ROOT / row["checkpoint"]
            assert sha(path) == row["checkpoint_sha256"]
            cases.append(dict(slot=len(cases), model=model, label=label,
                              size_native=float(row["domain_length"]), physical_length=float(row["physical_length"]),
                              core_width_native=float(row["old_pole_width"]), background_native=float(row["background_active"]),
                              source_checkpoint=row["checkpoint"], checkpoint_sha256=row["checkpoint_sha256"]))
    assert len(cases) == 15
    return cases, complete


def preflight(destination=None):
    cases, selection = select_cases()
    audit = []
    for case in cases:
        state, fields, x = load_state(ROOT/case["source_checkpoint"], case["model"])
        geom = initial_geometry(case, state, fields, x)
        for epsilon in EPSILONS:
            probe = kick(state, case["model"], geom["direction"], epsilon)
            if epsilon:
                assert np.any(active(state, case["model"]) != active(probe, case["model"]))
        audit.append(dict(slot=case["slot"], model=case["model"], size=case["size_native"],
                          checkpoint_hash_passed=True, all_kicks_nonnegative=True,
                          core_width_native=geom["width_native"], core_sample_centers=geom["core_sample_centers"],
                          coarse_core_flag=geom["core_sample_centers"] < 8,
                          zero_kick_exact_identity=True, other_fields_and_pool_transfer_checked=True))
    result = dict(generated_at=now(), passed=True, new_integrations=0, selected_cases=15,
                  selection=selection, dense_banded_storage_fixture=dense_wrapper_fixture(), cases=audit,
                  model_geometry="qualified established linear 1D basals only",
                  current_script_sha256=sha(Path(__file__)))
    if destination: write(destination, result)
    return result


def initialize(run):
    if run.exists(): raise ValueError("Choose a fresh pilot run directory")
    check = preflight(); cases, selection = select_cases()
    for folder in ("inputs", "frozen_sources", "workers", "logs", "results"):
        (run/folder).mkdir(parents=True)
    sources = [Path(__file__)]
    sources += sorted((ROOT/"code/core").glob("*.py"))
    sources += sorted((ROOT/"code/core").glob("*.json"))
    sources += sorted((ROOT/"code/core/models").glob("*.py"))
    sources += sorted((ROOT/"code/core/src/polarity1d").glob("*.py"))
    sources += [METRICS, METRICS.parent/"report_provenance.json", TOPIC/"DESIGN.md"]
    frozen = {}
    for source in sources:
        relative = source.relative_to(ROOT); target = run/"frozen_sources"/relative
        target.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(source, target)
        frozen[str(relative)] = sha(target)
    for case in cases:
        original = ROOT/case["source_checkpoint"]
        target = run/"inputs"/f'case_{case["slot"]:02d}.npz'
        shutil.copy2(original, target); assert sha(target) == case["checkpoint_sha256"]
        case["input"] = str(target.relative_to(run))
    manifest = dict(generated_at=now(), experiment="Bounded all-five incumbent-active-deactivation pilot",
                    geometry="linear 1D exact retained basals", cases=cases, selection=selection,
                    epsilons=list(EPSILONS), primary_epsilon=.005,
                    continuous_candidates=CANDIDATES.tolist(), debelly_candidates=DB_CANDIDATES.tolist(),
                    continuous_samples=np.unique(np.r_[0., np.logspace(-7, 0, 281), CANDIDATES]).tolist(),
                    debelly_samples=np.arange(0, 10.5, .5).tolist(),
                    continuous_tight=dict(rtol=1e-10, atol=1e-12, max_step=.01),
                    continuous_tighter=dict(rtol=1e-12, atol=1e-14, max_step=.005),
                    debelly=dict(step=.5, relaxation=.98, center_boundary_motion=False,
                                 reference="Independent dense linear solves at the exact same native discrete steps; not temporal convergence"),
                    gate=dict(error_max=1e-5, linearity_relative=.05, linearity_error_allowance=3.,
                              saturation_max=.2, roundoff_factor=ROUND_FACTOR, prefix=True,
                              minimum_effect_required=False,
                              E="max(2*abs(F_tight-F_tighter),2*abs(D_zero)/D0,roundoff_floor)",
                              wording="Empirical numerical sensitivity, not a certified error bound"),
                    stopping_rule="Choose one largest passing horizon per model across all 3 sizes and amplitudes, report and stop. No full grid launch.",
                    frozen_source_sha256=frozen, preflight=check)
    write(run/"manifest.json", manifest)
    write(run/"provenance.json", dict(manifest_sha256=sha(run/"manifest.json"), generated_at=now()))
    write(run/"results/status.json", dict(generated_at=now(), producer="initialized", complete=False, jobs=15))
    return manifest


def verify(run):
    manifest = read(run/"manifest.json")
    assert sha(run/"manifest.json") == read(run/"provenance.json")["manifest_sha256"]
    for relative, digest in manifest["frozen_source_sha256"].items():
        assert sha(run/"frozen_sources"/relative) == digest, relative
    for case in manifest["cases"]:
        assert sha(run/case["input"]) == case["checkpoint_sha256"]
    return manifest


def continuous_rhs(case, state):
    model = case["model"]; size = case["size_native"]
    if model == "spring": return n.sp.rhs(0., state, L0=size, d1=80., b=4., cue=0.)
    native = n.model_for(model); cells = state.shape[1]
    x = (np.arange(cells)+.5)*size/cells
    return (np.asarray(native.diffusivities)[:, None]*laplacian_no_flux(state, size/cells)
            + native.reaction(0., x, state, np.zeros(cells)))


def integrate(case, state, manifest, reference=False):
    model = case["model"]
    if model == "debelly":
        audit = AuditedBanded(dense=reference)
        original = n.ds.solve_banded
        parameters = replace(n.ds._parameters(), center_boundary_motion=False)
        native = n.ds.DeBelly(parameters)
        try:
            n.ds.solve_banded = audit
            trajectory = n.ds.simulate_debelly(native, protocol=NoStimulus(), length=case["size_native"],
                cells=len(state["rac"]), end_time=10., times=np.asarray(manifest["debelly_samples"]),
                step=.5, initial_checkpoint=copy.deepcopy(state))
        finally:
            n.ds.solve_banded = original
        fields = {name: trajectory.fields[name].copy() for name in n.DB_NAMES}
        np.testing.assert_array_equal(fields["rac"][0], state["rac"])
        for name in n.DB_NAMES:
            np.testing.assert_array_equal(fields[name][0], state[name])
        return dict(time=trajectory.time.copy(), fields=fields,
                    numerical=dict(kind="dense_native_discrete_reference" if reference else "banded_native_discrete",
                                   calls=audit.calls, max_scaled_linear_residual=audit.maximum_scaled_residual))
    settings = manifest["continuous_tighter" if reference else "continuous_tight"]
    shape = state.shape
    fun = lambda t, y: continuous_rhs(case, y.reshape(shape)).ravel()
    options = dict(settings)
    if model != "spring": options["jac_sparsity"] = _jacobian_sparsity(shape[0], shape[1])
    solution = solve_ivp(fun, (0., 1.), state.ravel(), method="LSODA" if model == "spring" else "BDF",
                         t_eval=np.asarray(manifest["continuous_samples"]), **options)
    if not solution.success: raise RuntimeError(solution.message)
    states = solution.y.T.reshape((len(solution.t),)+shape)
    if model == "spring":
        fields = dict(active=states[:, :-2:2], inactive=states[:, 1:-2:2],
                      length=states[:, -1]-states[:, -2],
                      left_endpoint=states[:, -2], right_endpoint=states[:, -1])
    else:
        fields = {name: states[:, i] for i, name in enumerate(n.model_for(model).state_names)}
    if not np.isfinite(states).all():
        raise FloatingPointError("Continuous native trajectory has nonfinite state")
    chemical = states[:, :-2] if model == "spring" else states
    if float(np.min(chemical)) < -max(50*settings["atol"], 1e-7):
        raise FloatingPointError("Continuous native trajectory has materially negative concentration")
    return dict(time=solution.t, fields=fields, states=states,
                numerical=dict(kind="tighter_continuous" if reference else "tight_continuous",
                               nfev=solution.nfev, njev=solution.njev, nlu=solution.nlu,
                               minimum_raw_chemical_value=float(np.min(chemical)),
                               clipping_applied=False))


def amounts(trial, case, geometry):
    fields = trial["fields"]; model = case["model"]
    name = "rac" if model == "debelly" else "rac_active" if model == "holmes_model3" else "active"
    a = fields[name]
    if model == "spring":
        length = fields["length"]
        core = (a @ geometry["core_weights"]) * length
        whole = np.trapezoid(a, np.linspace(0, 1, a.shape[1]), axis=1) * length
    else:
        core = a @ geometry["core_weights"]
        whole = np.sum(a, axis=1)*case["size_native"]/a.shape[1]
    return core, whole


def pools(trial, case):
    f = trial["fields"]; model = case["model"]
    groups = [("rac_active", "rac_inactive"), ("rho_active", "rho_inactive")] if model == "holmes_model3" else [("mca_bound", "mca_unbound")] if model == "debelly" else [("active", "inactive")]
    values = []
    for group in groups:
        density = sum(f[name] for name in group)
        if model == "spring":
            amount = np.trapezoid(density, np.linspace(0, 1, density.shape[1]), axis=1)*f["length"]
        else: amount = density.sum(axis=1)*case["size_native"]/density.shape[1]
        error = float(np.max(np.abs(amount/amount[0]-1)))
        if error > 1e-8: raise ArithmeticError(f"Tracked pool drift {group}: {error}")
        values.append(dict(fields=list(group), initial=float(amount[0]), max_relative_error=error))
    return values


def amount_rates(trial, case, geometry):
    core, whole = amounts(trial, case, geometry)
    if case["model"] == "debelly":
        core_rate = np.r_[np.nan, np.diff(core)/np.diff(trial["time"])]
        whole_rate = np.r_[np.nan, np.diff(whole)/np.diff(trial["time"])]
        return core_rate, whole_rate
    rates = []
    for state in trial["states"]:
        derivative = continuous_rhs(case, state)
        if case["model"] == "spring":
            length = state[-1]-state[-2]; length_rate = derivative[-1]-derivative[-2]
            a, adot = state[:-2:2], derivative[:-2:2]
            rates.append((float(adot @ geometry["core_weights"])*length + float(a @ geometry["core_weights"])*length_rate,
                          float(np.trapezoid(adot, np.linspace(0, 1, len(a))))*length + float(np.trapezoid(a, np.linspace(0, 1, len(a))))*length_rate))
        else:
            rates.append((float(derivative[0] @ geometry["core_weights"]), float(derivative[0].sum()*case["size_native"]/state.shape[1])))
    return np.asarray(rates).T


def worker(run, slot):
    manifest = verify(run); case = manifest["cases"][slot]
    state, fields, x = load_state(run/case["input"], case["model"])
    geometry = initial_geometry(case, state, fields, x)
    folder = run/"workers"/f"{slot:02d}"; folder.mkdir(exist_ok=True)
    write(folder/"status.json", dict(generated_at=now(), producer="running", complete=False, case=case))
    start = time.monotonic(); trials = {}
    for reference in (False, True):
        prefix = "tighter" if reference else "tight"
        # The zero-removal path executes independently of the unperturbed sham.
        trials[prefix+"_sham"] = integrate(case, copy.deepcopy(state), manifest, reference)
        for epsilon in EPSILONS:
            probe = kick(state, case["model"], geometry["direction"], epsilon)
            name = prefix+f"_eps{epsilon:g}"
            trials[name] = integrate(case, probe, manifest, reference)
    arrays = dict(x=x, direction=geometry["direction"], core_weights=geometry["core_weights"],
                  initial_core_edge_material=geometry["core_edge_material"])
    pool_checks, field_differences = {}, {}
    for name, trial in trials.items():
        core, whole = amounts(trial, case, geometry)
        cr, wr = amount_rates(trial, case, geometry)
        arrays[name+"_time"] = trial["time"]
        arrays[name+"_core_amount"] = core; arrays[name+"_whole_amount"] = whole
        arrays[name+"_core_rate"] = cr; arrays[name+"_whole_rate"] = wr
        for field, value in trial["fields"].items(): arrays[name+"_field_"+field] = value
        pool_checks[name] = pools(trial, case)
    for epsilon in EPSILONS:
        primary, reference = trials[f"tight_eps{epsilon:g}"], trials[f"tighter_eps{epsilon:g}"]
        field_differences[f"eps{epsilon:g}"] = {name: float(np.max(np.abs(primary["fields"][name]-reference["fields"][name]))) for name in primary["fields"]}
    for prefix in ("tight", "tighter"):
        sham = trials[prefix+"_sham"]
        zero = trials[prefix+"_eps0"]
        for field in sham["fields"]:
            np.testing.assert_array_equal(sham["fields"][field], zero["fields"][field])
    temporary = folder/"curves.tmp.npz"; np.savez_compressed(temporary, **arrays); temporary.replace(folder/"curves.npz")
    result = dict(generated_at=now(), producer="finished", complete=True, status="complete", case=case,
                  wall_seconds=time.monotonic()-start, curves_sha256=sha(folder/"curves.npz"),
                  geometry={k: v for k, v in geometry.items() if np.isscalar(v)},
                  zero_kick_full_field_identity=True, tracked_pool_checks=pool_checks,
                  primary_reference_full_field_maxabs=field_differences,
                  numerical={name: trial["numerical"] for name, trial in trials.items()})
    write(folder/"status.json", result)


def evaluate_case(case, archive, manifest):
    model = case["model"]; times = archive["tight_sham_time"]
    candidates = manifest["debelly_candidates" if model == "debelly" else "continuous_candidates"]
    curves, rows = [], []
    stats = {}
    for epsilon in EPSILONS[1:]:
        label = f"eps{epsilon:g}"
        As = archive["tight_sham_core_amount"]; Ap = archive["tight_"+label+"_core_amount"]
        Ars = archive["tighter_sham_core_amount"]; Arp = archive["tighter_"+label+"_core_amount"]
        D = As-Ap; Dr = Ars-Arp
        Dzero = As-archive["tight_eps0_core_amount"]
        F = 1-D/D[0]; Fr = 1-Dr/Dr[0]
        roundoff = ROUND_FACTOR*np.finfo(float).eps*(np.abs(As)+np.abs(Ap)+np.abs(As[0])+np.abs(Ap[0]))/D[0]
        E = np.maximum.reduce([2*np.abs(F-Fr), 2*np.abs(Dzero)/D[0], roundoff])
        whole = archive["tight_sham_whole_amount"]-archive["tight_"+label+"_whole_amount"]
        stats[epsilon] = dict(F=F, E=E)
        for i, t in enumerate(times):
            curves.append(dict(model=model, size_native=case["size_native"], label=case["label"], epsilon=epsilon,
                time_native=float(t), initial_core_damage=float(D[0]), core_deficit=float(D[i]),
                whole_active_deficit=float(whole[i]), outside_core_deficit=float(whole[i]-D[i]),
                restored_fraction=float(F[i]), empirical_numerical_sensitivity=float(E[i]),
                conditional_K_per_native_time=float(F[i]/t) if t else None,
                sham_core_change=float(As[i]-As[0]), probe_core_change=float(Ap[i]-Ap[0]),
                sham_core_net_rate=float(archive["tight_sham_core_rate"][i]) if np.isfinite(archive["tight_sham_core_rate"][i]) else None,
                probe_core_net_rate=float(archive["tight_"+label+"_core_rate"][i]) if np.isfinite(archive["tight_"+label+"_core_rate"][i]) else None,
                change_resolved_above_3E=bool(abs(F[i]) > 3*E[i])))
    for horizon in candidates:
        at = int(np.argmin(abs(times-horizon))); assert np.isclose(times[at], horizon, rtol=0, atol=1e-14)
        prefix = (times > 0) & (times <= times[at])
        error = max(float(s["E"][prefix].max()) for s in stats.values())
        saturation = max(float((abs(s["F"])+s["E"])[prefix].max()) for s in stats.values())
        linearity = True; excess = -np.inf
        for i, eps1 in enumerate(EPSILONS[1:]):
            for eps2 in EPSILONS[i+2:]:
                a, b = stats[eps1], stats[eps2]
                tolerance = .05*np.maximum(abs(a["F"]), abs(b["F"]))+3*(a["E"]+b["E"])
                difference = abs(a["F"]-b["F"])-tolerance
                excess = max(excess, float(difference[prefix].max()))
                linearity = linearity and bool(np.all(difference[prefix] <= 0))
        rows.append(dict(model=model, size_native=case["size_native"], label=case["label"], horizon_native=float(horizon),
                         prefix_error_max=error, prefix_saturation_max=saturation,
                         prefix_linearity_max_excess=excess, error_pass=error <= 1e-5,
                         saturation_pass=saturation <= .2, linearity_pass=linearity,
                         passes=bool(error <= 1e-5 and saturation <= .2 and linearity),
                         primary_F=float(stats[.005]["F"][at]), primary_K=float(stats[.005]["F"][at]/horizon),
                         primary_E=float(stats[.005]["E"][at]),
                         primary_change_resolved_above_3E=bool(abs(stats[.005]["F"][at]) > 3*stats[.005]["E"][at])))
    return curves, rows


def collect(run):
    manifest = verify(run); curves = []; candidates = []; statuses = []
    for case in manifest["cases"]:
        folder = run/"workers"/f'{case["slot"]:02d}'
        status = read(folder/"status.json"); assert status["complete"] and status["status"] == "complete"
        assert sha(folder/"curves.npz") == status["curves_sha256"]
        with np.load(folder/"curves.npz") as z:
            a, b = evaluate_case(case, z, manifest); curves += a; candidates += b
        statuses.append(status)
    selection = []
    for model in MODELS:
        horizons = manifest["debelly_candidates" if model == "debelly" else "continuous_candidates"]
        passes = [h for h in horizons if all(r["passes"] for r in candidates if r["model"] == model and r["horizon_native"] == h)]
        chosen = max(passes) if passes else None
        selection.append(dict(model=model, selected_horizon_native=chosen,
                              status="available" if chosen is not None else "no_common_pilot_clock",
                              passing_horizons=passes,
                              primary_results=[r for r in candidates if r["model"] == model and r["horizon_native"] == chosen]))
    csvwrite(run/"results/short_curves.csv", curves); csvwrite(run/"results/candidate_gates.csv", candidates)
    write(run/"results/clock_decision.json", dict(generated_at=now(), complete=True, models=selection,
          interpretation="One finite-window clock per model; neither an instantaneous derivative nor generic stability.",
          error_interpretation="Empirical tolerance/dense-arithmetic sensitivity, not a certified or temporal-discretization error bound."))
    write(run/"results/verification.json", dict(generated_at=now(), passed=True,
        completed_cases=len(statuses), zero_kick_all_fields_identity=True,
        maximum_pool_relative_drift=max(item["max_relative_error"] for s in statuses for checks in s["tracked_pool_checks"].values() for item in checks),
        maximum_debelly_scaled_linear_residual=max((v.get("max_scaled_linear_residual", 0) for s in statuses for v in s["numerical"].values()), default=0),
        statuses=statuses))
    write(run/"results/status.json", dict(generated_at=now(), complete=True, producer="finished", jobs=15,
        stopping_rule="Pilot complete; stop before full grid, await root review."))
    return selection


def launch(run):
    manifest = verify(run)
    if (run/"launch.json").exists(): raise ValueError("Pilot already launched")
    frozen = run/"frozen_sources/code/mechanisms"/Path(__file__).name
    command = shlex.join([sys.executable, str(frozen), "worker", "--run", str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    options = ["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=3G", "--time=00:30:00",
               "--job-name=incumbent-recovery-pilot", "--array=0-14", "--chdir="+str(run),
               "--output="+str(run/"logs/pilot_%A_%a.out"), "--error="+str(run/"logs/pilot_%A_%a.err")]
    job = subprocess.check_output(options+["--wrap", command], text=True).strip().split(";")[0]
    write(run/"launch.json", dict(generated_at=now(), job_id=job, tasks=15, frozen_script_sha256=sha(frozen),
          command=command, stop_after="Pilot report and fresh critique; no full grid."))
    write(run/"results/status.json", dict(generated_at=now(), complete=False, producer="submitted", job_id=job, jobs=15))
    return job


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("action", choices=("preflight", "initialize", "launch", "worker", "collect"))
    parser.add_argument("--run", type=Path, default=DEFAULT); parser.add_argument("--slot", type=int)
    args = parser.parse_args(); run = args.run.resolve()
    if args.action == "preflight": print(json.dumps(preflight(TOPIC/"supporting/metadata/ACTIVE_PILOT_PREFLIGHT.json"), indent=2))
    elif args.action == "initialize": print(json.dumps({"cases": len(initialize(run)["cases"]), "run": str(run)}, indent=2))
    elif args.action == "launch": print(launch(run))
    elif args.action == "worker": worker(run, args.slot)
    elif args.action == "collect": print(json.dumps(collect(run), indent=2))


if __name__ == "__main__":
    main()
