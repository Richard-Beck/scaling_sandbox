"""Frozen basal organization/gain/common active-loss recovery on Spring rings.

This bounded addition uses the exact shared qualifying Figure 3/4 basals and
the explicitly accepted mean-force stretch/dilution/diffusion closure. All
integration uses SLURM; collection independently recomputes saved observables.
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

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np
from scipy.integrate import solve_ivp
from spring_ring_model_20261019 import Closure, physical_raised_cosine

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "data/runs/revision_spring_ring_function_20261019"
BASALS = ROOT / "data/runs/revision_ring_spring_20261019"
TOPIC = ROOT / "figures/length_mechanisms_20261008/dev/revision_program_20261019/spring_function"
PYTHON = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"
TAU = .01778279410038923
EPS = (0., .001, .005, .01)


def now(): return datetime.now(timezone.utc).isoformat()
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path): return json.loads(Path(path).read_text())


def write(path, obj):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False, default=lambda x:x.item() if isinstance(x, np.generic) else str(x))+"\n")
    temporary.replace(path)


def csvwrite(path, records):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys(k for r in records for k in r))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys); writer.writeheader(); writer.writerows(records)


def geometry(model, state):
    """Frozen basal half-contrast core, FV edges and material taper.

    Core boundaries are interpolated in periodic material position. Fractional
    FV edge weights are frozen in material coordinates; amounts include each
    sham/probe trajectory's own evolving circumference.
    """
    u, _, C = model.fields(state)
    origin = float(model.material_positions[np.argmax(u)])
    arc = (model.material_positions-origin+.5) % 1.-.5
    order = np.argsort(arc); x, active = arc[order], u[order]
    background = float(u.min()); level = background+.5*float(np.ptp(u))
    inside = np.flatnonzero(active >= level)
    assert len(inside) > 1 and np.all(np.diff(inside) == 1)
    a, b = int(inside[0]), int(inside[-1]); assert a > 0 and b < len(x)-1
    left = float(x[a-1]+(level-active[a-1])/(active[a]-active[a-1])*(x[a]-x[a-1]))
    right = float(x[b]+(level-active[b])/(active[b+1]-active[b])*(x[b+1]-x[b]))
    weights = np.maximum(0., np.minimum(arc+.5/model.cells, right)-np.maximum(arc-.5/model.cells, left))
    np.testing.assert_allclose(weights.sum(), right-left, rtol=1e-13, atol=1e-15)
    radius = np.where(arc < 0., -left, right)
    ratio = abs(arc)/radius
    taper = np.where(ratio < 1., .5*(1.+np.cos(np.pi*np.minimum(ratio, 1.))), 0.)
    direction = taper*(u-background)
    return dict(origin=origin, arc=arc, background=background, threshold=level, left_material=left,
                right_material=right, weights_material=weights, direction=direction,
                width_native=C*(right-left), basal_physical_circumference=C)


def kick(model, state, direction, epsilon):
    n = model.cells
    result = state.copy(); result[:n] -= epsilon*direction; result[n:-1] += epsilon*direction
    assert result[:-1].min() >= 0
    np.testing.assert_allclose(result[:n]+result[n:-1], state[:n]+state[n:-1], rtol=0, atol=3e-14)
    assert result[-1] == state[-1]
    np.testing.assert_allclose(model.amount(result), model.amount(state), rtol=2e-15)
    if epsilon == 0.: np.testing.assert_array_equal(result, state)
    return result


def launch(run, basals, topic):
    if run.exists(): raise FileExistsError("Preserve retained outputs; use a fresh run directory")
    for name in ("inputs", "frozen_sources", "workers", "logs", "results"):
        (run/name).mkdir(parents=True)
    for name in ("vectors", "data", "metadata", "viewers"):
        (topic/"supporting"/name).mkdir(parents=True, exist_ok=True)
    sources = [Path(__file__), Path(__file__).with_name("spring_ring_model_20261019.py"),
        Path(__file__).with_name("spring_ring_closure_20261017.py"), ROOT/"code/core/models/spring.py", ROOT/"code/core/models/__init__.py"]
    frozen = {}
    for path in sources:
        relative = path.relative_to(ROOT); target = run/"frozen_sources"/relative
        target.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(path, target); frozen[str(relative)] = sha(target)
    cases = []
    for C0 in (1., 2., 3.):
        status_path = basals/"preparations"/f"C{C0:g}"/"status.json"
        record = read(status_path)
        assert record["complete"] and record["status"] == "qualified"
        original = Path(record["checkpoint"]); assert sha(original) == record["checkpoint_sha256"]
        target = run/"inputs"/f"C{C0:g}.npz"; shutil.copy2(original, target)
        n = int(record["job"]["cells"])
        with np.load(target) as z:
            model = Closure(C0, n); geom = geometry(model, z["state"])
            for epsilon in EPS: kick(model, z["state"], geom["direction"], epsilon)
        cases.append(dict(slot=len(cases), resting_circumference=C0, cells=n, input=str(target.relative_to(run)), input_sha256=sha(target),
            original_checkpoint=str(original), source_status=str(status_path), source_status_sha256=sha(status_path),
            basal_physical_circumference=geom["basal_physical_circumference"], core_width_native=geom["width_native"]))
    manifest = dict(generated_at=now(), experiment="Spring ring basal organization, native input gain and common active-loss recovery", cases=cases,
        topic=str(topic), frozen_source_sha256=frozen, primary_tau=TAU, epsilons=EPS, primary_epsilon=.005,
        sample_times=np.unique(np.r_[0., np.logspace(-7, np.log10(TAU), 121), TAU]).tolist(),
        tight=dict(rtol=1e-10, atol=1e-12, max_step=.01), tighter=dict(rtol=1e-12, atol=1e-14, max_step=.005),
        gate=dict(error_max=1e-5, linearity_relative=.05, linearity_error_allowance=3., saturation_max=.2, prefix=True,
            E="max(2*abs(F_tight-F_tighter),2*abs(D_zero)/D0,64*machine_epsilon*amounts/D0)",
            interpretation="Empirical numerical sensitivity; not a certified error bound"),
        gain=dict(full_physical_support=.2, peak_for_identity=1., actual_competition_peak=4.,
            prescribed_fractions=[.30,.35,.40,.45,.50], dense_requested_fractions=np.linspace(0., .5, 101).tolist(),
            placement="Nearest native material cell to requested positive offset from the exact basal maximum; keep realized separation"),
        scope="Three exact shared qualifying basals only. No challenged-state recovery campaign; no changes to established interval Fig2.")
    write(run/"manifest.json", manifest)
    write(run/"results/status.json", dict(generated_at=now(), complete=False, producer="submitted", cases=3))
    runner = run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    command = shlex.join([PYTHON, str(runner), "worker", "--run", str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    args = ["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=3G", "--time=00:30:00", "--array=0-2",
        "--job-name=Spring-ring-function", "--output="+str(run/"logs/%A_%a.out"), "--error="+str(run/"logs/%A_%a.err"), "--wrap", "exec "+command]
    job = subprocess.check_output(args, text=True).strip()
    write(run/"launch.json", dict(generated_at=now(), array_job_id=job, manifest_sha256=sha(run/"manifest.json"), command=args))
    print(json.dumps(dict(run=str(run), array_job_id=job, cases=cases), indent=2))


def verify(run):
    manifest = read(run/"manifest.json")
    for path, digest in manifest["frozen_source_sha256"].items(): assert sha(run/"frozen_sources"/path) == digest
    for case in manifest["cases"]: assert sha(run/case["input"]) == case["input_sha256"]
    return manifest


def integrate(model, state, times, settings):
    solution = solve_ivp(lambda t,y:model.rhs(y,0.), (0.,float(times[-1])), state, method="BDF", t_eval=times,
        jac=lambda t,y:model.jacobian(y,0.), **settings)
    if not solution.success: raise RuntimeError(solution.message)
    states = solution.y.T
    assert np.isfinite(states).all() and states[:,:-1].min() >= -5e-10
    amount = states[:,-1]*np.mean(states[:,:model.cells]+states[:,model.cells:-1], axis=1)
    drift = float(np.max(abs(amount/amount[0]-1.)))
    assert drift <= 1e-8
    return states, dict(nfev=solution.nfev,njev=solution.njev,nlu=solution.nlu,physical_amount_relative_drift=drift)


def worker(run, slot):
    if not os.environ.get("SLURM_JOB_ID"): raise RuntimeError("Native recovery integration requires SLURM")
    manifest = verify(run); case = manifest["cases"][slot]
    folder = run/"workers"/f"{slot:02d}"; folder.mkdir(parents=True, exist_ok=True)
    write(folder/"status.json", dict(generated_at=now(), complete=False, producer="running", slurm_job_id=os.environ["SLURM_JOB_ID"]))
    with np.load(run/case["input"]) as z: state = z["state"].copy()
    model = Closure(case["resting_circumference"], case["cells"]); geom = geometry(model, state)
    times = np.asarray(manifest["sample_times"])
    arrays = dict(time=times, basal_state=state, weights_material=geom["weights_material"], direction=geom["direction"],
        material_positions=model.material_positions)
    numerics = {}
    for accuracy in ("tight", "tighter"):
        labels = [("sham", None)]+[(f"eps{epsilon:g}",epsilon) for epsilon in EPS]
        for label, epsilon in labels:
            initial = state.copy() if epsilon is None else kick(model, state, geom["direction"], epsilon)
            states, record = integrate(model, initial, times, manifest[accuracy])
            arrays[accuracy+"_"+label] = states; numerics[accuracy+"_"+label] = record
        np.testing.assert_array_equal(arrays[accuracy+"_sham"], arrays[accuracy+"_eps0"])
    archive = folder/"recovery.npz"; np.savez_compressed(archive, **arrays)
    write(folder/"status.json", dict(generated_at=now(), complete=True, producer="finished", slurm_job_id=os.environ["SLURM_JOB_ID"],
        case=case, archive=str(archive), archive_sha256=sha(archive), manifest_sha256=sha(run/"manifest.json"),
        numerical=numerics, zero_kick_all_fields_identity=True))
    print(json.dumps(dict(slot=slot, complete=True, archive=str(archive))))


def evaluate(model, arrays, tau):
    time = arrays["time"]; w = arrays["weights_material"]
    amount = lambda states:states[:,-1]*(states[:,:model.cells]@w)
    rows, stats, history = [], {}, []
    for epsilon in EPS[1:]:
        As, Ap = amount(arrays["tight_sham"]), amount(arrays[f"tight_eps{epsilon:g}"])
        Ars, Arp = amount(arrays["tighter_sham"]), amount(arrays[f"tighter_eps{epsilon:g}"])
        D, Dr = As-Ap, Ars-Arp; assert D[0] > 0 and Dr[0] > 0
        F, Fr = 1-D/D[0], 1-Dr/Dr[0]
        Dzero = As-amount(arrays["tight_eps0"])
        roundoff = 64*np.finfo(float).eps*(abs(As)+abs(Ap)+abs(As[0])+abs(Ap[0]))/D[0]
        E = np.maximum.reduce([2*abs(F-Fr),2*abs(Dzero)/D[0],roundoff])
        stats[epsilon] = dict(F=F,E=E)
        whole_s = arrays["tight_sham"][:,-1]*arrays["tight_sham"][:,:model.cells].mean(axis=1)
        whole_p = arrays[f"tight_eps{epsilon:g}"][:,-1]*arrays[f"tight_eps{epsilon:g}"][:,:model.cells].mean(axis=1)
        row = dict(epsilon=epsilon,tau=tau,F=float(F[-1]),K=float(F[-1]/tau),D0=float(D[0]),Dtau=float(D[-1]),E=float(E[-1]),
            sham_core_amount_initial=float(As[0]),probe_core_amount_initial=float(Ap[0]),
            sham_core_change=float(As[-1]-As[0]),probe_core_change=float(Ap[-1]-Ap[0]),
            whole_active_deficit=float(whole_s[-1]-whole_p[-1]),outside_core_deficit=float(whole_s[-1]-whole_p[-1]-D[-1]),
            resolved_above_3E=bool(abs(F[-1]) > 3*E[-1]))
        rows.append(row)
        for i,t in enumerate(time):
            history.append(dict(epsilon=epsilon,time_native=float(t),F=float(F[i]),E=float(E[i]),deficit=float(D[i]),
                sham_core_amount=float(As[i]),probe_core_amount=float(Ap[i])))
    prefix = time > 0
    error = max(float(s["E"][prefix].max()) for s in stats.values())
    saturation = max(float((abs(s["F"])+s["E"])[prefix].max()) for s in stats.values())
    excess = max(float((abs(stats[a]["F"]-stats[b]["F"])-.05*np.maximum(abs(stats[a]["F"]),abs(stats[b]["F"]))-3*(stats[a]["E"]+stats[b]["E"]))[prefix].max())
        for a,b in ((.001,.005),(.001,.01),(.005,.01)))
    gates = dict(prefix_error_max=error,prefix_saturation_max=saturation,prefix_linearity_max_excess=excess,
        error_pass=bool(error<=1e-5),saturation_pass=bool(saturation<=.2),linearity_pass=bool(excess<=0),
        passes=bool(error<=1e-5 and saturation<=.2 and excess<=0))
    for row in rows: row.update(gates)
    return rows, history, gates


def collect(run):
    manifest = verify(run); topic = Path(manifest["topic"])
    summaries, curves, gains, organizations, profiles, audits = [], [], [], [], [], []
    missing = []
    for case in manifest["cases"]:
        status_path = run/"workers"/f'{case["slot"]:02d}'/"status.json"
        if not status_path.exists() or not read(status_path).get("complete"):
            missing.append(case["slot"]); continue
        status = read(status_path); assert status["manifest_sha256"] == sha(run/"manifest.json")
        archive = Path(status["archive"]); assert sha(archive) == status["archive_sha256"]
        with np.load(archive) as z: arrays = {k:z[k].copy() for k in z.files}
        with np.load(run/case["input"]) as z: np.testing.assert_array_equal(arrays["basal_state"],z["state"])
        C0,n = case["resting_circumference"],case["cells"]
        model = Closure(C0,n); state=arrays["basal_state"]; geom = geometry(model,state)
        np.testing.assert_array_equal(arrays["direction"],geom["direction"])
        np.testing.assert_array_equal(arrays["weights_material"],geom["weights_material"])
        for accuracy in ("tight","tighter"):
            np.testing.assert_array_equal(arrays[accuracy+"_sham"],arrays[accuracy+"_eps0"])
            for epsilon in EPS:
                expected = kick(model,state,geom["direction"],epsilon)
                np.testing.assert_array_equal(arrays[f"{accuracy}_eps{epsilon:g}"][0],expected)
                trajectory = arrays[f"{accuracy}_eps{epsilon:g}"]
                mass = trajectory[:,-1]*np.mean(trajectory[:,:n]+trajectory[:,n:-1],axis=1)
                assert np.max(abs(mass/mass[0]-1))<=1e-8
        measured, history, gates = evaluate(model, arrays, manifest["primary_tau"])
        for row in measured: summaries.append(dict(model="spring",resting_circumference=C0,cells=n,physical_circumference=state[-1],
            archive=str(archive),archive_sha256=status["archive_sha256"],**row))
        for row in history: curves.append(dict(model="spring",resting_circumference=C0,**row))
        u,v,C=model.fields(state)
        reaction_increment = model.rhs(state,1.)[:n]-model.rhs(state,0.)[:n]
        np.testing.assert_allclose(reaction_increment,v,rtol=2e-12,atol=2e-12)
        for dense, fractions in ((False,manifest["gain"]["prescribed_fractions"]),(True,manifest["gain"]["dense_requested_fractions"])):
            for fraction in fractions:
                offset=int(round(fraction*n));realized=offset/n
                peak_cell=int(np.argmax(u));center_cell=(peak_cell+offset)%n
                signal,_=physical_raised_cosine(model,state,model.material_positions[center_cell],1.,.2)
                integral=float(signal.sum()*C/n); assert integral>0
                gain=float((signal*v).sum()/signal.sum())
                gains.append(dict(model="spring",resting_circumference=C0,cells=n,physical_circumference=C,
                    requested_fraction=fraction,realized_fraction=realized,cue_center_cell=center_cell,
                    physical_short_route=C*min(realized,1-realized),normal_angle_degrees=360*min(realized,1-realized),
                    full_physical_support=.2,cue_weight_integral=integral,cue_weighted_gain_native=gain,
                    cue_core_overlap=bool(np.any((signal>0)&(geom["weights_material"]>0))),dense_grid=dense))
        organizations.append(dict(model="spring",resting_circumference=C0,cells=n,physical_circumference=C,
            core_width_native=geom["width_native"],core_excess_active=float(C*((u-geom["background"])@geom["weights_material"])),
            core_raw_active=float(C*(u@geom["weights_material"])),peak_concentration=float(u.max()),background=float(u.min()),
            whole_active=float(C*u.mean()),whole_inactive=float(C*v.mean()),physical_amount=model.amount(state),
            checkpoint=case["original_checkpoint"],checkpoint_sha256=case["input_sha256"]))
        order=np.argsort(geom["arc"])
        for i in order:
            profiles.append(dict(model="spring",resting_circumference=C0,material_position=float(model.material_positions[i]),
                physical_arc_from_incumbent=float(C*geom["arc"][i]),active=float(u[i]),inactive=float(v[i]),
                frozen_core_weight_material=float(geom["weights_material"][i]),tapered_active_loss_direction=float(geom["direction"][i])))
        audits.append(dict(resting_circumference=C0,archive_sha256=status["archive_sha256"],
            source_checkpoint_sha256=case["input_sha256"],native_increment_identity_max_error=float(abs(reaction_increment-v).max()),
            zero_kick_all_fields_identity=True,conservative_pointwise_transfer_verified=True,exact_starting_circumference_verified=True,**gates))
    if missing:
        print(json.dumps(dict(complete=False,missing_worker_slots=missing))); return
    assert len(summaries)==9 and len(organizations)==3
    data=topic/"supporting/data";meta=topic/"supporting/metadata"
    for name,records in (("spring_basal_recovery.csv",summaries),("spring_recovery_curves.csv",curves),
            ("spring_basal_gain.csv",gains),("spring_basal_organization.csv",organizations),("spring_native_basal_profiles.csv",profiles)):
        csvwrite(data/name,records)
    verification=dict(generated_at=now(),complete=True,passed=all(r["passes"] for r in audits),fixed_clock_native=manifest["primary_tau"],
        native_recovery_cases=3,amplitudes_per_case=3,trajectories_per_case=10,cases=audits,
        exact_shared_basals=True,source_manifest_sha256=sha(run/"manifest.json"),
        scope="Basal-only full coupled Spring ring; distinct mean-force geometry analogue. Existing interval Figure2 unchanged.")
    write(run/"results/verification.json",verification);write(meta/"verification.json",verification)
    write(run/"results/status.json",dict(generated_at=now(),complete=True,producer="finished",all_fixed_clocks_qualified=verification["passed"],cases=3))
    render(topic,organizations,gains,summaries,curves)
    (topic/"README.md").write_text(dict_to_text(organizations,summaries,verification))
    artifacts=[p for p in topic.rglob("*") if p.is_file() and p.name!="manifest.json"]
    write(meta/"manifest.json",dict(generated_at=now(),run=str(run),run_manifest_sha256=sha(run/"manifest.json"),
        frozen_source_sha256=manifest["frozen_source_sha256"],artifact_sha256={str(p.relative_to(topic)):sha(p) for p in artifacts},
        collector_source_path=str(Path(__file__)),collector_source_sha256=sha(Path(__file__)),
        producer="finished",new_native_trajectories=30))
    print(json.dumps(verification,indent=2))


def render(topic,organizations,gains,recovery,curves):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(11.7,4.8))
    fig.subplots_adjust(left=.07,right=.985,top=.80,bottom=.30,wspace=.4)
    C=[r["resting_circumference"] for r in organizations]
    for key,marker,label in (("core_excess_active","o","Core excess amount"),("core_width_native","s","Physical core width"),("peak_concentration","^","Peak concentration")):
        values=np.array([r[key] for r in organizations]);axes[0].plot(C,values/values[1],marker=marker,label=label,lw=1.3,ms=4)
    axes[0].set(title="A  Basal incumbent organization",xlabel="Resting circumference",ylabel="Value / own C₀2 value");axes[0].legend(frameon=False,fontsize=7)
    selected=[r for r in gains if not r["dense_grid"]]
    for c in C:
        group=[r for r in selected if r["resting_circumference"]==c]
        axes[1].plot([r["realized_fraction"] for r in group],[r["cue_weighted_gain_native"] for r in group],marker="o",lw=1.3,ms=4,label=f"C₀{c:g}")
        for r in group:
            if r["cue_core_overlap"]:axes[1].scatter(r["realized_fraction"],r["cue_weighted_gain_native"],marker="x",s=60,color="black",zorder=5)
    axes[1].set(title="B  Native opposing-input gain",xlabel="Realized short-route fraction",ylabel="Cue-weighted inactive concentration");axes[1].legend(frameon=False,fontsize=7)
    for epsilon in EPS[1:]:
        group=[r for r in recovery if r["epsilon"]==epsilon]
        axes[2].plot([r["resting_circumference"] for r in group],[r["F"] for r in group],marker="o",lw=1.2,ms=4,label=f"ε = {100*epsilon:g}%")
    axes[2].set(title="C  Common local active-loss recovery",xlabel="Resting circumference",ylabel="Restored sham-relative deficit fraction");axes[2].legend(frameon=False,fontsize=7)
    for axis in axes:axis.grid(axis="y",color="#eeeeee",lw=.6)
    fig.suptitle("Spring ring | Basal organization, input gain and specified local restoration",fontsize=12,y=.98)
    fig.text(.53,.045,"Exact shared C₀1/2/3 basals; full evolving mean-force/stretch/dilution/diffusion dynamics. Fixed native clock τ=0.0177827941.\n"
        "Full physical cue support 0.2; × marks overlap with the original half-contrast core. Native length/time units are uncalibrated.\n"
        "Restoration includes native reaction, transport and circumference response; it is not generic stability or a physical migration law.",ha="center",fontsize=8,linespacing=1.4)
    for extension in ("png","pdf","svg"):
        path=(topic if extension=="png" else topic/"supporting/vectors")/f"spring_basal_function.{extension}"
        fig.savefig(path,dpi=200 if extension=="png" else None,bbox_inches="tight")
    plt.close(fig)


def dict_to_text(organizations,summaries,verification):
    lines=["# Spring ring basal organization and response", "", "Generated "+now()+"; all producers finished. No timed refresh.", "",
        "The exact shared C₀1/2/3 basals qualify. Native inactive concentration is the additive activation gain. A conservative local active-loss perturbation measures sham-relative deficit restoration while the full coupled chemical/circumference state evolves. This ring analogue retains the explicitly accepted mean-force/stretch/dilution/diffusion closure and does not replace the native interval model.", "",
        "![Spring basal function](spring_basal_function.png)", "", "| Resting C₀ | Actual C | Core width | Core excess | Peak | Fτ (0.5%) | Kτ | Qualified |", "|---:|---:|---:|---:|---:|---:|---:|---|"]
    for o in organizations:
        r=next(r for r in summaries if r["resting_circumference"]==o["resting_circumference"] and r["epsilon"]==.005)
        lines.append(f'| {o["resting_circumference"]:g} | {o["physical_circumference"]:.8g} | {o["core_width_native"]:.8g} | {o["core_excess_active"]:.8g} | {o["peak_concentration"]:.8g} | {r["F"]:.8g} | {r["K"]:.8g} | {r["passes"]} |')
    lines += ["", "The original periodic half-contrast core has interpolated boundaries and fractional finite-volume edge weights. A smooth cosine taper of its basal excess is reduced by 0.1%, 0.5% or 1%; the same material enters inactive at each cell and initial circumference is unchanged. Sham and probe both evolve all native fields. Core weights remain fixed in material coordinates, and each trajectory's own physical circumference enters every amount.", "",
        r"\[D(t)=A_{\mathrm{sham}}(t)-A_{\mathrm{probe}}(t),\qquad F_{\tau}=\frac{D(0)-D(\tau)}{D(0)},\qquad K_{\tau}=\frac{F_{\tau}}{\tau}.\tag{1}\]", "",
        "The fixed clock is τ=0.01778279410038923 native units. Paired BDF histories use an analytic Jacobian and two numerical tolerances. Zero-kick identity, pointwise conservative transfer, physical amount conservation, prefix empirical numerical sensitivity ≤10⁻⁵, ≤20% restoration, and 5% amplitude linearity plus numerical allowance are checked. "+("All three cases pass." if verification["passed"] else "Some cases fail; see the explicit gate rows.")+" Every complete history and gate remains in supporting/data and supporting/metadata/verification.json.", "",
        "Gain uses a raised cosine spanning 0.2 native physical units at the same native cell positions as competition, requested offsets 0.30/0.35/0.40/0.45/0.50 of circumference. An additional 101 requested offsets cover 0–0.5. The tables retain requested and realized position, physical separation, normal angle and overlap. C₀1 at fraction 0.30 overlaps the original core and is excluded from eligible competition. A larger local restoring response does not establish greater actual retention under a challenger, and a fixed physical site displacement has different angular meaning at different circumferences.", "",
        "No challenged-state recovery campaign, duration optimization, generic-stability inference or velocity law is introduced. See [recovery](supporting/data/spring_basal_recovery.csv), [gain](supporting/data/spring_basal_gain.csv), [organization](supporting/data/spring_basal_organization.csv), [native profiles](supporting/data/spring_native_basal_profiles.csv) and [numerical verification](supporting/metadata/verification.json)."]
    return "\n".join(lines)+"\n"


def main():
    parser=argparse.ArgumentParser();parser.add_argument("action",choices=("launch","worker","collect"))
    parser.add_argument("--run",type=Path,default=RUN);parser.add_argument("--basals",type=Path,default=BASALS)
    parser.add_argument("--topic",type=Path,default=TOPIC);parser.add_argument("--slot",type=int)
    args=parser.parse_args()
    if args.action=="launch":launch(args.run.resolve(),args.basals.resolve(),args.topic.resolve())
    elif args.action=="worker":worker(args.run.resolve(),args.slot)
    else:collect(args.run.resolve())


if __name__=="__main__":main()
