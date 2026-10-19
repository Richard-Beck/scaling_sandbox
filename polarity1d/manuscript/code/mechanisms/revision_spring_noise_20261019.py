"""Frozen SLURM pilot/main noisy Spring rings in native units.

The imposed field is stationary in material coordinates, with declared basal
physical covariance width. Dynamic circumference changes its physical width.
No state clipping, realization normalization or cell velocity law is added.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
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
from scipy.signal import find_peaks
from spring_ring_model_20261019 import Closure

MANUSCRIPT = Path(__file__).resolve().parents[2]
PYTHON = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"


def now(): return datetime.now(timezone.utc).isoformat()
def read(path): return json.loads(Path(path).read_text())
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def array_sha(x): return hashlib.sha256(np.asarray(x, dtype="<f8", order="C").tobytes()).hexdigest()


def write(path, data):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    tmp.replace(p)


class MaterialNoise:
    def __init__(self, n, basal_C, end, seed, mean, kind="noise", knot=.5,
                 coupling=.25, tau=5., ell=.15, sigma=1., bridge_base=None):
        self.n, self.basal_C, self.end = int(n), float(basal_C), float(end)
        self.mean, self.kind, self.knot = float(mean), kind, float(knot)
        self.tau, self.ell, self.sigma = float(tau), float(ell), float(sigma)
        self.seed, self.coupling = int(seed), float(coupling)
        self.bridge_base = bridge_base
        self.preserved_base_sha = None
        stride = int(round(knot / coupling))
        assert stride >= 1 and np.isclose(stride * coupling, knot)
        dx = basal_C / n
        s = np.arange(n) * dx
        images = int(np.ceil(8 * ell / basal_C)) + 1
        kernel = sum(np.exp(-((s + j * basal_C) / ell)**2) for j in range(-images, images + 1))
        self.kernel = kernel / np.linalg.norm(kernel)
        self.covariance = np.fft.irfft(abs(np.fft.rfft(self.kernel))**2, n=n)
        assert abs(self.covariance[0] - 1) < 1e-12
        self.rho = float(np.exp(-knot / tau))
        self.knots = None
        self.latent_sha = None
        if kind == "noise":
            rng = np.random.Generator(np.random.PCG64(seed))
            generating_step = coupling if bridge_base is None else float(bridge_base)
            if bridge_base is not None:
                assert np.isclose(generating_step, 2*coupling)
            count = int(round(end / generating_step)) + 1
            white = rng.standard_normal((count, n))
            eta = np.fft.irfft(np.fft.rfft(white, axis=1) * np.fft.rfft(self.kernel), n=n, axis=1)
            fine = np.empty_like(eta)
            fine[0] = eta[0]
            rho = np.exp(-generating_step / tau)
            scale = np.sqrt(-np.expm1(-2 * generating_step / tau))
            for i in range(1, count): fine[i] = rho * fine[i-1] + scale * eta[i]
            if bridge_base is not None:
                self.preserved_base_sha = array_sha(fine)
                bridge_rng = np.random.default_rng(np.random.SeedSequence([seed, 20261019, 125]))
                bridge_white = bridge_rng.standard_normal((count-1,n))
                innovation = np.fft.irfft(np.fft.rfft(bridge_white,axis=1)*np.fft.rfft(self.kernel),n=n,axis=1)
                half_rho = np.exp(-coupling/tau)
                bridge = np.empty((2*count-1,n))
                bridge[::2] = fine
                bridge[1::2] = half_rho*(fine[:-1]+fine[1:])/(1+half_rho**2) + np.sqrt((1-half_rho**2)/(1+half_rho**2))*innovation
                fine=bridge
            self.latent_sha = array_sha(fine)
            self.knots = fine[::stride].copy()
        else:
            assert kind in ("mean", "basal")
        self.constant = np.full(n, 0. if kind == "basal" else mean)

    def gaussian(self, t):
        index = min(int(max(0., t) / self.knot), len(self.knots) - 2)
        a = min(max(float(t) / self.knot - index, 0.), 1.)
        variance = (1-a)**2 + a*a + 2*a*(1-a)*self.rho
        return ((1-a)*self.knots[index] + a*self.knots[index+1]) / np.sqrt(variance)

    def __call__(self, t):
        if self.knots is None: return self.constant
        return self.mean * np.exp(self.sigma * self.gaussian(t) - self.sigma**2 / 2)

    def audit(self):
        result = dict(kind=self.kind, seed=self.seed, mean=self.mean, sigma=self.sigma,
            basal_C=self.basal_C, basal_physical_ell=self.ell, tau=self.tau,
            knot=self.knot, coupling_knot=self.coupling, kernel_sha256=array_sha(self.kernel),
            latent_fine_knots_sha256=self.latent_sha,
            bridge_base_interval=self.bridge_base, preserved_base_knots_sha256=self.preserved_base_sha,
            covariance_zero=float(self.covariance[0]), spatial_zero_mode_retained=True,
            geometry="Material field; instantaneous physical ell(t)=basal_ell*C(t)/Cbasal",
            interpolation="Linear Gaussian knots divided by analytic marginal SD; lognormal input, no realized normalization")
        if self.knots is not None:
            z = self.knots
            I = self.mean * np.exp(self.sigma*z - self.sigma**2/2)
            shift = int(round(self.ell / (self.basal_C / self.n)))
            tshift = int(round(self.tau/self.knot))
            result.update(gaussian_knots_sha256=array_sha(z),
                empirical_gaussian_mean=float(z.mean()), empirical_gaussian_second_moment=float(np.mean(z*z)),
                empirical_input_mean_over_mu=float(I.mean()/self.mean),
                empirical_input_variance_over_mu2=float(I.var()/self.mean**2),
                expected_input_variance_over_mu2=float(np.expm1(self.sigma**2)),
                gaussian_covariance_at_ell=float(np.mean(z*np.roll(z, shift, axis=1))),
                expected_covariance_at_ell=float(self.covariance[shift]),
                gaussian_covariance_at_tau=float(np.mean(z[:-tshift]*z[tshift:])),
                expected_covariance_at_tau=float(np.exp(-tshift*self.knot/self.tau)),
                theoretical_integrated_input_variance_over_mu2=float(self.basal_C**2 / self.n * np.sum(np.expm1(self.sigma**2*self.covariance))),
                expected_first_angular_input_second_moment_over_mu2=float(self.basal_C**2 / self.n * np.sum(np.expm1(self.sigma**2*self.covariance)*np.cos(2*np.pi*np.arange(self.n)/self.n))))
        return result


def diagnostics(u, q):
    n = len(q)
    out = {k: [] for k in ("contrast", "current_half_contrast_regions", "peaks_prominence_ge_005", "dominant_region_fraction", "dominant_region_width_fraction")}
    previous = None
    for field in u:
        contrast = float(np.ptp(field))
        mask = field >= field.min() + .5*contrast
        components = []
        if contrast > 1e-12:
            if mask.all(): components = [np.arange(n)]
            else:
                for start in np.flatnonzero(mask & ~np.roll(mask, 1)):
                    component, index = [], int(start)
                    while mask[index]:
                        component.append(index)
                        index = (index+1) % n
                    components.append(np.asarray(component))
        records = []
        for cells in components:
            excess = field[cells] - field.min()
            angle = np.angle(np.sum(excess*np.exp(2j*np.pi*q[cells]))) / (2*np.pi) % 1
            records.append((excess.sum(), angle, len(cells)/n))
        if records:
            best = max(records, key=lambda r: r[0])
            position, width = float(best[1]), float(best[2])
        else:
            position, width = (float(q[np.argmax(field)]) if previous is None else previous), 0.
        pk, _ = find_peaks(np.tile(field, 3), prominence=.05)
        peaks = int(np.sum((pk >= n) & (pk < 2*n)))
        for k, value in zip(out, (contrast, len(components), peaks, position, width)):
            out[k].append(value)
        previous = position
    return {k: np.asarray(v) for k, v in out.items()}


def freeze(run, mode, inputs, knot=.5, coupling=.25, bridge_base=None):
    run = Path(run).resolve()
    if run.exists(): raise ValueError("Fresh run directory required")
    for name in ("inputs", "frozen_sources", "workers", "logs", "results"):
        (run/name).mkdir(parents=True)
    paths = [Path(__file__), Path(__file__).with_name("spring_ring_model_20261019.py"),
             Path(__file__).with_name("spring_ring_closure_20261017.py"), MANUSCRIPT/"code/core/models/spring.py"]
    if mode=="main":
        paths += [Path(__file__).with_name("revision_spring_noise_analysis_20261019.py"),
                  Path(__file__).with_name("vision_noise_analysis_stable_v2_20261017.py")]
    hashes = {}
    for path in paths:
        relative = path.relative_to(MANUSCRIPT)
        target = run/"frozen_sources"/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        hashes[str(relative)] = sha(target)
    inputs = json.loads(Path(inputs).read_text())
    basals = []
    for record in inputs:
        original = Path(record["path"])
        if record.get("sha256"): assert sha(original) == record["sha256"]
        with np.load(original) as z:
            state = z["state"].copy()
        C0, n = float(record["resting_circumference"]), (len(state)-1)//2
        closure = Closure(C0, n)
        closure.fields(state)
        assert np.max(abs(closure.rhs(state, 0)) / np.r_[np.full(2*n, 2.), C0])*1000 < 1e-4
        assert np.ptp(state[:n]) >= .05
        target = run/"inputs"/f"C{C0:g}.npz"
        np.savez_compressed(target, state=state, material_positions=closure.material_positions)
        basals.append(dict(resting_circumference=C0, basal_circumference=float(state[-1]), cells=n,
                           path=str(target.relative_to(run)), sha256=sha(target), original=str(original), original_sha256=sha(original)))
    if mode == "pilot": assert len(basals) == 1 and basals[0]["resting_circumference"] == 2
    else: assert len(basals) == 3
    jobs = []
    for basal in basals:
        for f in (.001, .01):
            for replicate in range(2 if mode == "pilot" else 12):
                jobs.append(dict(kind="noise", fraction=f, replicate=replicate, seed=202610190400 + len(jobs), variant="standard", **basal))
            jobs.append(dict(kind="mean", fraction=f, replicate=0, seed=0, variant="standard", **basal))
        jobs.append(dict(kind="basal", fraction=0., replicate=0, seed=0, variant="standard", **basal))
    if mode == "pilot":
        base = next(j for j in jobs if j["kind"] == "noise" and j["fraction"] == .01)
        jobs.append(dict(base, variant="tight_ode"))
        for candidate in ([j for j in jobs if j["kind"] == "noise" and j["variant"]=="standard"] if bridge_base is not None else [base]):
            jobs.append(dict(candidate, variant="fine_knots"))
    for i, job in enumerate(jobs): job["slot"] = i
    manifest = dict(generated_at=now(), mode=mode, horizon=400. if mode == "pilot" else 3200.,
        maxlag=100. if mode == "pilot" else 800., scalar_interval=.5, profile_interval=4.,
        noise=dict(ell=.15, tau=5., knot=float(knot), coupling=float(coupling), bridge_base=bridge_base, sigma=1., amplitude_reference=4.,
                   rationale="ell=.15 resolves approximately15 native .01 cells and is shorter than preparation support.2; tau5 is one tenth of global circumference relaxation50. No cross-model biological calibration.",
                   physical_law="stationary material field with basalphysical covarianceell.15; actualell scalesC(t)/Cbasal",
                   kernel="Periodized Gaussian exp(-(s+j*Cbasal)^2/ell^2), discreteL2 normalized. Circular convolution of iid Gaussian white noise; exact discrete covariance=kernel autocorrelation; zero mode retained.",
                   temporal_law="Exact stationary OU knots at declared coupling; declared coarse knot subsample; analytic variance-normalized interpolation; lognormalexp(Z-.5)."),
        solver=dict(method="BDF", rtol=1e-8, atol=1e-10, max_step=float(knot)), source_sha256=hashes, basals=basals, jobs=jobs,
        moment="M0=C(t)*mean(u); M1=C(t)*mean(u exp(i2piq)); P=M1/M1(0); q1=|M1|/M0",
        site_displacement="wrapped material deltaq; primary referencebasalarc=Cbasal*deltaq; midpointC*deltaq sensitivity; deltaC separately. No ancestry or material travel inference.",
        uncertainty="whole independent seeds, pooled endpoint products/energies or sums/counts; no time-origin independence",
        guardrails=["No tuning to force memory loss", "No cell movement law", "No class of abrupt reorientations", "Native stretching/dilution/diffusion retained", "Mean-force global extension closure has no local force propagation distance"],
        resources=dict(qos="small", cpus=1, memory="3G", time="01:00:00" if mode=="main" else "00:30:00"))
    if mode=="main":
        metadata=MANUSCRIPT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/supporting/metadata"
        evidence=[metadata/"spring_noise_refinement_plan.json",metadata/"spring_noise_bridge_audit.json",
                  MANUSCRIPT/"data/runs/revision_ring_spring_20261019/results/root_independent_preparation_holdout.json"]
        manifest_evidence=[]
        for source in evidence:
            target=run/"protocol"/source.name;target.parent.mkdir(exist_ok=True);shutil.copy2(source,target)
            manifest_evidence.append(dict(original=str(source),copy=str(target.relative_to(run)),sha256=sha(target)))
        assert read(metadata/"spring_noise_bridge_audit.json")["all_estimator_gates_pass"]
        manifest["scientific_evidence"]=manifest_evidence
        manifest["analysis_plan"]=dict(windows=["full","first_half","last_half"],maxlag=800.,bootstrap_draws=2000,
            representative="Within-condition IQRscaled medoid of finalReP,minR,meanR,multiregionfraction; L1 to median; unit zeroIQR fallback; slot tie-break",
            raw_signed_restricted_integrals=True,exponential_fits=False,adaptation_exclusion=False,
            interpolation_resolution=".125 selected after nested .25/.125 estimator sensitivity checks; finite interpolated field, not exact continuous OU")
    write(run/"manifest.json", manifest)
    return manifest


def submit(run):
    run = Path(run).resolve()
    assert not (run/"submission.json").exists()
    m = read(run/"manifest.json")
    command = shlex.join(["env", "PYTHONDONTWRITEBYTECODE=1", "OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", PYTHON,
        str(run/"frozen_sources/code/mechanisms"/Path(__file__).name), "worker", "--run", str(run)])
    args = ["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=3G", "--time="+m["resources"]["time"],
        "--array=0-"+str(len(m["jobs"])-1), "--job-name=spring-noise-"+m["mode"],
        "--output="+str(run/"logs/%A_%a.out"), "--error="+str(run/"logs/%A_%a.err"), "--wrap", "exec "+command]
    job = subprocess.check_output(args, text=True).strip().split(";")[0]
    result = dict(submitted_at=now(), job_id=job, command=args, manifest_sha256=sha(run/"manifest.json"))
    write(run/"submission.json", result)
    return result


def worker(run, slot):
    assert os.environ.get("SLURM_JOB_ID"), "Integrations require SLURM"
    run = Path(run).resolve()
    m = read(run/"manifest.json")
    for relative, digest in m["source_sha256"].items(): assert sha(run/"frozen_sources"/relative) == digest
    j = m["jobs"][slot]
    assert sha(run/j["path"]) == j["sha256"]
    directory = run/"workers"/f"{slot:04d}"
    directory.mkdir(exist_ok=False)
    status = dict(started_at=now(), producer="running", complete=False, status="running", job=j,
                  slurm_job_id=os.environ["SLURM_JOB_ID"], manifest_sha256=sha(run/"manifest.json"))
    write(directory/"status.json", status)
    start = time.monotonic()
    try:
        with np.load(run/j["path"]) as z: initial = z["state"].copy()
        C0, n, end = j["resting_circumference"], j["cells"], m["horizon"]
        closure = Closure(C0, n)
        noise = MaterialNoise(n, initial[-1], end, j["seed"], 4*j["fraction"], kind=j["kind"],
                              knot=m["noise"]["knot"]/2 if j["variant"] == "fine_knots" else m["noise"]["knot"],
                              coupling=m["noise"]["coupling"], tau=m["noise"]["tau"], ell=m["noise"]["ell"], sigma=m["noise"]["sigma"], bridge_base=m["noise"].get("bridge_base"))
        audit = noise.audit()
        write(directory/"noise_audit.json", audit)
        t = np.arange(0., end + .25, .5)
        rtol = 1e-10 if j["variant"] == "tight_ode" else 1e-8
        atol = 1e-12 if j["variant"] == "tight_ode" else 1e-10
        solution = solve_ivp(lambda t,y: closure.rhs(y, noise(t)), (0.,end), initial,
            method="BDF", t_eval=t, rtol=rtol, atol=atol, max_step=min(m["solver"]["max_step"], noise.knot),
            jac=lambda t,y: closure.jacobian(y, noise(t)))
        assert solution.success and len(solution.t) == len(t), solution.message
        y = solution.y.T
        u, v, C = y[:,:n], y[:,n:-1], y[:,-1]
        assert np.isfinite(y).all() and y[:,:-1].min() >= -5e-10
        amount = C * np.mean(u+v, axis=1)
        amount_error = float(np.max(abs(amount/amount[0]-1)))
        assert amount_error <= 1e-7
        q = closure.material_positions
        M0 = C*np.mean(u, axis=1)
        M1 = C*(u@np.exp(2j*np.pi*q))/n
        M2 = C*(u@np.exp(4j*np.pi*q))/n
        P = M1/M1[0]
        desc = diagnostics(u, q)
        inputmean = np.array([np.mean(noise(tt)) for tt in t])
        np.savez_compressed(directory/"trace.npz", time=t, circumference=C, P_real=P.real, P_imag=P.imag,
            M1_real=M1.real, M1_imag=M1.imag, M0=M0, q1=abs(M1)/M0, q2=abs(M2)/M0,
            amount=amount, local_input_mean=inputmean, material_positions=q, **desc)
        primary = np.flatnonzero(np.isclose(t/4, np.rint(t/4)))
        np.savez_compressed(directory/"primary.npz", time=t[primary], primary=u[primary], inactive=v[primary],
            circumference=C[primary], material_positions=q)
        full = np.flatnonzero(np.isclose(t/100, np.rint(t/100)))
        np.savez_compressed(directory/"full_states.npz", time=t[full], states=y[full], material_positions=q)
        np.savez_compressed(directory/"final_state.npz", state=y[-1], time=end)
        status.update(status="complete", native_time=end, scalar_samples=len(t), profile_samples=len(primary),
            all_finite=True, minimum_concentration=float(y[:,:-1].min()), max_amount_relative_error=amount_error,
            circumference_range=[float(C.min()), float(C.max())],
            physical_noise_ell_range=[float(.15*C.min()/C[0]),float(.15*C.max()/C[0])],
            final_R=float(abs(P[-1])), min_R=float(abs(P).min()),
            nfev=solution.nfev, njev=solution.njev, nlu=solution.nlu,
            output_sha256={p.name:sha(p) for p in directory.iterdir() if p.name != "status.json"})
    except Exception as e:
        status.update(status="error", error=f"{type(e).__name__}: {e}")
        raise
    finally:
        status.update(producer="terminal", complete=True, finished_at=now(), elapsed_seconds=time.monotonic()-start)
        write(directory/"status.json", status)
    return status


def collect(run):
    run = Path(run).resolve()
    m = read(run/"manifest.json")
    observed = subprocess.run(["squeue", "--noheader", "--jobs", read(run/"submission.json")["job_id"], "--format=%i|%T|%M|%R"], capture_output=True, text=True)
    rows = []
    for j in m["jobs"]:
        path = run/"workers"/f"{j['slot']:04d}"/"status.json"
        rows.append(read(path) if path.exists() else dict(job=j, status="not_started", complete=False))
    done = all(r.get("status") == "complete" and r.get("complete") for r in rows)
    result = dict(generated_at=now(), mode=m["mode"], complete=done, total=len(rows),
                  finished=sum(r.get("status")=="complete" for r in rows), statuses=rows,
                  scheduler=dict(returncode=observed.returncode, stdout=observed.stdout, stderr=observed.stderr))
    if done:
        for j, row in zip(m["jobs"], rows):
            folder = run/"workers"/f"{j['slot']:04d}"
            assert row["manifest_sha256"] == sha(run/"manifest.json")
            for name,digest in row["output_sha256"].items(): assert sha(folder/name) == digest
            assert row["max_amount_relative_error"] <= 1e-7
        if m["mode"] == "pilot":
            comparisons = []
            for job in [j for j in m["jobs"] if j["variant"]!="standard"]:
                variant=job["variant"]
                standard, = [j for j in m["jobs"] if j["kind"]=="noise" and j["seed"]==job["seed"] and j["variant"]=="standard"]
                with np.load(run/"workers"/f"{standard['slot']:04d}"/"trace.npz") as z:
                    p0=z["P_real"]+1j*z["P_imag"]; c0=z["circumference"].copy()
                audit0=read(run/"workers"/f"{standard['slot']:04d}"/"noise_audit.json")
                with np.load(run/"workers"/f"{standard['slot']:04d}"/"primary.npz") as z: u0=z["primary"].copy()
                folder=run/"workers"/f"{job['slot']:04d}"
                a=read(folder/"noise_audit.json")
                assert a["latent_fine_knots_sha256"]==audit0["latent_fine_knots_sha256"]
                if variant=="tight_ode": assert a["gaussian_knots_sha256"]==audit0["gaussian_knots_sha256"]
                with np.load(folder/"trace.npz") as z: p=z["P_real"]+1j*z["P_imag"]; c=z["circumference"].copy()
                with np.load(folder/"primary.npz") as z: u=z["primary"].copy()
                comparisons.append(dict(variant=variant, seed=job["seed"],fraction=job["fraction"],baseline_slot=standard["slot"],variant_slot=job["slot"],identical_latent_fine_knots=True,
                    max_P_difference=float(np.max(abs(p-p0))), max_C_difference=float(np.max(abs(c-c0))),
                    max_primary_difference=float(np.max(abs(u-u0))), max_primary_difference_over_initial_contrast=float(np.max(abs(u-u0))/np.ptp(u0[0]))))
            result["numerical_comparisons"]=comparisons
            result["validity_passed"]=all(c["max_P_difference"] < (1e-4 if c["variant"]=="tight_ode" else .05) for c in comparisons)
            result["numerical_gate_note"]="TightODE tolerance1e-4 onP; coupledinterpolation tolerance.05 onP; finitepilot sensitivity, not long-run convergence proof"
    write(run/"results/status.json",result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("freeze","submit","worker","collect"))
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--mode", choices=("pilot","main"))
    parser.add_argument("--inputs", type=Path)
    parser.add_argument("--slot", type=int)
    parser.add_argument("--noise-knot", type=float, default=.5)
    parser.add_argument("--noise-coupling", type=float, default=.25)
    parser.add_argument("--bridge-base", type=float)
    args=parser.parse_args()
    if args.command=="freeze": result=freeze(args.run,args.mode,args.inputs,args.noise_knot,args.noise_coupling,args.bridge_base)
    elif args.command=="submit": result=submit(args.run)
    elif args.command=="worker": result=worker(args.run,args.slot if args.slot is not None else int(os.environ["SLURM_ARRAY_TASK_ID"]))
    else: result=collect(args.run)
    print(json.dumps({k:v for k,v in result.items() if k not in ("jobs","statuses","source_sha256")},indent=2))


if __name__=="__main__": main()
