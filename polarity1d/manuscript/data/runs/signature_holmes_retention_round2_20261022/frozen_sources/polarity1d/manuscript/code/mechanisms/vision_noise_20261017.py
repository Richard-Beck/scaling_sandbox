"""Config-driven noisy native chemical rings for the vision-aligned checkpoint.

No native model is modified. ``describe`` and ``checks`` do not integrate any
model; ``freeze`` requires an explicit configuration and never submits work.
Only ``worker`` integrates, and it requires a SLURM allocation. The stochastic
input is a declared continuous interpolation of exact stationary OU knots,
not an exact continuous-time OU sample between those knots.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import math
import os
import shlex
import shutil
import signal
import subprocess
import sys
import time

for _thread_key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_thread_key] = "1"
sys.dont_write_bytecode = True

import numpy as np
import scipy
from scipy.integrate import BDF
from scipy.signal import find_peaks
import ring_distance_prototype_20261017 as ring

ROOT = Path(__file__).resolve().parents[2]
PYTHON = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"
BASE = ROOT / "data/runs/chemical_ring_functional_v2_20261017"
DEFAULT = ROOT / "data/runs/vision_noise_pilot_20261017"
TOPIC = ROOT / "figures/length_mechanisms_20261008/dev/vision_noise_20261017"
CUE_PEAKS = {"wave_pinning": .1, "goryachev": 6.4, "holmes_model3": 1.}
STOP = False


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_sha(array):
    """Hash canonical little-endian float64 bytes, not ZIP metadata."""
    array = np.asarray(array, dtype="<f8", order="C")
    return hashlib.sha256(array.tobytes(order="C")).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def json_value(value):
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    return value


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(json_value(value), indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def sample_times(end, interval):
    if not (np.isfinite(end) and np.isfinite(interval) and end > 0 and interval > 0):
        raise ValueError("Need finite positive horizon and output interval")
    result = np.arange(0., end, interval)
    return np.r_[result, float(end)]


def spatial_kernel(grid, ell):
    """Periodized Gaussian filter with Gaussian covariance length ``ell``.

    The filter Gaussian has SD ell/sqrt(2); its autocorrelation has SD ell.
    Discrete L2 normalization gives exactly unit marginal Gaussian variance.
    The zero spatial Fourier mode is retained: no samplewise mean subtraction.
    """
    if not np.isfinite(ell) or ell <= 0:
        raise ValueError("Need a positive physical spatial correlation length")
    offsets = np.arange(grid.cells) * grid.dx
    images = int(math.ceil(8 * ell / grid.circumference)) + 1
    kernel = sum(np.exp(-((offsets + image * grid.circumference) / ell) ** 2)
                 for image in range(-images, images + 1))
    return kernel / np.linalg.norm(kernel)


def smooth_white(white, kernel):
    return np.fft.irfft(np.fft.rfft(white, axis=-1) * np.fft.rfft(kernel),
                        n=len(kernel), axis=-1)


class NoiseInput:
    """Reproducible stationary Gaussian knots and positive native-rate input.

    Z_k = rho Z_(k-1) + sqrt(1-rho**2) eta_k, including a stationary Z_0.
    Linear interpolation at fraction a is divided by its analytic standard
    deviation. Consequently Z(x,t) is marginally N(0,1) at every t, and
    E[mu exp(sigma Z - sigma**2/2)] = mu without realized-field normalization.
    Its inter-knot temporal covariance is defined by this interpolation.
    """
    def __init__(self, grid, end, seed, mean, sigma=1., ell=1.5, tau=50., knot=5.,
                 coupling_knot=2.5, kind="noise"):
        if kind not in ("noise", "mean", "basal"):
            raise ValueError("kind must be noise, mean or basal")
        if any(not np.isfinite(v) for v in (end, mean, sigma, ell, tau, knot,coupling_knot)):
            raise ValueError("Noise parameters must be finite")
        if end <= 0 or mean < 0 or sigma < 0 or tau <= 0 or knot <= 0 or coupling_knot<=0:
            raise ValueError("Invalid noise parameters")
        stride=int(round(knot/coupling_knot))
        if stride<1 or not np.isclose(stride*coupling_knot,knot,rtol=0,atol=1e-12):
            raise ValueError("knot must be an integer multiple of coupling_knot")
        self.grid, self.end, self.seed = grid, float(end), int(seed)
        self.mean, self.sigma, self.ell = float(mean), float(sigma), float(ell)
        self.tau, self.knot, self.kind = float(tau), float(knot), kind
        self.coupling_knot=float(coupling_knot)
        self.rho = float(np.exp(-self.knot / self.tau))
        self.kernel = spatial_kernel(grid, self.ell)
        self.knots = None
        if kind == "noise":
            generator = np.random.Generator(np.random.PCG64(self.seed))
            count = int(np.ceil(self.end / self.knot))*stride + 1
            white = generator.standard_normal((count, grid.cells))
            innovations = smooth_white(white, self.kernel)
            fine = np.empty_like(innovations)
            fine[0] = innovations[0]
            base_rho=float(np.exp(-self.coupling_knot/self.tau))
            scale = float(np.sqrt(-np.expm1(-2 * self.coupling_knot / self.tau)))
            for index in range(1, count):
                fine[index] = base_rho * fine[index - 1] + scale * innovations[index]
            self.coupling_knots_sha256=array_sha(fine)
            self.knots=fine[::stride].copy()
        self.constant = np.full(grid.cells, 0. if kind == "basal" else self.mean)

    def gaussian(self, t):
        if self.knots is None:
            raise ValueError("Deterministic controls have no Gaussian realization")
        t = min(max(float(t), 0.), self.end)
        index = min(int(t / self.knot), len(self.knots) - 2)
        fraction = t / self.knot - index
        variance = ((1 - fraction) ** 2 + fraction ** 2
                    + 2 * fraction * (1 - fraction) * self.rho)
        return ((1 - fraction) * self.knots[index] + fraction * self.knots[index + 1]) / np.sqrt(variance)

    def __call__(self, t):
        if self.kind != "noise":
            return self.constant
        return self.mean * np.exp(self.sigma * self.gaussian(t) - .5 * self.sigma ** 2)

    def temporal_covariance(self,first,second):
        def coefficients(t):
            index=int(t/self.knot);a=t/self.knot-index
            variance=(1-a)**2+a*a+2*a*(1-a)*self.rho
            return ((index,(1-a)/np.sqrt(variance)),(index+1,a/np.sqrt(variance)))
        return float(sum(a*b*self.rho**abs(i-j) for i,a in coefficients(first) for j,b in coefficients(second)))

    def audit(self):
        covariance = np.fft.irfft(np.abs(np.fft.rfft(self.kernel)) ** 2,
                                  n=self.grid.cells)
        result = dict(kind=self.kind, seed=self.seed, generator="numpy.Generator(PCG64)",
            gaussian_kernel_sha256=array_sha(self.kernel), theoretical_variance=float(covariance[0]),
            physical_covariance_length=self.ell, covariance_definition="Gaussian SD ell; periodized discrete convolution",
            retained_spatial_zero_mode=True, knot_interval=self.knot, tau=self.tau, rho=self.rho,
            coupling_knot=self.coupling_knot,
            interpolation="linear Gaussian interpolation / analytic pointwise standard deviation",
            temporal_covariance_examples=[dict(start_phase=phase,lag=lag,
                gaussian_covariance=self.temporal_covariance(phase*self.knot,phase*self.knot+lag))
                for phase in (0.,.25,.5,.75) for lag in (0.,self.knot/2,self.knot,self.tau)],
            theoretical_input_mean=self.mean if self.kind != "basal" else 0.,
            theoretical_input_variance=(self.mean ** 2 * np.expm1(self.sigma ** 2)) if self.kind == "noise" else 0.)
        if self.knots is None:
            return result
        z = self.knots
        knot_input = self.mean * np.exp(self.sigma * z - .5 * self.sigma ** 2)
        spatial_shift = min(int(round(self.ell / self.grid.dx)), self.grid.cells // 2)
        temporal_shift = max(1, int(round(self.tau / self.knot)))
        midpoint = np.stack([self.gaussian((index + .5) * self.knot)
                             for index in range(len(z) - 1)])
        result.update(gaussian_knots_sha256=array_sha(z), knot_shape=list(z.shape),
            coupling_knots_sha256=self.coupling_knots_sha256,
            empirical_gaussian_mean=float(z.mean()), empirical_gaussian_second_moment=float(np.mean(z*z)),
            empirical_midpoint_mean=float(midpoint.mean()), empirical_midpoint_second_moment=float(np.mean(midpoint*midpoint)),
            empirical_input_mean=float(knot_input.mean()), empirical_input_second_moment=float(np.mean(knot_input*knot_input)),
            empirical_input_minimum=float(knot_input.min()), empirical_input_maximum=float(knot_input.max()),
            empirical_spatial_covariance_at_ell=float(np.mean(z*np.roll(z,spatial_shift,axis=1))),
            theoretical_spatial_covariance_at_ell=float(covariance[spatial_shift]),
            empirical_temporal_covariance_at_tau=(float(np.mean(z[:-temporal_shift]*z[temporal_shift:]))
                if temporal_shift < len(z) else None),
            theoretical_temporal_covariance_at_tau=float(self.rho ** temporal_shift),
            empirical_spatial_mean_sd=float(np.std(z.mean(axis=1))),
            finite_realization_warning="Empirical local moments/correlations fluctuate; they are diagnostics, not rescaling targets.")
        input_covariance=self.mean**2*np.expm1(self.sigma**2*covariance)
        raw_modes=self.mode_records()
        result.update(
            empirical_input_spatial_covariance_at_ell=float(np.mean((knot_input-self.mean)*np.roll(knot_input-self.mean,spatial_shift,axis=1))),
            theoretical_input_spatial_covariance_at_ell=float(input_covariance[spatial_shift]),
            empirical_input_temporal_covariance_at_tau=(float(np.mean((knot_input[:-temporal_shift]-self.mean)*(knot_input[temporal_shift:]-self.mean)))
                if temporal_shift<len(z) else None),
            theoretical_input_temporal_covariance_at_tau=float(self.mean**2*np.expm1(self.sigma**2*self.rho**temporal_shift)),
            empirical_integrated_input_mean=float(np.mean(raw_modes["input_M0"])),
            empirical_integrated_input_variance=float(np.var(raw_modes["input_M0"])),
            theoretical_integrated_input_mean=float(self.mean*self.grid.circumference),
            theoretical_integrated_input_variance=float(self.grid.circumference*self.grid.dx*input_covariance.sum()),
            empirical_first_angular_input_mode_second_moment=float(np.mean(raw_modes["input_M1_real"]**2+raw_modes["input_M1_imag"]**2)),
            theoretical_first_angular_input_mode_second_moment=float(self.grid.circumference*self.grid.dx*np.sum(input_covariance*np.cos(2*np.pi*np.arange(self.grid.cells)/self.grid.cells))))
        return result

    def mode_records(self):
        if self.knots is None:
            times=sample_times(self.end,self.knot)
            inputs=np.repeat(self.constant[None,:],len(times),axis=0)
        else:
            times=np.arange(len(self.knots))*self.knot
            inputs=self.mean*np.exp(self.sigma*self.knots-.5*self.sigma**2)
        first=inputs@np.exp(2j*np.pi*self.grid.x/self.grid.circumference)*self.grid.dx
        return dict(time=times,input_M0=inputs.sum(axis=1)*self.grid.dx,
            input_M1_real=first.real,input_M1_imag=first.imag,
            local_mean=inputs.mean(axis=1),local_variance=inputs.var(axis=1),
            local_minimum=inputs.min(axis=1),local_maximum=inputs.max(axis=1))


def retained_basals():
    manifest_path = BASE / "manifest.json"
    manifest = read(manifest_path)
    records = []
    for key in ring.KEYS:
        for circumference in (20, 40, 60):
            job = next(j for j in manifest["jobs"] if j["model"] == key and j["circumference"] == circumference)
            status_path = BASE / "workers" / f'{job["slot"]:02d}' / "status.json"
            status = read(status_path)
            if not status.get("complete") or status["status"] != "complete":
                raise ValueError(f"Unavailable retained basal: {key} C{circumference}")
            checkpoint = Path(status["checkpoint"])
            if sha(checkpoint) != status["checkpoint_sha256"]:
                raise ValueError(f"Changed checkpoint: {checkpoint}")
            with np.load(checkpoint) as archive:
                state, x = archive["state"], archive["x"]
                grid = ring.RingGrid(float(circumference), int(job["cells"]))
                np.testing.assert_array_equal(x, grid.x)
                assert state.shape == (len(ring.native_model(key).state_names), grid.cells)
                assert np.isfinite(state).all() and state.min() >= 0
                level, _ = ring.basal_core(state[0])
            records.append(dict(model=key, circumference=circumference, cells=grid.cells,
                dx=grid.dx, checkpoint=str(checkpoint), checkpoint_sha256=sha(checkpoint),
                status_path=str(status_path), status_sha256=sha(status_path),
                frozen_basal_half_contrast=level, model_peak_reference=CUE_PEAKS[key]))
    return records


def describe():
    return dict(generated_at=now(), basals=retained_basals(),
        proposed_noise=dict(ell=1.5, tau=50., knot=5., coupling_knot=2.5, sigma=1.,
            mean="native finite-cue peak times a declared model-specific fraction; not biologically calibrated",
            spatial_law="periodized Gaussian convolution, discrete L2-normalized, Gaussian covariance SD ell",
            temporal_law="stationary exact OU knots with variance-normalized continuous interpolation",
            local_stats="fixed within model across C; no samplewise mean/variance normalization",
            stationarity="stationary OU knots and unit Gaussian marginals from t=0; normalized interpolation has phase-dependent covariance; cell begins from the exact cue-free basal"),
        recommended_numerics=dict(rtol=1e-8, atol=1e-10, max_step=5.,
            audit="same frozen noise realization with tighter tolerances and max_step=1; no temporal stochastic-law change"),
        protocol_status="Proposed infrastructure only; no trajectory grid frozen or submitted by describe",
        configuration_schema=dict(noise=dict(ell=1.5,tau=50.,knot=5.,coupling_knot=2.5,sigma=1.),
            horizon=4000., output_interval=5., primary_interval=20., full_interval=200.,
            solver=dict(rtol=1e-8,atol=1e-10,max_step=5.),
            jobs=[dict(model="wave_pinning",circumference=20,kind="noise",mean_fraction=.001,seed=202610170001)]))


def freeze(run, config_path):
    run, config_path = Path(run).resolve(), Path(config_path).resolve()
    if run.exists():
        raise ValueError("Choose a fresh output run; existing runs are immutable")
    if not run.is_relative_to(ROOT / "data/runs"):
        raise ValueError("Run must be under manuscript/data/runs")
    config = read(config_path)
    if not config.get("jobs"):
        raise ValueError("Configuration must explicitly declare jobs")
    protocol_files=config.get("protocol_files",{})
    missing_labels=set(config.get("protocol_required_labels",[]))-set(protocol_files)
    if missing_labels:
        raise ValueError(f"Missing explicitly required scientific records: {sorted(missing_labels)}")
    for label,path in protocol_files.items():
        if Path(label).name!=label or label in (".","..") or not Path(path).is_file():
            raise ValueError(f"Invalid scientific protocol input: {label}: {path}")
    basals = retained_basals()
    for name in ("inputs", "frozen_sources", "workers", "logs", "results"):
        (run / name).mkdir(parents=True)
    shutil.copy2(config_path, run / "configuration.json")
    paths = [Path(__file__), Path(__file__).with_name("ring_distance_prototype_20261017.py"),
             Path(__file__).with_name("length_dev_output_layout.py")]
    for folder in (ROOT / "code/core/models", ROOT / "code/core/src/polarity1d"):
        paths.extend(sorted(folder.glob("*.py")))
    sources = {}
    for source in paths:
        relative = source.relative_to(ROOT)
        target = run / "frozen_sources" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        sources[str(relative)] = sha(target)
    external, jobs, inputs = {}, [], {}
    for slot, supplied in enumerate(config["jobs"]):
        job = dict(supplied)
        basal = next(b for b in basals if b["model"] == job["model"] and b["circumference"] == job["circumference"])
        identity = f'{basal["model"]}_C{basal["circumference"]}'
        target = run / "inputs" / f"{identity}.npz"
        if identity not in inputs:
            shutil.copy2(basal["checkpoint"], target)
            inputs[identity] = dict(path=str(target.relative_to(run)),sha256=sha(target),source=basal)
        external[basal["checkpoint"]] = basal["checkpoint_sha256"]
        external[basal["status_path"]] = basal["status_sha256"]
        job.update(slot=slot, input=inputs[identity]["path"], input_sha256=inputs[identity]["sha256"],
            cells=basal["cells"], dx=basal["dx"], mean=CUE_PEAKS[job["model"]]*float(job["mean_fraction"]))
        if job["kind"] not in ("noise","mean","basal") or job["mean"] < 0:
            raise ValueError("Invalid job kind or mean fraction")
        job["seed"] = int(job["seed"])
        job["condition"]="zero" if job["kind"]=="basal" else job["kind"]
        jobs.append(job)
    manifest = dict(created_at=now(),experiment="Vision-aligned native ring noise ensemble",
        configuration_sha256=sha(run/"configuration.json"),jobs=jobs,inputs=inputs,
        noise=config["noise"],horizon=float(config["horizon"]),
        output_interval=float(config.get("output_interval",5.)),primary_interval=float(config.get("primary_interval",20.)),full_interval=float(config.get("full_interval",200.)),
        solver=config.get("solver",dict(rtol=1e-8,atol=1e-10,max_step=5.)),
        diagnostics=dict(allowed_negative_roundoff=5e-10,pool_relative_tolerance=1e-7,
            minimum_contrast_native=.05,peak_prominence_native=.05),
        source_sha256=sources,external_sha256=external,original_configuration=str(config_path),
        time_budget=dict(slurm_seconds=3600,worker_soft_seconds=3300),
        interpretation="Activity polarity vectors and local patterns, not cell velocity or migration trajectories",
        environment=dict(executable=PYTHON,python=sys.version,numpy=np.__version__,scipy=scipy.__version__))
    if protocol_files:
        copied={}
        for label,path in protocol_files.items():
            source=Path(path).resolve();target=run/"protocol"/(label+source.suffix)
            target.parent.mkdir(exist_ok=True);shutil.copy2(source,target)
            copied[label]=dict(source=str(source),source_sha256=sha(source),
                copy=str(target.relative_to(run)),copy_sha256=sha(target))
        record=dict(frozen_at=now(),configuration_sha256=manifest["configuration_sha256"],
            analysis=config.get("analysis",{}),files=copied,
            scope="Scientific plan, prospective analysis/selection rules, reviews and authorization frozen before producer launch")
        write(run/"plan_freeze.json",record)
        manifest["plan_freeze_sha256"]=sha(run/"plan_freeze.json")
    if manifest["solver"]["max_step"]>manifest["noise"]["knot"]:
        raise ValueError("Maximum integration step must not exceed the noise knot spacing")
    for job in jobs:
        grid=ring.RingGrid(job["circumference"],job["cells"])
        # Constructing a short realization validates parameter semantics, not dynamics.
        NoiseInput(grid,min(manifest["horizon"],5.),job["seed"],job["mean"],kind=job["kind"],**manifest["noise"])
    write(run/"manifest.json",manifest)
    write(run/"status.json",dict(generated_at=now(),complete=False,producer="frozen_not_submitted"))
    return manifest


def submit(run, slots=None, label="initial"):
    run=Path(run).resolve()
    manifest=read(run/"manifest.json")
    runner=run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    if slots is None:
        slots=f'0-{len(manifest["jobs"])-1}' if len(manifest["jobs"])>1 else "0"
    command="exec "+shlex.join(["env","PYTHONDONTWRITEBYTECODE=1","OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1","MKL_NUM_THREADS=1",PYTHON,str(runner),"worker","--run",str(run)])
    command+=' --slot "$SLURM_ARRAY_TASK_ID"'
    arguments=["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1","--mem=3G",
        "--time=01:00:00","--signal=B:USR1@60","--array="+slots,"--job-name=vision-noise",
        "--output="+str(run/"logs/%A_%a.out"),"--error="+str(run/"logs/%A_%a.err"),"--wrap",command]
    job_id=subprocess.check_output(arguments,text=True).strip().split(";")[0]
    record=dict(generated_at=now(),job_id=job_id,slots=slots,command=arguments,manifest_sha256=sha(run/"manifest.json"))
    write(run/f"launch_{label}.json",record)
    return record


def profile_record(grid, primary, basal, previous_position, previous_unwrapped):
    """Local descriptors retain ambiguities; dominant-region labels imply no ancestry."""
    background=float(primary.min()); contrast=float(np.ptp(primary))
    initial_level=float(basal.min()+.5*np.ptp(basal))
    current_level=background+.5*contrast
    dynamic_regions=ring.circular_components(primary>=current_level) if contrast>1e-12 else []
    fixed_regions=ring.circular_components(primary>=initial_level)
    peaks, properties=find_peaks(np.tile(primary,3),prominence=0.)
    keep=(peaks>=grid.cells)&(peaks<2*grid.cells)
    peak_cells=(peaks[keep]-grid.cells).astype(int)
    prominences=properties["prominences"][keep]
    records=[]
    for component in dynamic_regions:
        excess=primary[component]-background
        mass=float(np.sum(excess)*grid.dx)
        local_moment=np.sum(excess*np.exp(2j*np.pi*grid.x[component]/grid.circumference))
        position=float(np.angle(local_moment)%(2*np.pi)*grid.circumference/(2*np.pi))
        records.append(dict(cells=int(len(component)),width=float(len(component)*grid.dx),
            excess_amount=mass,center=position,peak_value=float(primary[component].max())))
    dominant=max(records,key=lambda r:r["excess_amount"]) if records else None
    position=dominant["center"] if dominant is not None else previous_position
    unwrapped=previous_unwrapped+float(grid.arc(position,previous_position))
    resolved=int(np.sum(prominences>=.05))
    scalar=dict(primary_minimum=background,primary_maximum=float(primary.max()),contrast=contrast,
        frozen_threshold_regions=len(fixed_regions),current_half_contrast_regions=len(dynamic_regions),
        peaks_prominence_ge_005=resolved,dominant_region_position=position,
        dominant_region_unwrapped_position=unwrapped,
        dominant_region_width=dominant["width"] if dominant else 0.,
        dominant_region_excess_amount=dominant["excess_amount"] if dominant else 0.,
        dominant_region_unique=bool(len(dynamic_regions)==1 and contrast>=.05),
        dominant_region_tracking_warning=bool(len(dynamic_regions)!=1 or resolved>1 or contrast<.05))
    detail=dict(regions=records,peak_cells=peak_cells.tolist(),
        peak_positions=grid.x[peak_cells].tolist(),peak_values=primary[peak_cells].tolist(),
        peak_prominences=prominences.tolist())
    return scalar,detail,position,unwrapped


def integrate(model,grid,initial,forcing,manifest,deadline,directory,status):
    times=sample_times(manifest["horizon"],manifest["output_interval"])
    primary_times=sample_times(manifest["horizon"],manifest["primary_interval"])
    full_times=sample_times(manifest["horizon"],manifest["full_interval"])
    pools=ring.pool_amounts(model,grid,initial)
    initial_peak_position=float(grid.x[int(np.argmax(initial[0]))])
    _,_,initial_position,_=profile_record(grid,initial[0],initial[0],initial_peak_position,initial_peak_position)
    phase=np.exp(2j*np.pi*grid.x/grid.circumference)*grid.dx
    phase2=np.exp(4j*np.pi*grid.x/grid.circumference)*grid.dx
    initial_moment=complex(initial[0]@phase)
    if abs(initial_moment)<=1e-12*float(initial[0].sum()*grid.dx):
        raise ValueError("Initial first moment is too small to define reference polarity")
    full_states=[initial.copy()]; collected_times=[]; primary=[initial[0].copy()]; scalars=[]; regions=[]
    last_position,last_unwrapped=initial_position,initial_position
    minimum=float(initial.min());worst_pool=0.;noise_input_min=math.inf;noise_input_max=0.
    def inspect_native(state):
        nonlocal minimum,worst_pool
        if not np.isfinite(state).all():
            raise FloatingPointError("Nonfinite native fields")
        minimum=min(minimum,float(state.min()))
        worst_pool=max(worst_pool,float(np.max(np.abs(ring.pool_amounts(model,grid,state)/pools-1))))
        if minimum < -manifest["diagnostics"]["allowed_negative_roundoff"] or worst_pool>manifest["diagnostics"]["pool_relative_tolerance"]:
            raise FloatingPointError("Native concentration or pool-drift tolerance exceeded")
    def record(t,state):
        nonlocal last_position,last_unwrapped,noise_input_min,noise_input_max
        inspect_native(state)
        moment=complex(state[0]@phase);moment2=complex(state[0]@phase2)
        p=moment/initial_moment
        amount=float(state[0].sum()*grid.dx)
        moment_denominator_defined=bool(amount>0)
        coherence=abs(moment)/amount if moment_denominator_defined else 0.
        local,detail,last_position,last_unwrapped=profile_record(grid,state[0],initial[0],last_position,last_unwrapped)
        stimulus=forcing(t)
        noise_input_min=min(noise_input_min,float(stimulus.min()));noise_input_max=max(noise_input_max,float(stimulus.max()))
        local.update(time=float(t),polarity_x=float(p.real),polarity_y=float(p.imag),
            polarity_magnitude=float(abs(p)),active_amount=amount,coherence=float(coherence),
            q1=float(coherence),q2=float(abs(moment2)/amount) if moment_denominator_defined else 0.,
            normalized_moment_denominator_defined=moment_denominator_defined,
            M1_real=moment.real,M1_imag=moment.imag,M2_real=moment2.real,M2_imag=moment2.imag,
            dominant_region_displacement_from_initial=float(grid.arc(last_position,initial_position)),
            dominant_region_unwrapped_displacement=float(last_unwrapped-initial_position),
            largest_peak_position=float(grid.x[int(np.argmax(state[0]))]),
            local_input_mean=float(stimulus.mean()),local_input_minimum=float(stimulus.min()),local_input_maximum=float(stimulus.max()))
        scalars.append(local);regions.append(dict(time=float(t),**detail));collected_times.append(float(t))
    record(0.,initial)
    solver=BDF(lambda t,y:ring.native_rhs(model,grid,y.reshape(initial.shape),forcing(t)).ravel(),
        0.,initial.ravel(),manifest["horizon"],jac_sparsity=ring.periodic_sparsity(len(model.state_names),grid.cells),
        **manifest["solver"])
    wanted=1;wanted_primary=1;wanted_full=1;last_heartbeat=time.monotonic();accepted_steps=0
    try:
        while solver.status=="running":
            if STOP or time.monotonic()>=deadline:
                raise TimeoutError("SLURM advance signal or worker soft deadline")
            message=solver.step()
            if solver.status=="failed":
                raise RuntimeError(message)
            state=solver.y.reshape(initial.shape)
            inspect_native(state)
            dense=solver.dense_output()
            while wanted<len(times) and times[wanted]<=solver.t:
                record(times[wanted],dense(times[wanted]).reshape(initial.shape));wanted+=1
            while wanted_primary<len(primary_times) and primary_times[wanted_primary]<=solver.t:
                sampled=dense(primary_times[wanted_primary]).reshape(initial.shape);inspect_native(sampled)
                primary.append(sampled[0].copy());wanted_primary+=1
            while wanted_full<len(full_times) and full_times[wanted_full]<=solver.t:
                sampled=dense(full_times[wanted_full]).reshape(initial.shape);inspect_native(sampled)
                full_states.append(sampled.copy());wanted_full+=1
            accepted_steps+=1
            if time.monotonic()-last_heartbeat>=30.:
                status.update(generated_at=now(),producer="running",native_time=float(solver.t),
                    accepted_steps=accepted_steps,minimum_concentration=minimum,max_pool_relative_error=worst_pool)
                write(directory/"status.json",status);last_heartbeat=time.monotonic()
        termination="complete";message=""
    except Exception as error:
        termination="timeout" if isinstance(error,TimeoutError) else "numerical_error"
        message=f"{type(error).__name__}: {error}"
    # Partial output and last accepted state remain available even on failure.
    np.savez_compressed(directory/"primary.npz",time=primary_times[:len(primary)],primary=np.asarray(primary),x=grid.x,circumference=np.asarray(grid.circumference))
    np.savez_compressed(directory/"full_states.npz",time=full_times[:len(full_states)],states=np.asarray(full_states),x=grid.x)
    np.savez_compressed(directory/"final_state.npz",state=solver.y.reshape(initial.shape),time=np.asarray(solver.t),x=grid.x)
    with (directory/"scalars.csv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(scalars[0]));writer.writeheader();writer.writerows(scalars)
    with (directory/"local_patterns.jsonl").open("w") as handle:
        for row in regions:
            handle.write(json.dumps(json_value(row),allow_nan=False)+"\n")
    trace={name:np.asarray([row[name] for row in scalars]) for name in scalars[0]}
    trace.update(P_real=trace["polarity_x"],P_imag=trace["polarity_y"],M0=trace["active_amount"],
        circumference=np.asarray(grid.circumference),x=grid.x)
    np.savez_compressed(directory/"trace.npz",**trace)
    return dict(status=termination,message=message,complete=termination=="complete",native_time=float(solver.t),
        samples=len(scalars),full_samples=len(full_states),minimum_concentration=minimum,
        max_pool_relative_error=worst_pool,nfev=int(solver.nfev),njev=int(solver.njev),nlu=int(solver.nlu),
        accepted_steps=accepted_steps,sampled_input_minimum=noise_input_min,sampled_input_maximum=noise_input_max,
        initial_moment_real=initial_moment.real,initial_moment_imag=initial_moment.imag,
        interpretation="Dominant-region displacement is descriptive and can change identity; flags and full profiles retain ambiguity.")


def worker(run,slot):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Native integrations require a SLURM allocation")
    signal.signal(signal.SIGUSR1,lambda signum,frame:globals().__setitem__("STOP",True))
    signal.signal(signal.SIGTERM,lambda signum,frame:globals().__setitem__("STOP",True))
    run=Path(run).resolve();manifest=read(run/"manifest.json")
    for path,digest in manifest["source_sha256"].items():
        if sha(run/"frozen_sources"/path)!=digest:
            raise ValueError(f"Changed frozen source: {path}")
    if sha(run/"configuration.json")!=manifest["configuration_sha256"]:
        raise ValueError("Changed frozen configuration")
    if "plan_freeze_sha256" in manifest:
        if sha(run/"plan_freeze.json")!=manifest["plan_freeze_sha256"]:
            raise ValueError("Changed frozen scientific plan record")
        for item in read(run/"plan_freeze.json")["files"].values():
            if sha(run/item["copy"])!=item["copy_sha256"]:
                raise ValueError(f'Changed frozen scientific protocol: {item["copy"]}')
    job=manifest["jobs"][slot]
    if sha(run/job["input"])!=job["input_sha256"]:
        raise ValueError("Changed copied native basal")
    directory=run/"workers"/f"{slot:04d}"
    directory.mkdir()
    status=dict(generated_at=now(),producer="running",status="starting",complete=False,
        job=job,slurm_job_id=os.environ["SLURM_JOB_ID"],manifest_sha256=sha(run/"manifest.json"))
    write(directory/"status.json",status)
    begin=time.monotonic()
    try:
        grid=ring.RingGrid(job["circumference"],job["cells"]);model=ring.native_model(job["model"])
        with np.load(run/job["input"]) as archive:
            initial=archive["state"].copy();np.testing.assert_array_equal(archive["x"],grid.x)
        forcing=NoiseInput(grid,manifest["horizon"],job["seed"],job["mean"],kind=job["kind"],**manifest["noise"])
        write(directory/"noise_audit.json",forcing.audit())
        modes=forcing.mode_records()
        with (directory/"noise_modes.csv").open("w",newline="") as handle:
            writer=csv.writer(handle);writer.writerow(list(modes));writer.writerows(zip(*modes.values()))
        result=integrate(model,grid,initial,forcing,manifest,
            begin+manifest["time_budget"]["worker_soft_seconds"],directory,status)
        status.update(result)
    except Exception as error:
        status.update(status="worker_error",message=f"{type(error).__name__}: {error}",complete=False)
    status.update(generated_at=now(),producer="terminal",elapsed_seconds=time.monotonic()-begin)
    status["output_sha256"]={str(p.relative_to(directory)):sha(p) for p in sorted(directory.iterdir()) if p.is_file() and p.name!="status.json"}
    write(directory/"status.json",status)
    if not status["complete"]:
        raise RuntimeError(status.get("message",status["status"]))
    return status


def checks():
    """Algebraic/statistical preflight only; no native trajectory integration."""
    results=[]
    rng=np.random.Generator(np.random.PCG64(20261017001))
    for circumference in (20,40,60):
        grid=ring.RingGrid(circumference,int(circumference/.05));kernel=spatial_kernel(grid,1.5)
        np.testing.assert_allclose(kernel@kernel,1.,rtol=0.,atol=5e-16)
        covariance=np.fft.irfft(np.abs(np.fft.rfft(kernel))**2,n=grid.cells)
        displacement=grid.distance(np.arange(grid.cells)*grid.dx,0.)
        images=sum(np.exp(-.5*((np.arange(grid.cells)*grid.dx+k*circumference)/1.5)**2) for k in range(-2,3))
        expected=images/images[0]
        covariance_error=float(np.max(np.abs(covariance-expected)))
        assert covariance_error<2e-14
        white=rng.standard_normal((32,grid.cells));smooth=smooth_white(white,kernel)
        rotation_errors=[]
        for shift in (1,17,grid.cells//2):
            error=float(np.max(np.abs(smooth_white(np.roll(white,shift,axis=1),kernel)-np.roll(smooth,shift,axis=1))))
            assert error<5e-15;rotation_errors.append(error)
        forcing=NoiseInput(grid,100.,12345,.002)
        duplicate=NoiseInput(grid,100.,12345,.002)
        np.testing.assert_array_equal(forcing.knots,duplicate.knots)
        refined=NoiseInput(grid,100.,12345,.002,knot=2.5)
        np.testing.assert_array_equal(forcing.knots,refined.knots[::2])
        assert np.isfinite(forcing(12.3)).all() and forcing(12.3).min()>0
        interpolation_errors=[]
        for index in (1,5,10):
            np.testing.assert_allclose(forcing.gaussian(index*5),forcing.knots[index],atol=2e-15)
            left=forcing.gaussian(index*5-1e-8);right=forcing.gaussian(index*5+1e-8)
            interpolation_errors.append(float(np.max(np.abs(left-right))))
        assert max(interpolation_errors)<1e-7
        rho=forcing.rho
        for a in (.1,.5,.9):
            variance=((1-a)**2+a*a+2*a*(1-a)*rho)
            coefficients=np.array([1-a,a])/np.sqrt(variance)
            np.testing.assert_allclose(coefficients@np.array([[1,rho],[rho,1]])@coefficients,1.,atol=5e-16)
        # Independent fields, not correlated knots, provide an empirical variance check.
        independent=smooth_white(rng.standard_normal((1024,grid.cells)),kernel)
        empirical_mean=float(independent.mean());empirical_second=float(np.mean(independent**2))
        for kind,expected_mean in (("mean",.002),("basal",0.)):
            control=NoiseInput(grid,100.,12345,.002,kind=kind)
            np.testing.assert_array_equal(control(33.),np.full(grid.cells,expected_mean))
        results.append(dict(circumference=circumference,passed=True,covariance_error=covariance_error,
            coupled_refinement_coarse_knots_identical=True,
            coupled_refinement_midpoint_gaussian_rms=float(np.sqrt(np.mean((forcing.gaussian(12.5)-refined.gaussian(12.5))**2))),
            translation_max_error=max(rotation_errors),continuity_max_difference=max(interpolation_errors),
            independent_empirical_mean=empirical_mean,independent_empirical_second_moment=empirical_second))
    basal_checks=[]
    for basal in retained_basals():
        grid=ring.RingGrid(basal["circumference"],basal["cells"])
        model=ring.native_model(basal["model"])
        with np.load(basal["checkpoint"]) as archive:state=archive["state"].copy()
        phase=np.exp(2j*np.pi*grid.x/grid.circumference)*grid.dx
        moment=complex(state[0]@phase);amount=float(state[0].sum()*grid.dx)
        assert abs(moment)>1e-12*amount
        np.testing.assert_allclose(moment/moment,1+0j,rtol=0,atol=2e-16)
        local,_,position,_=profile_record(grid,state[0],state[0],grid.x[int(np.argmax(state[0]))],grid.x[int(np.argmax(state[0]))])
        assert local["dominant_region_unique"]
        stimulus=np.full(grid.cells,CUE_PEAKS[basal["model"]]*.001)
        delta=model.reaction(0.,grid.x,state,stimulus)-model.reaction(0.,grid.x,state,np.zeros(grid.cells))
        expected=stimulus*state[1]/(7.5 if basal["model"]=="holmes_model3" else 1.)
        np.testing.assert_allclose(delta[0],expected,rtol=1e-10,atol=5e-15)
        np.testing.assert_allclose(delta[1],-expected,rtol=1e-10,atol=5e-15)
        poolrate=ring.pool_amounts(model,grid,ring.native_rhs(model,grid,state,stimulus))
        assert np.max(np.abs(poolrate))<1e-10
        basal_checks.append(dict(model=basal["model"],circumference=grid.circumference,
            basal_sha256=basal["checkpoint_sha256"],initial_moment_magnitude=abs(moment),
            initial_active_amount=amount,dominant_region_center=position,
            native_uniform_input_channel_max_error=float(np.max(np.abs(delta[0]-expected))),
            native_pool_rate_max=float(np.max(np.abs(poolrate)))))
    return dict(generated_at=now(),passed=True,tests=results,basal_checks=basal_checks,
        scope="Gaussian field algebra/statistics, reproducibility, continuity, covariance; no native model integration")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command",choices=("describe","checks","freeze","submit","worker"))
    parser.add_argument("--run",type=Path,default=DEFAULT)
    parser.add_argument("--config",type=Path)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--slot",type=int)
    parser.add_argument("--slots")
    parser.add_argument("--label",default="initial")
    args=parser.parse_args()
    if args.command=="describe":result=describe()
    elif args.command=="checks":result=checks()
    elif args.command=="freeze":
        if args.config is None:parser.error("freeze requires --config")
        result=freeze(args.run,args.config)
    elif args.command=="submit":result=submit(args.run,args.slots,args.label)
    else:
        if args.slot is None:parser.error("worker requires --slot")
        result=worker(args.run,args.slot)
    if args.output:write(args.output,result)
    print(json.dumps(json_value(result),indent=2,allow_nan=False))


if __name__=="__main__":
    main()
