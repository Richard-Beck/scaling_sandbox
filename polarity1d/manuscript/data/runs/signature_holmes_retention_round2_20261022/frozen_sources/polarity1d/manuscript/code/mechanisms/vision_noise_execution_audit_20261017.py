"""Reproducible archive verification and separated noisy-ring numerical audits."""
from pathlib import Path
import argparse
import importlib.util
import json
import os
import subprocess
import sys

for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key]="1"
sys.dont_write_bytecode=True
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/"data/runs/vision_noise_pilot_20261017"
AUDIT=ROOT/"data/runs/vision_noise_numerical_audit_20261017"


def load_frozen(run):
    path=run/"frozen_sources/code/mechanisms/vision_noise_20261017.py"
    sys.path.insert(0,str(path.parent))
    spec=importlib.util.spec_from_file_location("vision_noise_frozen_audit",path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def allocated(path):
    return int(subprocess.check_output(["du","-s","-B1",str(path)],text=True).split()[0])


def verify(run,v):
    manifest=v.read(run/"manifest.json")
    for path,digest in manifest["source_sha256"].items():
        assert v.sha(run/"frozen_sources"/path)==digest,path
    assert v.sha(run/"configuration.json")==manifest["configuration_sha256"]
    cases=[]
    for job in manifest["jobs"]:
        directory=run/"workers"/f'{job["slot"]:04d}'
        status=v.read(directory/"status.json")
        assert status["complete"] and status["status"]=="complete" and status["producer"]=="terminal",job
        assert status["native_time"]==manifest["horizon"]
        assert v.sha(run/job["input"])==job["input_sha256"]
        for path,digest in status["output_sha256"].items():
            assert v.sha(directory/path)==digest,path
        with np.load(run/job["input"]) as archive:initial=archive["state"].copy()
        grid=v.ring.RingGrid(job["circumference"],job["cells"]);model=v.ring.native_model(job["model"])
        with np.load(directory/"trace.npz") as trace,np.load(directory/"primary.npz") as profile,np.load(directory/"full_states.npz") as full,np.load(directory/"final_state.npz") as final:
            np.testing.assert_array_equal(trace["time"],v.sample_times(manifest["horizon"],manifest["output_interval"]))
            np.testing.assert_array_equal(profile["time"],v.sample_times(manifest["horizon"],manifest["primary_interval"]))
            np.testing.assert_array_equal(full["time"],v.sample_times(manifest["horizon"],manifest["full_interval"]))
            np.testing.assert_array_equal(profile["x"],grid.x)
            for key in trace.files:assert np.isfinite(trace[key]).all(),key
            assert np.isfinite(profile["primary"]).all() and np.isfinite(full["states"]).all()
            np.testing.assert_array_equal(full["states"][0],initial)
            np.testing.assert_allclose(full["states"][-1],final["state"],rtol=0,atol=1e-13)
            np.testing.assert_allclose(profile["primary"][-1],final["state"][0],rtol=0,atol=1e-13)
            assert float(final["time"])==manifest["horizon"]
            indices=np.searchsorted(trace["time"],profile["time"])
            np.testing.assert_array_equal(trace["time"][indices],profile["time"])
            phase=np.exp(2j*np.pi*grid.x/grid.circumference)*grid.dx
            raw1=profile["primary"]@phase;raw2=profile["primary"]@(np.exp(4j*np.pi*grid.x/grid.circumference)*grid.dx)
            stored1=trace["M1_real"][indices]+1j*trace["M1_imag"][indices]
            stored2=trace["M2_real"][indices]+1j*trace["M2_imag"][indices]
            raw_error=float(max(np.max(np.abs(raw1-stored1)),np.max(np.abs(raw2-stored2))))
            np.testing.assert_allclose(raw1,stored1,rtol=1e-12,atol=2e-12)
            np.testing.assert_allclose(raw2,stored2,rtol=1e-12,atol=2e-12)
            np.testing.assert_allclose(raw1/raw1[0],trace["P_real"][indices]+1j*trace["P_imag"][indices],rtol=1e-12,atol=2e-12)
            np.testing.assert_allclose(profile["primary"].sum(axis=1)*grid.dx,trace["M0"][indices],rtol=1e-12,atol=2e-12)
            assert trace["normalized_moment_denominator_defined"].all()
            assert abs(trace["dominant_region_unwrapped_displacement"][0])<1e-12
            pools=v.ring.pool_amounts(model,grid,initial)
            pool_error=float(max(np.max(np.abs(v.ring.pool_amounts(model,grid,state)/pools-1)) for state in full["states"]))
            assert pool_error<manifest["diagnostics"]["pool_relative_tolerance"]
            minimum=float(min(profile["primary"].min(),full["states"].min()))
            assert minimum>=-manifest["diagnostics"]["allowed_negative_roundoff"]
        forcing=v.NoiseInput(grid,manifest["horizon"],job["seed"],job["mean"],kind=job["kind"],**manifest["noise"])
        saved_noise=v.read(directory/"noise_audit.json")
        if job["kind"]=="noise":
            assert v.array_sha(forcing.knots)==saved_noise["gaussian_knots_sha256"]
            assert forcing.coupling_knots_sha256==saved_noise["coupling_knots_sha256"]
        cases.append(dict(slot=job["slot"],model=job["model"],circumference=job["circumference"],kind=job["kind"],
            mean_fraction=job["mean_fraction"],passed=True,raw_moment_max_absolute_error=raw_error,
            sampled_pool_relative_error=pool_error,accepted_step_pool_relative_error=status["max_pool_relative_error"],
            sampled_minimum_concentration=minimum,accepted_step_minimum_concentration=status["minimum_concentration"],
            runtime_seconds=status["elapsed_seconds"],noise_regeneration_exact=True))
    result=dict(generated_at=v.now(),passed=True,run=str(run),manifest_sha256=v.sha(run/"manifest.json"),
        verifier_source_sha256=v.sha(Path(__file__)),frozen_module=str(v.__file__),source_hashes=len(manifest["source_sha256"]),
        cases=cases,allocated_storage_bytes=allocated(run),
        scope="All completed job/output hashes, prescribed time grids, full final states, raw activity moments and P, sampled native budgets, exact Gaussian noise regeneration")
    v.write(run/"results/execution_verification.json",result)
    return result


def compare_pair(reference,alternative,basal,v):
    with np.load(reference/"trace.npz") as a,np.load(alternative/"trace.npz") as b,np.load(reference/"primary.npz") as pa,np.load(alternative/"primary.npz") as pb:
        np.testing.assert_array_equal(a["time"],b["time"]);np.testing.assert_array_equal(pa["time"],pb["time"])
        p=a["P_real"]+1j*a["P_imag"];q=b["P_real"]+1j*b["P_imag"]
        raw=pa["primary"];other=pb["primary"];contrast=float(np.ptp(basal[0]));assert contrast>0
        circumference=float(pa["circumference"]);cells=raw.shape[1]
        shifts=[];aligned=[]
        for first,second in zip(raw,other):
            correlation=np.fft.irfft(np.conj(np.fft.rfft(second))*np.fft.rfft(first),n=cells)
            shift=int(np.argmax(correlation));shift=shift if shift<=cells//2 else shift-cells
            shifts.append(shift*circumference/cells)
            aligned.append(float(np.max(np.abs(first-np.roll(second,shift)))/contrast))
        displacement=(a["dominant_region_position"]-b["dominant_region_position"]+circumference/2)%circumference-circumference/2
        return dict(vector_max_absolute_difference=float(np.max(np.abs(p-q))),vector_rms_difference=float(np.sqrt(np.mean(np.abs(p-q)**2))),
            vector_final_difference=float(abs(p[-1]-q[-1])),maximum_primary_profile_difference_over_basal_contrast=float(np.max(np.abs(raw-other))/contrast),
            rms_primary_profile_difference_over_basal_contrast=float(np.sqrt(np.mean((raw-other)**2))/contrast),
            maximum_integer_aligned_profile_difference_over_basal_contrast=max(aligned),
            maximum_best_integer_alignment_physical_shift=float(np.max(np.abs(shifts))),
            maximum_dominant_region_physical_position_difference=float(np.max(np.abs(displacement))),
            maximum_total_activity_relative_difference=float(np.max(np.abs(a["M0"]-b["M0"]))/a["M0"][0]),
            maximum_q1_difference=float(np.max(np.abs(a["q1"]-b["q1"]))),maximum_q2_difference=float(np.max(np.abs(a["q2"]-b["q2"]))),
            fraction_current_region_count_disagreement=float(np.mean(a["current_half_contrast_regions"]!=b["current_half_contrast_regions"])),
            fraction_prominent_peak_count_disagreement=float(np.mean(a["peaks_prominence_ge_005"]!=b["peaks_prominence_ge_005"])))


def compare(pilot,audit,v):
    m=v.read(pilot/"manifest.json");ode=v.read(audit/"ode/manifest.json");temporal=v.read(audit/"temporal/manifest.json")
    rows=[]
    for tight in ode["jobs"]:
        original=m["jobs"][tight["pilot_slot"]]
        fine=next(j for j in temporal["jobs"] if j["model"]==tight["model"])
        assert tight["seed"]==fine["seed"]==original["seed"]
        assert tight["input_sha256"]==fine["input_sha256"]==original["input_sha256"]
        primary=pilot/"workers"/f'{original["slot"]:04d}';coarse=audit/"ode/workers"/f'{tight["slot"]:04d}';refined=audit/"temporal/workers"/f'{fine["slot"]:04d}'
        for directory in (primary,coarse,refined):assert v.read(directory/"status.json")["complete"]
        noise0=v.read(primary/"noise_audit.json");noise1=v.read(coarse/"noise_audit.json");noise2=v.read(refined/"noise_audit.json")
        assert noise0["gaussian_knots_sha256"]==noise1["gaussian_knots_sha256"]
        assert noise0["coupling_knots_sha256"]==noise1["coupling_knots_sha256"]==noise2["coupling_knots_sha256"]
        with np.load(pilot/original["input"]) as archive:basal=archive["state"].copy()
        rows.append(dict(model=tight["model"],circumference=tight["circumference"],mean_fraction=tight["mean_fraction"],seed=tight["seed"],
            fixed_forcing_ode_sensitivity=compare_pair(primary,coarse,basal,v),
            coupled_interpolation_sensitivity=compare_pair(coarse,refined,basal,v)))
    result=dict(generated_at=v.now(),complete=True,passed_execution=True,pilot_manifest_sha256=v.sha(pilot/"manifest.json"),
        ode_manifest_sha256=v.sha(audit/"ode/manifest.json"),temporal_manifest_sha256=v.sha(audit/"temporal/manifest.json"),
        verifier_source_sha256=v.sha(Path(__file__)),cases=rows,
        interpretation="Fixed-forcing ODE sensitivity and changed stochastic interpolation are reported separately; no automatic scientific-convergence claim or favorable-result gate.")
    v.write(audit/"NUMERICAL_COMPARISON.json",result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("command",choices=("verify","compare"));parser.add_argument("--run",type=Path,default=DEFAULT);parser.add_argument("--audit",type=Path,default=AUDIT)
    args=parser.parse_args();run=args.run.resolve();v=load_frozen(run)
    result=verify(run,v) if args.command=="verify" else compare(run,args.audit.resolve(),v)
    print(json.dumps(v.json_value(result),indent=2,allow_nan=False))


if __name__=="__main__":main()
