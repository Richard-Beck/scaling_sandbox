"""Independent archive and acceptance audit for Spring ring validation.

Reads immutable worker arrays and the run's frozen dynamical source. Does not
rely on worker success labels to establish scientific or numerical acceptance.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import subprocess
import sys

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def times(end):
    return np.unique(np.r_[np.arange(0., min(10., end), .1),
        np.arange(10., min(100., end), 1.), np.arange(100., end, 5.), end])


def vector(y, n):
    return (y[:, :n] @ np.exp(2j*np.pi*np.arange(n)/n)) / np.sum(y[:, :n], axis=1)


def shape(state, n):
    u = state[:n]
    contrast = float(np.ptp(u))
    mask = u >= u.min() + contrast/2.
    components = 0 if contrast < 1e-12 else (1 if mask.all() else int(np.sum(mask & ~np.roll(mask, 1))))
    peaks, _ = find_peaks(np.tile(u, 3), prominence=.05)
    peaks = peaks[(peaks >= n) & (peaks < 2*n)] - n
    vv = vector(state[None], n)[0]
    return dict(contrast=contrast, half_contrast_components=components,
        resolved_peaks=len(peaks), peak_cells=peaks.tolist(), vector_magnitude=float(abs(vv)),
        vector_direction_fraction=float(np.angle(vv)/(2*np.pi) % 1))


def check(run):
    manifest, launch = read(run / "manifest.json"), read(run / "launch.json")
    assert sha(run / "manifest.json") == launch["manifest_sha256"]
    assert sha(run / "inputs/basal.npz") == manifest["input_sha256"] == manifest["original_preparation_sha256"]
    for relative, digest in manifest["source_sha256"].items():
        assert sha(run / "frozen_sources" / relative) == digest
    sys.path.insert(0, str(run / "frozen_sources/code/mechanisms"))
    from spring_ring_model_20261019 import Closure, physical_raised_cosine
    with np.load(run / "inputs/basal.npz") as z:
        retained = z["state"].copy()
    old = ROOT / "data/runs/spring_ring_distance_pilot_20261017"
    assert sha(old / "manifest.json") == manifest["original_distance_manifest_sha256"]
    old_status = read(old / "workers/01/status.json")
    original = {}
    for phase in old_status["phases"]:
        assert sha(old / phase["path"]) == phase["sha256"]
        with np.load(old / phase["path"]) as z:
            original[phase["phase"]] = (z["time"].copy(), z["states"].copy())

    checks, arrays, numerical = [], {}, []
    for job in manifest["jobs"]:
        folder = run / "workers" / f'{job["slot"]:02d}'
        record = read(folder / "status.json")
        assert record["complete"] and record["status"] == "complete"
        assert record["job"] == job
        n = job["cells"]
        closure = Closure(2., n)
        previous = closure.initial() if job["fresh_preparation"] else retained.copy()
        if job["shift_cells"]:
            shift = job["shift_cells"]
            previous = np.r_[np.roll(previous[:n], shift), np.roll(previous[n:-1], shift), previous[-1]]
        if job["perturb"]:
            delta = np.random.default_rng(20261019).uniform(-.001, .001, n)
            previous[:n] += delta
            previous[n:-1] -= delta
            np.testing.assert_allclose(previous[:n]+previous[n:-1],
                retained[:n]+retained[n:-1], rtol=0., atol=8*np.finfo(float).eps)
        expected = []
        if job["fresh_preparation"]:
            expected.extend([("preparation", 10., 4.), ("basal", 13000., 0.)])
        expected.extend([("challenge", job["cue_duration"], 4.), ("recovery", 13000., 0.)]
                        if job["cue_duration"] else [("basal_perturbed", 13000., 0.)])
        assert [p["kind"] for p in record["phases"]] == [p[0] for p in expected]
        result = dict(name=job["name"], phase_checks=[], input_chain_exact=True)
        arrays[job["name"]] = {}
        center = float(closure.material_positions[job["shift_cells"] % n])
        for phase, (kind, end, peak) in zip(record["phases"], expected):
            archive = run / phase["path"]
            assert sha(archive) == phase["sha256"]
            with np.load(archive) as z:
                t, y = z["time"].copy(), z["states"].copy()
                np.testing.assert_array_equal(z["material_positions"], closure.material_positions)
                assert float(z["cue_peak"]) == peak
                assert float(z["total_physical_support"]) == .2
                cue_center = float(z["center_fraction"])
            np.testing.assert_array_equal(t, times(end))
            np.testing.assert_array_equal(y[0], previous)
            assert np.isfinite(y).all() and np.all(np.diff(t) > 0)
            minimum = float(y[:, :-1].min())
            amount = y[:, -1] * np.mean(y[:, :n] + y[:, n:-1], axis=1)
            amount_error = float(np.max(np.abs(amount/4. - 1)))
            assert minimum >= -manifest["negative_tolerance"]
            assert amount_error <= manifest["amount_relative_tolerance"]
            assert float(y[:, -1].min()) >= 2.-1e-10 and float(y[:, -1].max()) <= 2.1+1e-10
            scales = np.r_[np.full(2*n, 2.), 2.]
            tail = t >= end - 1000.
            tail_range = float(np.max(np.ptp(y[tail], axis=0)/scales))
            cue = physical_raised_cosine(closure, y[-1], cue_center, peak, .2)[0]
            rhs = float(np.max(np.abs(closure.rhs(y[-1], cue))/scales)*1000.)
            phase_check = dict(kind=kind, archive_sha256=phase["sha256"],
                complete_expected_time_grid=True, exact_initial_state=True,
                minimum_concentration=minimum, max_amount_relative_error=amount_error,
                final_1000_scaled_full_field_and_C_range=tail_range,
                final_rhs_scaled_to_1000=rhs, final_shape=shape(y[-1], n))
            if peak == 0:
                assert tail_range <= manifest["tail_full_field_scaled_tolerance"]
                assert rhs <= manifest["stationary_rhs_scaled_tolerance"]
                assert phase_check["final_shape"]["contrast"] > .05
                assert phase_check["final_shape"]["half_contrast_components"] == 1
                assert phase_check["final_shape"]["resolved_peaks"] == 1
                phase_check["full_field_and_C_stationary_single_site"] = True
            if kind == "basal":
                windows = []
                for start in (1000., 2000.):
                    use = (t >= start) & (t <= start + 1000.)
                    window = float(np.max(np.ptp(y[use], axis=0)/scales))
                    drift = float(np.max(np.abs(y[use][-1]-y[use][-2])/scales)*1000./5.)
                    assert max(window, drift) <= 1e-4
                    windows.append(dict(start=start, stop=start+1000., range=window, projected_drift=drift))
                hold = t >= 3000.
                excursion = float(np.max(np.abs(y[hold] - y[np.flatnonzero(t == 3000.)[0]])/scales))
                assert excursion <= .001
                phase_check["fresh_basal_qualification"] = dict(two_windows=windows,
                    holdout_start=3000., holdout_duration=10000., max_full_field_and_C_excursion=excursion,
                    passed=True)
            result["phase_checks"].append(phase_check)
            previous = y[-1].copy()
            arrays[job["name"]][kind] = (t, y)
        stderr = run / "logs" / (launch["array_job_id"] + f'_{job["slot"]}.err')
        result["stderr"] = stderr.read_text()
        result["warning_free"] = result["stderr"] == ""
        assert result["warning_free"]
        if job["name"] == "finite_0p1":
            assert record["reusable_protocol_smoke"]["passed"]
            result["reusable_protocol_smoke"] = record["reusable_protocol_smoke"]
        checks.append(result)

    def diff(label, left, right, n=200, grid=False, shift=0):
        ta, a = left
        tb, b = right
        np.testing.assert_allclose(ta, tb, rtol=0, atol=1e-12)
        if grid:
            position = np.arange(401)/400.
            target = np.arange(200)/200.
            b = np.c_[np.asarray([np.interp(target, position, np.r_[row[:400], row[0]]) for row in b]),
                      np.asarray([np.interp(target, position, np.r_[row[400:-1], row[400]]) for row in b]), b[:, -1]]
        if shift:
            b = np.c_[np.roll(b[:, :n], -shift, axis=1), np.roll(b[:, n:-1], -shift, axis=1), b[:, -1]]
        result = dict(label=label, max_chemical_field_difference=float(np.max(np.abs(a[:, :-1]-b[:, :-1]))),
            max_circumference_difference=float(np.max(np.abs(a[:, -1]-b[:, -1]))),
            final_chemical_field_difference=float(np.max(np.abs(a[-1, :-1]-b[-1, :-1]))),
            max_vector_magnitude_difference=float(np.max(np.abs(np.abs(vector(a,n))-np.abs(vector(b,n))))))
        rule = manifest["spatial_comparison_rule"]
        result["within_predeclared_bounded_sensitivity_limits"] = bool(
            result["max_chemical_field_difference"] <= rule["max_field_difference"] and
            result["max_circumference_difference"] <= rule["max_circumference_difference"] and
            result["max_vector_magnitude_difference"] <= rule["max_vector_magnitude_difference"])
        numerical.append(result)

    for phase in ("preparation", "basal", "challenge", "recovery"):
        diff("grid_200_vs_400_" + phase, arrays["fresh_200"][phase], arrays["fresh_400"][phase], grid=True)
        diff("fresh_tolerance_" + phase, arrays["fresh_200"][phase], arrays["fresh_200_tight"][phase])
    for name in ("retained_radau", "retained_tight", "rotation_73"):
        for phase, oldphase in (("challenge", "maintained_cue"), ("recovery", "recovery")):
            diff(name + "_" + phase, original[oldphase], arrays[name][phase], shift=73 if name == "rotation_73" else 0)
    grid_pass = all(v["within_predeclared_bounded_sensitivity_limits"] for v in numerical if v["label"].startswith("grid"))
    sensitivity_pass = all(v["within_predeclared_bounded_sensitivity_limits"] for v in numerical)
    scheduler = subprocess.check_output(["sacct", "-j", launch["array_job_id"],
        "--format=JobID,State,Elapsed,ExitCode", "-n", "-P"], text=True).strip()
    outcome = dict(generated_at=datetime.now(timezone.utc).isoformat(),
        producer="finished", passed=sensitivity_pass,
        status="passed_bounded_validation" if sensitivity_pass else "additional_refinement_required",
        original_scope="C0=2 single-site preparation, imposed cue response and bounded numerical sensitivity",
        manifest_sha256=sha(run / "manifest.json"), auditor_sha256=sha(Path(__file__)),
        jobs=checks, numerical_comparisons=numerical, grid_200_400_pass=grid_pass,
        need_800_cell_check=not grid_pass, scheduler=scheduler,
        interpretation=["Two grids quantify bounded sensitivity and do not prove asymptotic convergence.",
            "One small conserved perturbation is one local recovery observation, not a basin or stability proof.",
            "Warning-free fresh trajectories do not establish the cause of the historical BDF warning.",
            "The force-average closure remains a new global-extension analogue, not the original 2D mechanics."])
    target = run / "results/audit.json"
    target.write_text(json.dumps(outcome, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k:v for k,v in outcome.items() if k not in ("jobs", "scheduler")}, indent=2))
    return outcome


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, default=ROOT / "data/runs/spring_ring_validation_20261019")
    args = parser.parse_args()
    check(args.run.resolve())
