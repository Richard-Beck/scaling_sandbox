"""Read-only first-angular-moment audit of a declared retained ring subset.

No solver is called. Scientific inputs and prior figures are never modified.
"""
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

MANUSCRIPT = Path(__file__).resolve().parents[2]
REPO = MANUSCRIPT.parents[1]
OUT = MANUSCRIPT / "figures/length_mechanisms_20261008/dev/global_polarity_vector_20261017"
MATRIX = MANUSCRIPT / "data/runs/chemical_ring_functional_v2_20261017/results/matrix.csv"
CASES = [
    ("wave_pinning", 40, .4, "common_C40_d16"),
    ("goryachev", 40, .4, "common_C40_d16"),
    ("holmes_model3", 40, .4, "common_C40_d16_and_Holmes_location_control"),
    ("holmes_model3", 40, .3, "Holmes_location_control"),
    ("holmes_model3", 30, .4, "Holmes_cue_single_example"),
    ("goryachev", 40, .5, "Goryachev_antipodal_size_control"),
    ("goryachev", 80, .5, "Goryachev_antipodal_size_control"),
    ("holmes_model3", 40, .5, "Holmes_antipodal_cancellation_diagnostic"),
]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def moment(primary, x, circumference, order=1):
    dx = circumference / len(x)
    theta = 2 * np.pi * (x - x[0]) / circumference
    return dx * (primary @ np.exp(1j * order * theta))


def analytic_checks():
    n = 4096
    circumference = 2 * np.pi
    x = (np.arange(n) + .5) * circumference / n
    theta = 2 * np.pi * (x - x[0]) / circumference
    delta = .7
    baseline = 3 + np.cos(theta)
    shifted = 3 + np.cos(theta - delta)
    two_lobes = 3 + np.cos(2 * theta)
    ambiguous = baseline + np.cos(2 * theta)
    m0 = moment(baseline, x, circumference)
    checks = {
        "uniform_first_moment_abs": float(abs(moment(np.ones(n), x, circumference))),
        "shifted_analytic_moment_error": float(abs(moment(shifted, x, circumference) - np.pi * np.exp(1j * delta))),
        "uniform_addition_moment_error": float(abs(moment(baseline + 8, x, circumference) - m0)),
        "rotated_unit_vector_error": float(abs(moment(shifted, x, circumference) / m0 - np.exp(1j * delta))),
        "two_lobe_first_moment_abs": float(abs(moment(two_lobes, x, circumference))),
        "two_lobe_second_moment_error": float(abs(moment(two_lobes, x, circumference, 2) - np.pi)),
        "different_peak_structure_same_first_moment_error": float(abs(moment(ambiguous, x, circumference) - m0)),
        "scale_invariance_error": float(abs(moment(7 * shifted, x, circumference) / moment(7 * baseline, x, circumference) - np.exp(1j * delta))),
        "nonnegative_test_profile_min": float(min(baseline.min(), shifted.min(), two_lobes.min(), ambiguous.min())),
    }
    assert max(v for k, v in checks.items() if k.endswith(("error", "abs"))) < 1e-12
    return checks


def load(row, key, offset):
    with np.load(row[key]) as z:
        time = z["time"].copy() + offset
        x = z["x"].copy()
        primary = z["states"][:, 0].copy() if "states" in z else z["primary"].copy()
        legacy = (z["harmonic_real"] + 1j * z["harmonic_imag"]) if "harmonic_real" in z else None
    return time, x, primary, legacy


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    for child in ("vectors", "data", "metadata", "viewers"):
        (OUT / "supporting" / child).mkdir(parents=True, exist_ok=True)
    with MATRIX.open() as stream:
        matrix = list(csv.DictReader(stream))
    summaries, traces, sources = [], [], []
    for model, circumference, fraction, role in CASES:
        row, = [r for r in matrix if r["model"] == model and float(r["circumference"]) == circumference and np.isclose(float(r["fraction"]), fraction)]
        assert row["status"] == "complete"
        tc, x, pc, lc = load(row, "cue_path", 0)
        tr, xr, pr, lr = load(row, "recovery_path", 1000)
        np.testing.assert_array_equal(x, xr)
        np.testing.assert_array_equal(pc[-1], pr[0])
        t = np.r_[tc, tr[1:]]
        primary = np.concatenate((pc, pr[1:]), axis=0)
        assert t[0] == 0 and t[-1] == 14000 and np.all(np.diff(t) > 0)
        dx = circumference / len(x)
        total = primary.sum(axis=1) * dx
        m1 = moment(primary, x, circumference)
        m2 = moment(primary, x, circumference, order=2)
        assert abs(m1[0]) / total[0] > 1e-8
        vector = m1 / m1[0]
        magnitude = np.abs(vector)
        angle = np.angle(vector)
        coherence = abs(m1) / total
        legacy_error = None
        if lc is not None and lr is not None:
            legacy = np.r_[lc, lr[1:]]
            legacy_error = float(np.max(abs(abs(legacy) - coherence)))
            assert legacy_error < 1e-12
        summary = dict(model=model, circumference=circumference, fraction=fraction,
                       shorter_route=float(row["short_route"]), declared_role=role,
                       final_geographic_class=row["final_location_association"],
                       frames=len(t), max_recovery_sample_gap=float(np.diff(tr).max()),
                       basal_first_moment_abs=float(abs(m1[0])), basal_total_active=float(total[0]),
                       basal_mass_normalized_magnitude=float(coherence[0]),
                       minimum_relative_magnitude=float(magnitude.min()),
                       time_at_minimum_magnitude=float(t[np.argmin(magnitude)]),
                       final_relative_magnitude=float(magnitude[-1]),
                       final_vector_distance_from_initial=float(abs(vector[-1] - 1)),
                       minimum_primary=float(primary.min()), legacy_magnitude_max_error=legacy_error)
        for gate in (.05, .1, .2):
            eligible = magnitude >= gate
            label = str(gate).replace(".", "p")
            summary[f"sampled_max_angle_deg_gate_{label}"] = float(np.rad2deg(abs(angle[eligible])).max())
            summary[f"final_angle_deg_gate_{label}"] = float(np.rad2deg(angle[-1])) if eligible[-1] else None
            summary[f"masked_frames_gate_{label}"] = int((~eligible).sum())
        for time in (10, 1000, 2000, 14000):
            indices, = np.where(t == time)
            assert len(indices) == 1
            idx = int(indices[0])
            summary[f"magnitude_t{time}"] = float(magnitude[idx])
            summary[f"angle_deg_t{time}_gate_0p1"] = float(np.rad2deg(angle[idx])) if magnitude[idx] >= .1 else None
        summaries.append(summary)
        for i, time in enumerate(t):
            traces.append(dict(model=model, circumference=circumference, fraction=fraction, time=float(time),
                               vector_x=float(vector[i].real), vector_y=float(vector[i].imag),
                               relative_magnitude=float(magnitude[i]),
                               angle_deg_gate_0p1=float(np.rad2deg(angle[i])) if magnitude[i] >= .1 else None,
                               mass_normalized_first_magnitude=float(coherence[i]),
                               mass_normalized_second_magnitude=float(abs(m2[i]) / total[i]),
                               total_active=float(total[i])))
        for key in ("cue_path", "recovery_path"):
            path = Path(row[key])
            sources.append(dict(model=model, circumference=circumference, fraction=fraction, phase=key,
                                path=str(path.relative_to(REPO)), bytes=path.stat().st_size, sha256=sha(path)))
    write_csv(OUT / "supporting/data/global_polarity_summary.csv", summaries)
    write_csv(OUT / "supporting/data/global_polarity_traces.csv", traces)
    write_csv(OUT / "supporting/data/global_polarity_sources.csv", sources)
    audit = dict(generated_at=datetime.now(timezone.utc).isoformat(), producer="completed; retained-data analysis only",
                 cases=len(CASES), sampled_frames=len(traces), analytic_checks=analytic_checks(),
                 no_new_integration=True, matrix_sha256=sha(MATRIX), generator_sha256=sha(Path(__file__)),
                 angle_gate_relative_to_initial_moment=.1, angle_gate_sensitivities=[.05, .2],
                 figure_selection="Not a new matrix-wide analysis; the eight cases were declared before computing moments.")
    audit["outputs"] = {str(p.relative_to(OUT)): sha(p) for p in (OUT / "supporting/data").glob("global_polarity_*.csv")}
    (OUT / "supporting/metadata/global_polarity_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
