"""Render the revised Figure 1 from immutable audited native comparison inputs.

No integration or retained-input mutation. Each chemical pool is assembled before
its mass-normalized reflection Wasserstein distance is computed; the plotted
burden is the sum of those pool distances, not a sum of distances for activation
states. Supporting exports retain the amount-weighted counterpart separately.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, MaxNLocator, NullFormatter
import numpy as np
from scipy.stats import wasserstein_distance

ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "data/runs/debelly_corrected_three_test_20261008"
OUT = ROOT / "figures/length_mechanisms_20261008/revision_20261010/figure1"
STEM = "length_responses"
MODELS = ["wave_pinning", "goryachev", "debelly", "holmes_model3", "spring"]
LABELS = dict(zip(MODELS, ["Wave-pinning", "Goryachev", "deBelly", "Holmes 3", "Spring"]))
CUES = [
    ("fixed_physical", "Fixed width", "#0072B2", "-"),
    ("fixed_relative", "20% width", "#D55E00", "-"),
    ("fixed_slope", "Constant slope", "#8B3F82", "--"),
]
POOL_DEFINITIONS = {
    "wave_pinning": {"total_active_inactive": ["active", "inactive"]},
    "goryachev": {"total_active_inactive": ["active", "inactive"]},
    "debelly": {"Rac": ["rac"], "Rho": ["rho"], "total_MCA": ["mca_bound", "mca_unbound"]},
    "holmes_model3": {"total_Rac": ["rac_active", "rac_inactive"], "total_Rho": ["rho_active", "rho_inactive"]},
    "spring": {"total_active_inactive": ["active", "inactive"]},
}
SOURCES: dict[str, str] = {}
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.labelsize": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": .7, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 10, "pdf.fonttype": 42, "svg.fonttype": "none",
})


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def register(path: Path, expected: str | None = None) -> None:
    digest = sha(path)
    if expected is not None:
        assert digest == expected, f"Input hash mismatch: {path}"
    SOURCES[str(path.relative_to(ROOT))] = digest


def read_json(path: Path) -> dict:
    register(path)
    return json.loads(path.read_text())


def write_csv(name: str, records: list[dict]) -> None:
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def boolean(value: str | bool | None) -> bool | None:
    if value in [None, "", "None", "null"]:
        return None
    if isinstance(value, bool):
        return value
    assert value.lower() in ["true", "false", "1", "0"]
    return value.lower() in ["true", "1"]


def duration_inputs(manifest: dict) -> tuple[list[dict], dict]:
    """Independently reselect minima from all stored native trial summaries."""
    status = read_json(CAMPAIGN / "results/status.json")
    original_audit = read_json(CAMPAIGN / "results/independent_audit.json")
    assert status["complete"] and status["execution_errors"] == 0
    assert original_audit["passed"] and original_audit["complete"]
    assert not original_audit["errors"]
    path = CAMPAIGN / "results/endpoints.csv"
    register(path)
    with path.open() as handle:
        rows = [row for row in csv.DictReader(handle) if row["definition"] == "full_polarity"]
    assert len(rows) == len(manifest["cases"]) == 330
    path = CAMPAIGN / "results/records.json.gz"
    register(path)
    with gzip.open(path, "rt") as handle:
        native_records = json.load(handle)
    cases = {case["case_index"]: case for case in native_records["cases"]}
    independently_audited = {
        record["case_index"]: record for record in original_audit["selected_readouts"]
        if record["definition"] == "full_polarity"
    }
    output = []
    counts = {"success": 0, "censored": 0, "ineligible": 0}
    for row in rows:
        index = int(row["case_index"])
        case = cases[index]
        assert case["model"] == row["model"]
        assert float(row["length"]) == case["length"]
        assert case["condition"] == row["condition"]
        audited = independently_audited[index]
        assert audited["status"] == row["status"]
        assert row["status"] in counts
        endpoint = float(row["endpoint"]) if row["endpoint"] else None
        field_flag = boolean(row["selected_field_tail_pass"])
        if case["eligible"]:
            assert case["preparation_valid"]
            passed = [test for test in case["tests"] if test["duration"] > 0 and test["full_polarity_success"]]
            selected = min(passed, key=lambda test: test["duration"]) if passed else None
            if selected:
                assert row["status"] == "success"
                assert selected["duration"] == endpoint == audited["endpoint"]
                assert selected["field_tail_pass"] == field_flag
                assert selected["orientation_success"]
                assert selected["minimum_contrast_last_500"] >= .05
                assert selected["orientation_range_last_500"] < .005
                assert selected["maximum_relative_polarity_error_last_500"] <= .01
                assert selected["final_orientation"] > .05
            else:
                assert row["status"] == "censored" and endpoint is None
            assert not any(test["full_polarity_success"] for test in case["tests"] if test["duration"] == 0)
        else:
            assert row["status"] == "ineligible" and endpoint is None
        counts[row["status"]] += 1
        output.append({
            "case_index": index, "model": row["model"], "length": case["length"],
            "condition": row["condition"], "gate": "full_polarity",
            "status": row["status"], "cue_duration": endpoint,
            "tested_limit": float(row["maximum_tested_duration"]) if row["maximum_tested_duration"] else None,
            "field_tail_pass": field_flag, "checkpoint": case["checkpoint"],
            "checkpoint_sha256": case["checkpoint_sha256"],
        })
    return output, {"passed": True, "readouts": len(output), "status_counts": counts,
                    "selection": "Shortest positive native tested duration passing canonical 1%-strength gate; native field flag never filters."}


def chemical_measurements(manifest: dict) -> tuple[list[dict], list[dict], list[dict], dict]:
    pool_rows, total_rows, profiles = [], [], []
    largest_error = 0.0
    checkpoint_count = 0
    for model in MODELS:
        cases = sorted([case for case in manifest["cases"]
                        if case["model"] == model and case["condition"] == "fixed_relative" and case["eligible"]],
                       key=lambda case: case["length"])
        for case in cases:
            siblings = [sibling for sibling in manifest["cases"]
                        if sibling["model"] == model and sibling["length"] == case["length"]]
            assert len(siblings) == 3
            assert all(sibling["checkpoint_sha256"] == case["checkpoint_sha256"] for sibling in siblings)
            path = ROOT / case["checkpoint"]
            register(path, case["checkpoint_sha256"])
            with np.load(path) as archive:
                fields = {key[6:]: archive[key][-1].copy() for key in archive.files if key.startswith("field_")}
                native_x = archive["x"].copy()
                active = archive["active"][-1].copy()
                if model == "spring":
                    state_active, state_inactive = archive["state"][:-2:2], archive["state"][1:-2:2]
                    np.testing.assert_array_equal(fields["active"], state_active)
                    np.testing.assert_array_equal(fields["inactive"], state_inactive)
                elif model == "debelly":
                    for name in ["rac", "rho", "mca_bound", "mca_unbound"]:
                        np.testing.assert_array_equal(fields[name], archive["state_" + name])
                else:
                    names = ["rac_active", "rac_inactive", "rho_active", "rho_inactive"] if model == "holmes_model3" else ["active", "inactive"]
                    for field_index, name in enumerate(names):
                        np.testing.assert_array_equal(fields[name], archive["state"][field_index])
            checkpoint_count += 1
            L0 = float(case["length"])
            L = float(fields["length"]) if model == "spring" else L0
            n = len(native_x)
            if model == "spring":
                x = np.linspace(0, L, n)
                weights = np.full(n, L / (n - 1))
                weights[[0, -1]] *= .5
                # Stored contract x is the prescribed resting grid; chemistry
                # lives at material nodes mapped onto the final physical domain.
                np.testing.assert_allclose(native_x / L0, np.linspace(0, 1, n), atol=1e-12)
            else:
                x = (np.arange(n) + .5) * L / n
                weights = np.full(n, L / n)
                np.testing.assert_allclose(native_x, x, atol=1e-12)
            orientation = float(np.dot(active * weights, 2*x/L - 1) / np.dot(active, weights))
            np.testing.assert_allclose(orientation, case["initial_orientation"], atol=1e-11, rtol=1e-9)
            sum_distance = sum_cost = 0.0
            for pool, names in POOL_DEFINITIONS[model].items():
                concentration = sum((fields[name] for name in names), start=np.zeros(n))
                assert concentration.shape == x.shape
                assert np.isfinite(concentration).all() and concentration.min() >= 0
                mass = concentration * weights
                target_mass = mass[::-1]
                amount = float(mass.sum())
                assert amount > 0
                # Exact W1 for the native weighted point-mass measures. Native
                # FV cell averages and Spring node quadrature are unchanged.
                distance = float(np.dot(np.abs(np.cumsum((mass-target_mass)/amount))[:-1], np.diff(x)))
                independently_computed = float(wasserstein_distance(x, x, u_weights=mass, v_weights=target_mass))
                error = abs(distance - independently_computed)
                largest_error = max(largest_error, error)
                np.testing.assert_allclose(distance, independently_computed, rtol=3e-12, atol=3e-12)
                # Reflection preserves amount; the centroid displacement is a
                # lower bound on W1 and independently verifies the coordinates.
                source_centroid = float(np.dot(mass, x)/amount)
                target_centroid = float(np.dot(target_mass, x)/amount)
                assert abs(target_centroid - source_centroid) <= distance + 1e-11
                np.testing.assert_allclose(source_centroid + target_centroid, L, atol=1e-11)
                pool_rows.append({
                    "model": model, "length": L0, "distance_domain_length": L,
                    "chemical_pool": pool, "summed_fields": "+".join(names),
                    "amount": amount, "reflection_W1": distance,
                    "amount_times_W1": amount*distance,
                    "initial_centroid": source_centroid, "reflected_centroid": target_centroid,
                    "independent_scipy_W1": independently_computed,
                    "checkpoint": case["checkpoint"], "checkpoint_sha256": case["checkpoint_sha256"],
                })
                sum_distance += distance
                sum_cost += amount * distance
                for position, value, weight in zip(x, concentration, weights):
                    profiles.append({"model": model, "length": L0, "chemical_pool": pool,
                                     "position": float(position), "total_concentration": float(value),
                                     "quadrature_weight": float(weight)})
            total_rows.append({"model": model, "length": L0, "distance_domain_length": L,
                               "number_chemical_pools": len(POOL_DEFINITIONS[model]),
                               "sum_pool_W1": sum_distance, "sum_amount_times_W1": sum_cost,
                               "checkpoint": case["checkpoint"], "checkpoint_sha256": case["checkpoint_sha256"]})
    assert checkpoint_count == 97
    trends = []
    for model in MODELS:
        rows = [row for row in total_rows if row["model"] == model]
        deltas = np.diff([row["sum_pool_W1"] for row in rows])
        assert np.all(deltas > 0), (model, deltas)
        first, last = rows[0], rows[-1]
        trends.append({"model": model, "n_lengths": len(rows), "short_length": first["length"], "long_length": last["length"],
                       "short_sum_pool_W1": first["sum_pool_W1"], "long_sum_pool_W1": last["sum_pool_W1"],
                       "long_short_ratio": last["sum_pool_W1"] / first["sum_pool_W1"],
                       "every_adjacent_interval_increases": True})
    return pool_rows, total_rows, profiles, {
        "passed": True, "checkpoints": checkpoint_count, "pools": len(pool_rows),
        "largest_absolute_scipy_W1_discrepancy": largest_error,
        "checkpoint_hashes_verified": True, "shared_three_cue_starts_verified": True,
        "fields_match_full_checkpoint_state": True, "native_initial_orientation_recomputed": True,
        "reflection_amount_and_centroid_verified": True, "pool_definitions": POOL_DEFINITIONS,
        "coordinate_convention": "Native FV nominal cell centres; Spring current physical node positions with trapezoidal weights.",
        "excluded_fields": "Mechanics, endpoint positions, derived activation fields; deBelly has no explicit inactive Rac/Rho pools. LEGI and Otsuji excluded from the native reversal comparison.",
        "trends": trends,
    }


def render(duration_rows: list[dict], distance_rows: list[dict], manifest: dict) -> None:
    fig, axes = plt.subplots(2, 5, figsize=(15.3, 6.8), sharex="col", gridspec_kw={"height_ratios": [1.23, 1]})
    fig.subplots_adjust(left=.066, right=.987, top=.835, bottom=.12, hspace=.28, wspace=.34)
    for col, model in enumerate(MODELS):
        top, bottom = axes[:, col]
        top.set_title(LABELS[model], pad=12)
        valid, limits = [], []
        for cue, _, color, linestyle in CUES:
            rows = sorted([row for row in duration_rows if row["model"] == model and row["condition"] == cue], key=lambda row: row["length"])
            x = np.array([row["length"] for row in rows])
            y = np.array([row["cue_duration"] if row["status"] == "success" else np.nan for row in rows])
            top.plot(x, y, color=color, ls=linestyle, lw=1.9, zorder=3)
            for row in rows:
                if row["status"] == "success":
                    valid.append(row["cue_duration"])
                    top.plot(row["length"], row["cue_duration"], "o", color=color, ms=3.5,
                             mfc="white" if row["field_tail_pass"] is False else color,
                             mew=.8, zorder=4)
                elif row["status"] == "censored":
                    limits.append(row["tested_limit"])
                    top.plot(row["length"], row["tested_limit"], "^", color=color, ms=5,
                             mfc="white", mew=1.1, zorder=4)
        top.set_yscale("linear" if model == "debelly" else "log")
        top.set_ylim((0, max(valid + limits)*1.12) if model == "debelly"
                     else (min(valid)/1.35, max(valid + limits)*1.45))
        lo, hi = top.get_ylim()
        if model != "debelly" and hi/lo >= 20:
            top.yaxis.set_major_locator(LogLocator(base=10, numticks=4))
            top.yaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value:g}"))
        else:
            top.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=2))
            top.yaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value:g}"))
        top.yaxis.set_minor_formatter(NullFormatter())
        top.tick_params(axis="y", which="minor", length=2, color="#bbbbbb")
        # One shared eligibility strip per model, not three overprinted marks.
        excluded = sorted({case["length"] for case in manifest["cases"] if case["model"] == model and not case["eligible"]})
        for length in excluded:
            top.plot(length, -.085, "x", color="#9a9a9a", ms=4, mew=.9,
                     transform=top.get_xaxis_transform(), clip_on=False)
        rows = [row for row in distance_rows if row["model"] == model]
        x = [row["length"] for row in rows]
        y = [row["sum_pool_W1"] for row in rows]
        bottom.fill_between(x, y, color="#d5e8e4", alpha=.85, linewidth=0)
        bottom.plot(x, y, color="#007D73", lw=2.2, marker="o", ms=3.5, mew=0, zorder=3)
        bottom.set_ylim(0, max(y)*1.12)
        bottom.yaxis.set_major_locator(MaxNLocator(nbins=4))
        bottom.yaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value:g}"))
        all_lengths = [case["length"] for case in manifest["cases"] if case["model"] == model]
        span = max(all_lengths) - min(all_lengths)
        bottom.set_xlim(min(all_lengths)-.035*span, max(all_lengths)+.035*span)
        bottom.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
        bottom.set_xlabel("Resting length L0" if model == "spring" else "Cell length L")
        for ax in [top, bottom]:
            ax.grid(axis="y", color="#e1e5e8", lw=.55, zorder=0)
            ax.tick_params(direction="out", length=3, width=.7, color="#65727c")
            ax.spines["left"].set_color("#879198")
            ax.spines["bottom"].set_color("#879198")
        if col == 0:
            top.set_ylabel("Reversal cue duration")
            bottom.set_ylabel("Chemical reorganization\ndistance")
    # Exactly two letters refer to the complete faceted rows.
    fig.text(.018, .843, "A", fontsize=16, fontweight="bold")
    fig.text(.018, .405, "B", fontsize=16, fontweight="bold")
    cue_handles = [Line2D([], [], color=color, lw=2, ls=linestyle, label=label) for _, label, color, linestyle in CUES]
    fig.legend(handles=cue_handles, loc="upper center", bbox_to_anchor=(.52, .982),
               ncol=3, frameon=False, handlelength=2.8, columnspacing=2.8)
    status_handles = [
        Line2D([], [], color="#626e77", marker="o", mfc="white", ls="none", ms=5, label="Field flag fails"),
        Line2D([], [], color="#626e77", marker="^", mfc="white", ls="none", ms=5, label="Tested limit"),
        Line2D([], [], color="#9a9a9a", marker="x", ls="none", ms=5, label="Ineligible start"),
    ]
    fig.legend(handles=status_handles, loc="lower center", bbox_to_anchor=(.52, .018),
               ncol=3, frameon=False, fontsize=9, columnspacing=2.6)
    for extension in ["png", "pdf", "svg"]:
        fig.savefig(OUT / f"{STEM}.{extension}", dpi=240, facecolor="white")
    plt.close(fig)


CAPTION = (
    "*Figure 1* Larger cells have a greater chemical reorganization burden even when their reversal cues become shorter. "
    "A) Shortest tested finite cue duration that restores reversed scalar polarity within 1% of its initial magnitude, "
    "for all three cue geometries in each model. Fixed physical support is 3 native length units (Spring, 0.2); "
    "the relative patch covers 20% of the cue domain; the gradient has constant physical slope. The two patch peaks "
    "remain fixed with length, whereas the gradient peak grows with length. Duration axes use native time; deBelly "
    "is linear and the other four are logarithmic. Filled/open circles retain passing/failing native-field tail flags without changing "
    "the scalar selection; open triangles mark censoring at the tested limit; gray crosses mark preparation-ineligible "
    "lengths. Each pulse is followed by 10000 native recovery units; that recovery time is excluded from the cue duration. "
    "B) Chemical reorganization distance from the exact initial basal state to its spatially reflected target. "
    "Active and inactive concentrations are added pointwise before computing the mass-normalized Wasserstein distance "
    "for each chemical pool; the plotted distance is the sum across pools. Wave-pinning, Goryachev and Spring each "
    "supply one active+inactive pool. Holmes supplies total Rac and total Rho. deBelly supplies Rac, Rho and total "
    "bound+unbound MCA; it has no explicit inactive Rac or Rho. Mechanical fields and derived reaction rates are excluded. "
    "Each pool retains its basal background. Distances use native nominal coordinates for the fixed-domain models "
    "and the steady physical domain for Spring; Spring's x axis shows resting length L0. Units differ between models. "
    "The sum increases at every adjacent eligible length in all five models. It measures spatial reorganization of "
    "chemical totals, not elapsed time or an amount-weighted transport cost."
)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = read_json(CAMPAIGN / "manifest.json")
    assert manifest["models"] == MODELS and manifest["center_boundary_motion"] is False
    assert {"legi", "otsuji"}.issubset(set(manifest["excluded_models"]))
    durations, duration_audit = duration_inputs(manifest)
    pools, totals, profiles, chemical_audit = chemical_measurements(manifest)
    write_csv("cue_duration_1percent.csv", durations)
    write_csv("chemical_pool_distances.csv", pools)
    write_csv("chemical_reorganization_distance.csv", totals)
    write_csv("chemical_total_profiles.csv", profiles)
    write_csv("chemical_reorganization_length_summary.csv", chemical_audit["trends"])
    render(durations, totals, manifest)
    (OUT / f"{STEM}_caption.txt").write_text(CAPTION + "\n")
    stamp = datetime.now(timezone.utc).isoformat()
    audit = {"generated_at": stamp, "producer": "finished", "passed": True,
             "duration_audit": duration_audit, "chemical_audit": chemical_audit}
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    artifacts = {path.name: sha(path) for path in sorted(OUT.iterdir())
                 if path.is_file() and path.name != f"{STEM}.manifest.json"}
    record = {
        "generated_at": stamp, "producer": "finished", "passed": True,
        "renderer": str(Path(__file__).resolve().relative_to(ROOT)), "renderer_sha256": sha(Path(__file__)),
        "artifacts": artifacts, "sources": dict(SOURCES), "legend": CAPTION,
        "description": "Two-row five-model Figure 1: three cue curves under only the canonical 1%-strength gate; summed chemical-pool reflection W1.",
        "methods": "Native checkpoint hash and chemical state identity audit; exact shared three-cue starts; independent all-trial-summary threshold reselection; native initial polarity; CDF W1 independently checked with scipy.stats.wasserstein_distance; native FV or Spring trapezoidal weights.",
        "metric": "sum over chemical pools of W1(normalized total pool, reflected normalized total pool)",
        "pool_definitions": POOL_DEFINITIONS,
        "reproduction": "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python " + str(Path(__file__).resolve()),
    }
    (OUT / f"{STEM}.manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({"output": str(OUT), "passed": True, "duration_audit": duration_audit,
                      "chemical_trends": chemical_audit["trends"]}, indent=2))


if __name__ == "__main__":
    main()
