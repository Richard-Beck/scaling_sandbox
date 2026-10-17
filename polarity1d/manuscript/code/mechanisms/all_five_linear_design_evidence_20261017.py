"""Algebraic cue-channel evidence for the all-five linear assay design.

This program reads retained checkpoints and evaluates native reactions only.
It does not integrate, perturb, replace or modify any retained model/state.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code/core"))
import native as n

OUT = ROOT / "figures/length_mechanisms_20261008/dev/all_five_linear_function_20261017"
SOURCE = ROOT / "data/runs/initial_polarity_length_corrected_20261008"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    records = list(csv.DictReader((SOURCE / "steady_state_metrics.csv").open()))
    provenance = json.loads((SOURCE / "report_provenance.json").read_text())
    assert sha(SOURCE / "steady_state_metrics.csv") == provenance["summary_csv_sha256"]
    rows, checks = [], []
    for row in records:
        key = row["model"]
        checkpoint = ROOT / row["checkpoint"]
        assert sha(checkpoint) == row["checkpoint_sha256"]
        with np.load(checkpoint) as z:
            fields = {name[6:]: z[name][-1].copy() for name in z.files if name.startswith("field_")}
            state = None if key == "debelly" else z["state"].copy()
            x = z["x"].copy()
        physical = float(row["physical_length"])
        if key == "spring":
            x = np.linspace(0, physical, len(fields["active"]))
            weights = np.full(len(x), physical / (len(x) - 1))
            weights[[0, -1]] *= .5
        else:
            weights = np.full(len(x), physical / len(x))
        if key == "debelly":
            rho = fields["rho"]
            coefficient = .04 / (1 + rho**2)
            # Rac source is linear in alpha; the external cue enters alpha additively.
            parameters = n.ds._parameters()
            alpha0, _ = n.ds._derived_rates(parameters, fields["mca_bound"], fields["tension"], np.zeros(len(x)))
            alpha1, _ = n.ds._derived_rates(parameters, fields["mca_bound"], fields["tension"], np.ones(len(x)))
            actual = parameters.rac_deactivation * (alpha1 - alpha0) / (1 + rho**2)
            channel = "Rac activation alpha; no explicit inactive Rac pool"
        elif key == "spring":
            u, v = fields["active"], fields["inactive"]
            coefficient = v
            delta = 3 + 80 * (physical - float(row["domain_length"]))
            actual = n.sp.auth["f"](u, v, 5., 5, 6, 2, delta) - n.sp.auth["f"](u, v, 4., 5, 6, 2, delta)
            channel = "basal activation b; conserved active/inactive transfer"
        else:
            model = n.model_for(key)
            coefficient = state[1] / (model.rac_total if key == "holmes_model3" else 1.)
            actual = model.reaction(0., x, state, np.ones(len(x)))[0] - model.reaction(0., x, state, np.zeros(len(x)))[0]
            channel = "Rac GEF activation" if key == "holmes_model3" else "activation of explicit inactive pool"
        error = float(np.max(np.abs(actual - coefficient)))
        np.testing.assert_allclose(actual, coefficient, rtol=1e-11, atol=1e-12)
        width = .2 if key == "spring" else 3.
        distance = physical - x
        profile = np.where(distance <= width, .5 * (1 + np.cos(np.pi * np.minimum(distance / width, 1))), 0.)
        denominator = float(profile @ weights)
        assert denominator > 0
        gain = float((profile * coefficient) @ weights / denominator)
        rows.append(dict(model=key, size_native=float(row["domain_length"]), physical_length=physical,
                         end_to_zero_support=width, opposite_cell_gain_native=float(coefficient[-1]),
                         cue_weighted_gain_native=gain, channel=channel,
                         checkpoint=row["checkpoint"], checkpoint_sha256=row["checkpoint_sha256"],
                         native_increment_identity_max_error=error))
    pairs = {"wave_pinning": (19, 39), "goryachev": (19, 39), "debelly": (19, 39),
             "holmes_model3": (20, 40), "spring": (1, 2)}
    for key, pair in pairs.items():
        group = [r for r in rows if r["model"] == key]
        a, b = [next(r for r in group if r["size_native"] == size) for size in pair]
        checks.append(dict(model=key, cases=len(group), size_small=pair[0], size_large=pair[1],
                           point_gain_small=a["opposite_cell_gain_native"], point_gain_large=b["opposite_cell_gain_native"],
                           point_gain_large_over_small=b["opposite_cell_gain_native"] / a["opposite_cell_gain_native"],
                           cue_weighted_gain_small=a["cue_weighted_gain_native"], cue_weighted_gain_large=b["cue_weighted_gain_native"],
                           cue_weighted_gain_large_over_small=b["cue_weighted_gain_native"] / a["cue_weighted_gain_native"]))
    folder = OUT / "supporting/metadata"
    folder.mkdir(parents=True, exist_ok=True)
    data = OUT / "supporting/data"
    data.mkdir(exist_ok=True)
    with (data / "native_gain_design_evidence.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    summary = dict(generated_at=datetime.now(timezone.utc).isoformat(), new_integrations=0,
                   retained_checkpoint_hashes_checked=len(rows), cases=checks,
                   maximum_native_increment_identity_error=max(r["native_increment_identity_max_error"] for r in rows),
                   source_sha256={str(p.relative_to(ROOT)): sha(p) for p in (
                       Path(__file__), SOURCE / "steady_state_metrics.csv", SOURCE / "report_provenance.json",
                       ROOT / "code/core/models/wave_pinning.py", ROOT / "code/core/models/goryachev.py",
                       ROOT / "code/core/models/holmes_model3.py", ROOT / "code/core/models/spring.py",
                       ROOT / "code/core/models/debelly.py")})
    (folder / "native_gain_design_evidence.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "source_sha256"}, indent=2))


if __name__ == "__main__":
    main()
