"""Replay frozen S0 criteria with round-2 rate ratios and noise tolerances.

No model integration, gate retuning, or main-figure selection occurs here.
H3 and the explicitly direct-tau A6 gate remain their original estimands.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import signature_figure4_final_statistics_20261020 as F

csv.field_size_limit(50_000_000)


def read(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def encode(value):
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: encode(v) if isinstance(v, (dict, list, tuple)) else v for k, v in row.items()})


def key(row):
    return row["arm"], row["view"], F.boolean(row["forced_inclusive"]), row["model"], row["criterion"]


def mapped_ratio(row):
    """Rename the estimand only; retain forcing/coverage and frozen truth rules."""
    result = dict(row)
    for suffix in ("", "_lower", "_upper", "_ci_lower", "_ci_upper", "_ci_upper_unbounded"):
        result["tau_ratio" + suffix] = row.get("rho_k" + suffix, "")
    return result


def evaluate(old, ratios, tolerances):
    old_index = {key(row): row for row in old}
    expected = {(a, v, inc, m, c) for a in ("common", "longest") for v in F.VIEWS for inc in (False, True)
                for m, c in [(m, c) for m in F.MODELS for c in ("H1_model", "H2_model")]
                + [("all", c) for c in ("H1_overall", "H2_overall", "H3_overall")]}
    if len(old) != 308 or set(old_index) != expected:
        raise ValueError("Original S0 does not contain exactly the 308 required rows")
    rkeys = [(r["arm"], r["model"], r["view"], float(r["f"])) for r in ratios]
    if len(rkeys) != len(set(rkeys)):
        raise ValueError("Duplicate new ratio condition")
    results = []
    for arm in ("common", "longest"):
        for view in F.VIEWS:
            for inclusive in (False, True):
                h1s, h2s = [], []
                for model in F.MODELS:
                    these = [mapped_ratio(r) for r in ratios if (r["arm"], r["view"], r["model"]) == (arm, view, model)]
                    ts = [r for r in tolerances if (r["arm"], r["view"], r["model"], F.boolean(r["forced_inclusive"])) == (arm, view, model, inclusive)]
                    if not these or len(ts) < 2:
                        raise ValueError(f"Missing designed ratio or tolerance coverage: {arm}/{view}/{model}/{inclusive}")
                    sizes = [float(r["C_label"]) for r in ts]
                    if len(sizes) != len(set(sizes)):
                        raise ValueError("Duplicate new tolerance condition")
                    small = next(r for r in ts if float(r["C_label"]) == min(sizes))
                    large = next(r for r in ts if float(r["C_label"]) == max(sizes))
                    h1, h2 = F.h1_model(these, small, large, inclusive=inclusive), F.h2_model(these, inclusive=inclusive)
                    h1s.append(h1)
                    h2s.append(h2)
                    results.extend(dict(arm=arm, view=view, forced_inclusive=inclusive, model=model, **r) for r in (h1, h2))
                previous_h2 = old_index[(arm, view, inclusive, "all", "H2_overall")]
                old_components = json.loads(previous_h2["components"])
                a6 = next(c for c in old_components if c["criterion"] == "H2_A6_predictor_gate")
                h3 = old_index[(arm, view, inclusive, "all", "H3_overall")]
                overall = [F.at_least("H1_overall", h1s, 3),
                           F.conjunction("H2_overall", [F.at_least("H2_model_count", h2s, 2), a6]),
                           {k: json.loads(v) if k in ("pairs", "components") and v else v
                            for k, v in h3.items() if k not in ("arm", "view", "forced_inclusive", "model") and v != ""}]
                results.extend(dict(arm=arm, view=view, forced_inclusive=inclusive, model="all", **r) for r in overall)
    if {key(r) for r in results} != expected or len(results) != 308:
        raise ValueError("New S0 lost a required forcing/window/arm criterion")
    paired = []
    for row in results:
        previous = old_index[key(row)]
        scope = ("legacy tau-based H3 retained; no new-rate H3 call" if row["criterion"] == "H3_overall" else
                 "updated rho_k model count with legacy direct-tau A6 gate; no new-rate predictor validation" if row["criterion"] == "H2_overall" else
                 "rho_k and fixed-lag fhalf substituted; unchanged frozen criteria")
        paired.append(dict(row, old_tau_status=previous["status"], changed=row["status"] != previous["status"],
                           old_tau_record=previous,
                           estimand_scope=scope))
    return paired


def old_replay(old, old_ratios, old_tolerances):
    """Check the reused arithmetic against all original H1/H2 calls."""
    ratios = []
    for source in old_ratios:
        row = dict(source)
        for suffix in ("", "_lower", "_upper", "_ci_lower", "_ci_upper", "_ci_upper_unbounded"):
            row["rho_k" + suffix] = source.get("tau_ratio" + suffix, "")
        ratios.append(row)
    rows = evaluate(old, ratios, old_tolerances)
    differences = [key(r) for r in rows if r["changed"]]
    if differences:
        raise ValueError(f"Frozen old-estimand replay differs: {differences}")
    return len(rows)


def render(rows, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import ListedColormap
    import numpy as np

    names = {"wave_pinning": "Wave-pinning", "goryachev": "Goryachev", "holmes_model3": "Holmes", "spring_ring": "Spring", "all": "All models"}
    statuses = {"fail": 0, "inconclusive": 1, "pass": 2}
    texts = {"fail": "not met", "inconclusive": "?", "pass": "met"}
    order = [(m, c) for m in F.MODELS for c in ("H1_model", "H2_model")] + [("all", c) for c in ("H1_overall", "H2_overall", "H3_overall")]
    vectors = out / "supporting/vectors"
    vectors.mkdir(parents=True, exist_ok=True)
    with PdfPages(out / "figure4_S0_round2_old_new.pdf") as book:
        for arm in ("common", "longest"):
            for inclusive in (False, True):
                selected = [r for r in rows if r["arm"] == arm and r["forced_inclusive"] == inclusive]
                lookup = {(r["model"], r["criterion"], r["view"]): r for r in selected}
                values = np.array([[statuses[lookup[m, c, v]["status"]] for v in F.VIEWS] for m, c in order])
                fig, ax = plt.subplots(figsize=(11.7, 8.3))
                ax.imshow(values, cmap=ListedColormap(("#efc1b8", "#ededed", "#b8d8c0")), vmin=0, vmax=2, aspect="auto")
                for i, (m, c) in enumerate(order):
                    for j, v in enumerate(F.VIEWS):
                        r = lookup[m, c, v]
                        ax.text(j, i, texts[r["old_tau_status"]] + " → " + texts[r["status"]], ha="center", va="center", fontsize=8)
                ax.set_xticks(range(7), [v.replace("analysis_", "").replace("_", " ") for v in F.VIEWS], rotation=35, ha="right")
                ax.set_yticks(range(11), [names[m] + " · " + c.replace("_model", "").replace("_overall", " overall") for m, c in order])
                ax.set_title(f"S0 · {arm} · {'with' if inclusive else 'without'} forced states\nOld persistence-time call → new rate/noise-tolerance call", loc="left", fontsize=12)
                fig.subplots_adjust(left=.24, right=.97, top=.87, bottom=.26)
                fig.text(.24, .10, "Colour shows new call. met: criterion met; not met: criterion not met; ?: inconclusive.\nFrozen criteria text is unchanged. H3 and the A6 predictor gate retain their explicitly direct-τ definitions.\nGeometric normalization is a reported outcome, not a gate. Full components and old records are in the CSV.", fontsize=9)
                book.savefig(fig)
                stem = f"S0_{arm}_{'inclusive' if inclusive else 'primary'}"
                fig.savefig(out / (stem + ".png"), dpi=160)
                fig.savefig(vectors / (stem + ".svg"))
                plt.close(fig)


def interim_summary(rows, ratios):
    main = [r for r in rows if r["view"] == "main" and r["forced_inclusive"] is False]
    h1 = [r for r in main if r["criterion"] == "H1_overall"]
    case = "A" if len(h1) == 2 and all(r["status"] == "pass" for r in h1) else "B"
    models = ("wave_pinning", "holmes_model3", "spring_ring")
    def example(model, f):
        return next(r for r in ratios if r["arm"] == "common" and r["view"] == "main" and r["model"] == model and float(r["f"]) == f)
    weak = [example(m, .003) for m in models]
    strong = [example(m, .05) for m in models]
    fmt = lambda rr, field: ", ".join(f"{float(r[field]):.2f}" for r in rr)
    statuses = "; ".join(r["arm"] + ": " + r["status"] for r in h1)
    paragraph = (
        f"Case {case}: the complete frozen H1 remains {h1[0]['status']} in the common and longest main-view arms "
        "without forced conditions. Its unchanged requirement is a ratio CI strictly above one at every evaluable input, "
        "no supported reversal, and non-overlapping size noise-tolerance CIs; weak-to-intermediate gains alone do not satisfy it. "
        f"At input 0.003, the rate ratios are {fmt(weak, 'rho_k')} for Wave-pinning, Holmes and Spring; at 0.05 they are "
        f"{fmt(strong, 'rho_k')}, with geometric-normalized ratios {fmt(strong, 'rho_geo')} and multiple-front majorities. "
        "Thus resolved size gains can coexist with erosion below the single-front geometric expectation, while Spring's "
        "0.05 ratio CI includes one. Main-view size noise-tolerance CIs separate in these three models. "
        "Goryachev remains coverage/lag/forcing qualified. Updated H2 model calls are reported, but H2 overall retains "
        "its legacy direct-time A6 predictor gate, and H3 retains its legacy time-bound definition; neither is a fully "
        "revalidated rate-based hypothesis. Stop for review before rebuilding the Figure 4 main."
    )
    return dict(case=case, summary=paragraph, frozen_H1_main=statuses,
                qualifications=[
                    "All 308 S0 rows cover both arms, all seven views, and both forcing choices; old calls appear alongside updated calls.",
                    "Frozen H1 applies to every evaluable input, with at least three of four models required overall; no weak-input-only gate was introduced.",
                    "The authorized fhalf substitution uses absolute heading autocorrelation 0.5 at the fixed reference lag, replacing the historical weakest-input-relative tolerance; the original authority text is retained verbatim as provenance.",
                    "H2 overall is a mixed-estimand replay with the frozen legacy direct-tau A6 predictor gate, not a new-rate predictor validation.",
                    "H3 is retained as its legacy tauStrong/tauWeak empirical-bound result; a new-rate H3 was not defined or validated.",
                    "Geometric normalization is relative to the single-front expectation, is not site diffusivity, and is not a decision gate.",
                    "Non-exponential lag dependence, unavailable exact saved-time lags, forced states, and one-sided bounds remain explicit in the metric tables.",
                    "The association between multiple fronts and diminished geometric advantage is descriptive, not an independently validated causal or predictive claim.",
                    "The modeled activity heading is distinct from observed migration summaries and conditional OU persistence estimates."
                ],
                pending=["User review of the Figure 4 interim (mandatory stop).", "Build/recompose the Figure 4 main only after review.",
                         "Figure 3 decomposition/recomposition and any conditional extension-run count remain subsequent work.",
                         "Figure 2 targeted Holmes retention runs, class verification, and recomposition remain subsequent work."])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--old", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    complete = json.loads((args.old / "complete.json").read_text())
    for name in ("criteria.csv", "size_ratios.csv", "noise_tolerance.csv"):
        if hashlib.sha256((args.old / name).read_bytes()).hexdigest() != complete["outputs"][name]:
            raise ValueError(f"Accepted original table changed: {name}")
    old = read(args.old / "criteria.csv")
    old_ratios, old_tolerances = read(args.old / "size_ratios.csv"), read(args.old / "noise_tolerance.csv")
    replayed = old_replay(old, old_ratios, old_tolerances)
    data = args.out / "supporting/data"
    ratios, tolerances = read(data / "size_ratios.csv"), read(data / "noise_tolerance.csv")
    rkey = lambda r: (r["arm"], r["model"], r["view"], float(r["f"]))
    tkey = lambda r: (r["arm"], r["model"], r["view"], float(r["C_label"]), F.boolean(r["forced_inclusive"]))
    if {rkey(r) for r in ratios} != {rkey(r) for r in old_ratios}:
        raise ValueError("New ratios omit or add designed input comparisons")
    if {tkey(r) for r in tolerances} != {tkey(r) for r in old_tolerances}:
        raise ValueError("New tolerance grid differs from original arms/views/sizes/forcing views")
    rows = evaluate(old, ratios, tolerances)
    write(data / "S0_old_new_criteria.csv", rows)
    write(data / "S0_old_tau_criteria.csv", old)
    components = [dict(arm=r["arm"], view=r["view"], forced_inclusive=r["forced_inclusive"], model=r["model"], parent=r["criterion"], **c)
                  for r in rows for c in r.get("components", [])]
    write(data / "S0_new_components.csv", components)
    authority = args.old.parents[1] / "inputs/20261020_figure4_final_statistics_decisions.json"
    metadata = args.out / "supporting/metadata"
    metadata.mkdir(parents=True, exist_ok=True)
    (metadata / "S0_unchanged_criteria_text.json").write_text(authority.read_text())
    render(rows, args.out)
    (metadata / "interim_summary.json").write_text(json.dumps(interim_summary(rows, ratios), indent=2) + "\n")
    tracked = [args.old / "criteria.csv", args.old / "size_ratios.csv", args.old / "noise_tolerance.csv", authority,
               data / "size_ratios.csv", data / "noise_tolerance.csv", Path(__file__), Path(F.__file__)]
    receipt = dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True, rows=len(rows),
                   old_estimand_exact_replay_rows=replayed, changed_rows=sum(r["changed"] for r in rows),
                   both_arms=True, all_seven_views=True, both_forcing_views=True, geometric_normalization_is_not_gate=True,
                   preserved_legacy_dependencies=["H2 overall retains frozen A6 direct-tau predictor gate", "H3 retains frozen tauStrong/tauWeak empirical-bound test"],
                   inputs={str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in tracked},
                   overall_main=[r for r in rows if r["view"] == "main" and r["model"] == "all" and r["criterion"] != "H3_overall"])
    (metadata / "S0_decision_replay_audit.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(encode({k: v for k, v in receipt.items() if k not in ("inputs", "overall_main")}))


if __name__ == "__main__":
    main()
