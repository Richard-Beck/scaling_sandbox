"""Display additive geometric qualifiers without changing hypothesis truth."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import signature_figure4_decisions_render_20261020 as D
import signature_figure4_interim_render_20261020 as P
import signature_figure4_package_20261020 as Q
import signature_figure4_stage_render_20261020 as S


def load_geometry(bundle, digest):
    bundle = Path(bundle).resolve(); root = bundle / "results/geometric_qualifiers"
    status_path = root / "status.json"
    Q.verified_file(dict(path=str(status_path), sha256=digest), root)
    status = json.loads(status_path.read_text())
    for key in ("passed", "complete", "compatibility_is_not_equivalence",
                "unbounded_intervals_do_not_establish_geometric_labels", "conflicting_time_mobility_conclusions_retained"):
        if status.get(key) is not True:
            raise ValueError("The complete separate geometric qualifier and its scientific qualifications are required")
    complete_path, audit_path = Path(status["source_complete_path"]), Path(status["source_independent_audit_path"])
    for path, key in ((complete_path, "source_complete_sha256"), (audit_path, "source_independent_audit_sha256")):
        Q.verified_file(dict(path=str(path), sha256=status[key]), bundle)
    complete, audit = (json.loads(path.read_text()) for path in (complete_path, audit_path))
    if (not complete.get("passed") or not complete.get("complete") or not audit.get("passed") or
            audit.get("producer_complete_sha256") != status["source_complete_sha256"] or
            audit.get("final_outputs_sha256") != complete["outputs"] or
            complete["source_freeze_sha256"] != status["source_final_freeze_sha256"]):
        raise ValueError("The geometry projection must bind the independently accepted unchanged final producer")
    active_path = bundle / "ACTIVE_FINAL_ANALYSIS.json"; active = json.loads(active_path.read_text())
    if active["active_source_freeze_sha256"] != status["source_final_freeze_sha256"]:
        raise ValueError("The geometry projection uses a stale final source closure")
    freeze_path = bundle / "geometry_projection_freeze.json"; freeze = json.loads(freeze_path.read_text())
    projection_path = root / Path(freeze["source"]).name
    for path in (projection_path, Path(freeze["source"])):
        Q.verified_file(dict(path=str(path), sha256=freeze["source_sha256"]), bundle)
    if freeze["source_sha256"] != status["projection_source_sha256"]:
        raise ValueError("The geometric projection differs from its reviewed separate source")
    ratio_path = Path(status["source_size_ratios_path"])
    Q.verified_file(dict(path=str(ratio_path), sha256=status["source_size_ratios_sha256"]), bundle)
    if complete["outputs"][ratio_path.name] != status["source_size_ratios_sha256"]:
        raise ValueError("The geometry projection substituted different raw ratio measurements")
    criteria_path = complete_path.parent / "criteria.csv"
    Q.verified_file(dict(path=str(criteria_path), sha256=complete["outputs"][criteria_path.name]), bundle)
    D.validate_criteria(P.read_csv(criteria_path))
    tables = {}
    sources = dict(geometry_status=dict(path=str(status_path), sha256=digest),
        projection_freeze=dict(path=str(freeze_path), sha256=P.sha(freeze_path)),
        unchanged_final_producer=dict(path=str(complete_path), sha256=status["source_complete_sha256"]),
        independent_final_decisions=dict(path=str(audit_path), sha256=status["source_independent_audit_sha256"]),
        original_size_ratios=dict(path=str(ratio_path), sha256=status["source_size_ratios_sha256"]))
    for name, expected in status["outputs"].items():
        path = root / name; Q.verified_file(dict(path=str(path), sha256=expected), bundle)
        tables[name] = P.read_csv(path); sources[name] = dict(path=str(path), sha256=expected)
    conditions, models = tables["geometric_conditions.csv"], tables["geometric_models.csv"]
    if len(conditions) != status["condition_rows"] or len(models) != status["model_rows"]:
        raise ValueError("The geometric qualifier table census changed")
    expected = {(arm, view, model, inclusive) for arm in D.ARMS for view in D.VIEWS for model in P.MODELS for inclusive in (False, True)}
    actual = [(row["arm"], row["view"], row["model"], P.boolean(row["forced_inclusive"])) for row in models]
    if len(actual) != len(set(actual)) or set(actual) != expected:
        raise ValueError("Every model needs both geometric forcing views and seven windows in both arms")
    return tables, sources


def ci_mark(ax, x, row, point, prefix, *, color, qualified=True):
    """Show actual points and finite CI endpoints without inventing midpoints."""
    value = P.num(row, point); lo, hi = P.num(row, prefix + "_ci_lower"), P.num(row, prefix + "_ci_upper")
    if P.boolean(row.get(prefix + "_ci_upper_unbounded")) is True:
        hi = math.inf
    if math.isfinite(value):
        ax.plot(x, value, "o", ms=3, color=color, mfc=color if qualified else "none")
    else:
        ax.plot(x, .025, "x", color=color, ms=3, transform=ax.get_xaxis_transform())
    if math.isfinite(lo) and math.isfinite(hi):
        ax.vlines(x, lo, hi, color=color, lw=.75)
    elif math.isfinite(lo):
        ax.plot(x, lo, "^", mfc="none", ms=4, color=color)
    elif math.isfinite(hi):
        ax.plot(x, hi, "v", mfc="none", ms=4, color=color)


def summary_pages(book, rows):
    for arm in D.ARMS:
        for inclusive in (False, True):
            selected = [row for row in rows if row["arm"] == arm and P.boolean(row["forced_inclusive"]) is inclusive]
            main = [row for row in selected if row["view"] == "main"]
            for start in (0, 2):
                batch = P.MODELS[start:start + 2]
                fig, axes = plt.subplots(2, 1, figsize=(11.7, 8.3))
                for ax, model in zip(axes, batch):
                    ax.axis("off"); row = next(row for row in main if row["model"] == model)
                    text = (P.LABELS[model] + ": " + D.words(row["per_model_pattern"]) + "\n" +
                        f"Included {row['included_inputs']} / {row['designed_inputs']} designed inputs; finite qualified time intervals {row['finite_time_intervals']}; "
                        f"finite qualified mobility intervals {row['finite_mobility_intervals']}; conflicts {row['conflicting_inputs']}.\n" +
                        "Time: " + D.words(D.nested(row["time_conclusion_counts"])) + ".\n" +
                        "Mobility: " + D.words(D.nested(row["mobility_conclusion_counts"])) + ".\n" +
                        "Input strengths by qualifier: " + D.words(D.nested(row["qualifier_input_strengths"])) + ".")
                    lines = "\n".join(textwrap.fill(line, 140) for line in text.splitlines())
                    ax.text(0, .98, lines, va="top", fontsize=8.5, linespacing=1.4)
                fig.subplots_adjust(left=.06, right=.97, top=.84, bottom=.22, hspace=.15)
                book.save(fig, f"S0_geometry_{arm}_{'inclusive' if inclusive else 'primary'}_{start//2+1}", "Per-model geometric qualification · " + arm + " · " + D.forcing_label(inclusive), [row for row in selected if row["model"] in batch],
                    "Main-view summaries are shown; source CSV retains all seven windows/cadences. Time excess uses direct persistence divided by actual C-ratio squared. Dx ratios convert short-lag global activity-heading fluctuations using actual C; they do not measure one tracked site's diffusion in multi-region states. Qualified single-region literal A3 heading/peak MSDs remain separate. A CI including flat Dx establishes compatibility only, never equivalence. Unbounded/undefined inference stays unavailable; increased Dx, angular persistence advantages, conflicts and heterogeneous inputs all remain visible. These labels do not change the 308 H1–H3 criteria.")


def condition_pages(book, rows):
    for arm in D.ARMS:
        for inclusive in (False, True):
            for view in D.VIEWS:
                selected = [row for row in rows if row["arm"] == arm and row["view"] == view and P.boolean(row["forced_inclusive"]) is inclusive]
                fig, axes = plt.subplots(4, 3, figsize=(11.7, 8.3), squeeze=False)
                for i, model in enumerate(P.MODELS):
                    group = sorted(P.subset(selected, which_model=model), key=lambda row: P.num(row, "f"))
                    for j, (point, prefix) in enumerate((("tau_ratio", "tau_ratio"), ("excess_persistence", "excess_persistence"), ("Dx_large_over_small", "Dx_ratio"))):
                        ax = axes[i, j]
                        for row in group:
                            qualified = row["mobility_geometry_status"] == "available" if j == 2 else row["time_geometry_status"] == "available"
                            color = "#888888" if row["geometric_qualification"] != "included" else "#222222"
                            ci_mark(ax, P.num(row, "f"), row, point, prefix, color=color, qualified=qualified)
                        if j == 0:
                            P.curve(ax, group, "f", ("actual_geometry_ratio",), color="#34974a", label="actual C ratio²")
                            ax.axhline(9., ls="--", color=".6", lw=.7, label="resting reference 9")
                        else:
                            ax.axhline(1., ls=":", color=".6", lw=.7)
                        if i == 0:
                            ax.set_title(("Direct persistence ratio", "Direct excess over geometry", "Global-heading Dx proxy ratio")[j], fontsize=8)
                        P.polish(ax, xlabel="Input fraction f" if i == 3 else None, ylabel=P.LABELS[model] if j == 0 else None, logx=True)
                        ax.set_yscale("symlog", linthresh=.5)
                        if i == 0 and j == 0:
                            ax.legend(frameon=False, fontsize=5.5)
                fig.subplots_adjust(left=.08, right=.97, top=.84, bottom=.24, hspace=.65, wspace=.35)
                book.save(fig, f"S4_S5_geometry_{arm}_{'inclusive' if inclusive else 'primary'}_{view}",
                    "Direct time and short-lag mobility · " + arm + " · " + D.words(view) + " · " + D.forcing_label(inclusive), selected,
                    "Columns keep direct time, direct excess-persistence CI and short-lag global-heading Dx-proxy CI separate. Dx converts global activity-heading fluctuations using actual C; in multi-region states it is not one tracked site's diffusion. Qualified single-region literal A3 heading/peak MSDs remain separate. Green is actual C-ratio²; grey dashed9 is the resting reference. Filled points have finite qualified inference; open circles retain raw points without it. Open triangles retain a finite confidence endpoint when the other is unbounded; × is unidentified, not zero. Symlog preserves signed raw estimates. CI including1 means compatibility, not equivalence. Increased Dx, persistence advantages and conflicts remain visible; this descriptor does not alter H1–H3.")


def render(bundle, digest, output, name):
    tables, sources = load_geometry(bundle, digest)
    sources["geometry_display_source"] = dict(path=str(Path(__file__).resolve()), sha256=P.sha(__file__))
    book = S.Book(output, name, "Separate geometric qualification", sources)
    try:
        summary_pages(book, tables["geometric_models.csv"])
        condition_pages(book, tables["geometric_conditions.csv"])
        return book.close()
    except Exception:
        book.pdf.close()
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--status-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--name", required=True)
    args = parser.parse_args()
    print(json.dumps(render(args.bundle, args.status_sha256, args.output, args.name)["pdf"], indent=2))
