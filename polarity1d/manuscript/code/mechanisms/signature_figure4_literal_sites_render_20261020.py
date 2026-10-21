"""Display literal material-origin A3 coordinates from audited native tables.

The inherited midpoint-increment and track-reanchored estimators remain in
their immutable source displays as sensitivities. No coordinate, estimand,
bootstrap, native integration or hypothesis decision is computed here.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import signature_figure4_interim_render_20261020 as P
import signature_figure4_package_20261020 as Q
import signature_figure4_stage_render_20261020 as S

METHODS = ("heading_literal_actualC", "profile_peak_literal_actualC")
POLICIES = ("split", "plus_pi", "minus_pi")
VIEWS = ("main", "analysis_first_half", "analysis_second_half", "late_quarter", "cadence_x2", "cadence_x4", "cadence_x8")


def identity(row):
    return P.model(row), P.num(row, "C_label"), P.num(row, "f"), P.num(row, "sigma"), P.num(row, "native_horizon")


def load_literal(root, digest):
    root = Path(root).resolve(); status_path = root / "results/status.json"
    Q.verified_file(dict(path=str(status_path), sha256=digest), root)
    status = json.loads(status_path.read_text())
    manifest_path = root / "manifest.json"; manifest = json.loads(manifest_path.read_text())
    if (status.get("passed") is not True or status.get("complete") is not True or
            status.get("manifest_sha256") != P.sha(manifest_path) or
            status.get("condition_groups") != len(manifest["groups"]) or
            status.get("no_native_integrations") is not True or status.get("no_H1_H2_H3_estimand_change") is not True or
            set(status.get("methods", [])) != set(METHODS) or set(manifest["views"]) != set(VIEWS) or
            status.get("source_version") != 2 or manifest.get("source_version") != 2 or
            set(status.get("antipodal_policies", [])) != set(POLICIES) or set(manifest.get("antipodal_policies", [])) != set(POLICIES) or
            status.get("antipodal_atol_radians") != 1e-12 or manifest.get("antipodal_atol_radians") != 1e-12 or
            status.get("all_qualified_native_frames_preserved") is not True or
            status.get("route_choices_are_not_independent_seeds") is not True):
        raise ValueError("Literal A3 displays require complete audited additive native measurements")
    for path, expected in {**manifest["source_sha256"], **manifest["input_receipts"]}.items():
        Q.verified_file(dict(path=path, sha256=expected), root)
    groups = {int(group["task"]): group for group in manifest["groups"]}
    tasks = [int(row["task"]) for row in status["receipts"]]
    if len(groups) != len(manifest["groups"]) or len(tasks) != len(set(tasks)) or set(tasks) != set(groups):
        raise ValueError("Every literal native group needs its own independent receipt")
    sources = dict(literal_status=dict(path=str(status_path), sha256=digest),
                   literal_manifest=dict(path=str(manifest_path), sha256=P.sha(manifest_path)))
    results = {}
    for receipt in status["receipts"]:
        task = int(receipt["task"]); folder = root / f"workers/{task:04d}"
        result_path, audit_path = folder / "result.json", root / f"audits/{task:04d}.json"
        for path, key in ((result_path, "result_sha256"), (audit_path, "audit_sha256")):
            Q.verified_file(dict(path=str(path), sha256=receipt[key]), root)
        result, audit = (json.loads(path.read_text()) for path in (result_path, audit_path))
        if (result.get("passed") is not True or audit.get("passed") is not True or
                result.get("task") != task or audit.get("task") != task or audit.get("source_version") != 2 or
                result.get("manifest_sha256") != P.sha(manifest_path) or audit.get("manifest_sha256") != P.sha(manifest_path) or
                audit.get("result_sha256") != receipt["result_sha256"] or
                audit.get("arrays_sha256") != result.get("arrays_sha256") or
                {(row["id"], row["scalar_sha256"], row["profile_sha256"]) for row in result["cases"]} !=
                    {(row["id"], row["scalar_sha256"], row["profile_sha256"]) for row in groups[task]["noise"]} or
                result.get("whole_seeds") != len(groups[task]["noise"]) or result.get("bootstrap_resamples") != 2000):
            raise ValueError("Literal A3 result/audit identities differ from the frozen native group")
        Q.verified_file(dict(path=str(folder / "arrays.npz"), sha256=result["arrays_sha256"]), root)
        results[task] = result
    expected_groups = {}
    for group in groups.values():
        row = dict(group, native_horizon=group["noise"][0]["native_horizon"])
        current = identity(row)
        if current in expected_groups:
            raise ValueError("Literal group identities are duplicated")
        expected_groups[current] = group
    path = root / "results/site_curves.csv.gz"
    Q.verified_file(dict(path=str(path), sha256=status["outputs"][path.name]), root)
    sources["literal_site_curves"] = dict(path=str(path), sha256=status["outputs"][path.name])
    qualifiers = {}
    expected_ambiguities = {}
    for current, group in expected_groups.items():
        result = results[int(group["task"])]
        expected = {(view, method, policy) for view in VIEWS for method in METHODS for policy in POLICIES}
        actual = {(row["view"], row["method"], row["antipodal_policy"]) for row in result["qualifications"]}
        if len(result["qualifications"]) != len(expected) or actual != expected:
            raise ValueError("Literal ambiguity policies lack complete audited qualifications")
        for row in result["qualifications"]:
            qualifiers[(current, row["view"], row["method"], row["antipodal_policy"])] = row
        cases = {row["id"] for row in group["noise"]}
        expected_cases = {(case, view, method) for case in cases for view in VIEWS for method in METHODS}
        actual_cases = {(row["id"], row["view"], row["method"]) for row in result["case_ambiguities"]}
        if len(result["case_ambiguities"]) != len(expected_cases) or actual_cases != expected_cases:
            raise ValueError("Every native case/window/method needs an ambiguity record")
        for row in result["case_ambiguities"]:
            expected_ambiguities[(current, row["id"], row["view"], row["method"])] = row
    rows, count, seen = [], 0, set()
    with gzip.open(path, "rt", newline="") as stream:
        for row in csv.DictReader(stream):
            count += 1; current = identity(row)
            if current not in expected_groups:
                raise ValueError("A literal curve changed its native group or horizon")
            group = expected_groups[current]
            if (row["native_group_sha256"] != group["native_group_sha256"] or
                    row["method"] not in METHODS or row["view"] not in VIEWS or row["antipodal_policy"] not in POLICIES or
                    P.boolean(row["primary_A3"]) is not (row["method"] == METHODS[0] and row["antipodal_policy"] == "split") or
                    P.num(row, "n_seeds") != len(group["noise"]) or P.num(row, "bootstrap_resamples") != 2000 or
                    not math.isfinite(P.num(row, "cadence")) or P.num(row, "cadence") <= 0):
                raise ValueError("A literal curve lost primary/sensitivity or whole-seed provenance")
            qual = qualifiers[(current, row["view"], row["method"], row["antipodal_policy"])]
            if any(P.num(row, name) != float(qual[name]) for name in
                   ("cadence", "antipodal_steps", "seeds_with_antipodal_steps", "qualified_frames")):
                raise ValueError("A literal curve changed its audited ambiguity qualification")
            seen.add((current, row["view"], row["method"], row["antipodal_policy"]))
            if row["view"] == "main":
                rows.append(row)
    if count != status["rows"] or seen != {(key, view, method, policy) for key in expected_groups for view in VIEWS for method in METHODS for policy in POLICIES}:
        raise ValueError("A declared literal method/window curve is missing")
    ambiguity_path = root / "results/case_ambiguities.csv.gz"
    Q.verified_file(dict(path=str(ambiguity_path), sha256=status["outputs"][ambiguity_path.name]), root)
    sources["literal_case_ambiguities"] = dict(path=str(ambiguity_path), sha256=status["outputs"][ambiguity_path.name])
    ambiguities, seen_cases = [], set()
    with gzip.open(ambiguity_path, "rt", newline="") as stream:
        for row in csv.DictReader(stream):
            current = identity(row); key = current, row["id"], row["view"], row["method"]
            if key not in expected_ambiguities or key in seen_cases:
                raise ValueError("A native ambiguity case is missing, duplicated or relabelled")
            native = expected_ambiguities[key]
            if row["native_group_sha256"] != expected_groups[current]["native_group_sha256"] or any(
                    P.num(row, name) != float(value) for name, value in native.items() if name not in ("id", "view", "method")):
                raise ValueError("An ambiguity sidecar changed an audited native case value")
            seen_cases.add(key); ambiguities.append(row)
    if seen_cases != set(expected_ambiguities):
        raise ValueError("A native case/window ambiguity sidecar is incomplete")
    return rows, ambiguities, sources


def literal_pages(book, rows):
    for sigma in sorted({P.num(row, "sigma") for row in rows}):
        current = [row for row in rows if P.num(row, "sigma") == sigma]
        strengths = sorted({P.num(row, "f") for row in current})
        for start in range(0, len(strengths), 3):
            batch = strengths[start:start + 3]
            for method, kind in tuple((method, kind) for kind in ("MSD", "coverage", "routes") for method in METHODS):
                fig, axes = plt.subplots(4, 3, figsize=(11.7, 8.3), squeeze=False); source = []
                for i, model in enumerate(P.MODELS):
                    for j in range(3):
                        ax = axes[i, j]
                        if j >= len(batch):
                            ax.axis("off"); continue
                        selected = P.subset(current, which_model=model, strength=batch[j])
                        selected = [row for row in selected if row["method"] == method and
                                    (kind == "routes" or row["antipodal_policy"] == "split")]
                        if not selected:
                            ax.text(.5, .5, "Outside declared grid", transform=ax.transAxes, ha="center", fontsize=7)
                        for size in sorted({P.num(row, "C_label") for row in selected}):
                            by_size = P.subset(selected, size=size)
                            color = P.size_color(model, size, P.forced(by_size[0]) is True)
                            cadences = {P.num(row, "cadence") for row in by_size}
                            if len(cadences) != 1:
                                raise ValueError("A displayed native-method curve has inconsistent observation cadence")
                            ties = {P.num(row, "antipodal_steps") for row in by_size}
                            if len(ties) != 1:
                                raise ValueError("A displayed native curve changes its antipodal step count")
                            label=f"C {size:g}, Δt={next(iter(cadences)):g}, ties={next(iter(ties)):g}"
                            series = [row for row in by_size if row["antipodal_policy"] == "split"]
                            if kind == "coverage":
                                P.curve(ax, series, "lag", ("coverage",), color=color, label=label)
                                P.curve(ax, series, "lag", ("coverage_lost_to_ambiguity_split",), color=color, label=None, ls=":")
                            elif kind == "routes":
                                for policy, style in zip(POLICIES, ("-", "--", ":")):
                                    route = [row for row in by_size if row["antipodal_policy"] == policy]
                                    P.curve(ax, route, "lag", ("MSD",), color=color, ls=style,
                                            label=label if policy == "split" else None)
                            else:
                                P.curve(ax, series, "lag", ("MSD",), color=color, label=label,
                                        lowkeys=("MSD_ci_lower",), highkeys=("MSD_ci_upper",))
                                if not any(math.isfinite(P.num(row, "MSD")) for row in series):
                                    ax.text(.02, .04 + .08 * sorted({P.num(row, 'C_label') for row in selected}).index(size),
                                                f"C {size:g}: no eligible pairs", transform=ax.transAxes, fontsize=5.5, color=color)
                        ax.set_title(f"{P.LABELS[model]} · f={batch[j]:g}", fontsize=8)
                        if kind == "coverage":
                            ax.set_ylim(-.02, 1.03)
                        P.polish(ax, xlabel="Lag (native time)" if i == 3 else None,
                                 ylabel=("Eligible / lost\npair fraction" if kind == "coverage" else "Physical position\nMSD") if j == 0 else None)
                        if ax.get_legend_handles_labels()[0]:
                            ax.legend(frameon=False, fontsize=4.8, ncol=2, handlelength=1.4, labelspacing=.15, columnspacing=.6)
                        source.extend(selected)
                title = ("Literal heading position · A3 primary split" if method == METHODS[0] else
                         "Literal profile-maximum position · sensitivity")
                if kind == "coverage":
                    title += " · pair coverage and ambiguity loss"
                elif kind == "routes":
                    title = title.replace("A3 primary split", "paired route sensitivities") + " · split / +π / −π"
                fig.subplots_adjust(left=.08, right=.94, top=.89, bottom=.25, hspace=.65, wspace=.36)
                book.save(fig, f"S5_literal_sigma{sigma:g}_{method}_{kind}_{start//3+1}", title + f" · σ={sigma:g}", source,
                    "Literal x=θC(t)/(2π) uses native material origin q=0 and section initial angle in [0,2π). Primary tracks split at unresolved antipodal steps (tolerance 10⁻¹² rad); both native frames remain. Native M₁ phase uses R=|M₁|/|M₁,basal|, contrast and single-region qualification. Solid lines split tracks; dotted pair coverage shows ambiguity loss. Route panels compare solid split, dashed +π and dotted −π: paired deterministic sensitivities, not measured directions, independent seeds or uncertainty bounds. Legends give actual cadence Δt and tie count; peak sampling may equal scalar cadence. Shading is a 2000-whole-seed interval; unavailable lags remain gaps. Grey marks qualified forced conditions when canonical flags are supplied; checkpoints remain descriptive. These qualified single-region site MSDs differ from the global-heading-derived Dθ→Dx proxy, which is not one tracked site's diffusion in multi-region states. No migration law is assumed.")


def ambiguity_pages(book, rows):
    """Expose native tie counts in all declared windows and cadence checks."""
    for sigma in sorted({P.num(row, "sigma") for row in rows}):
        for view in VIEWS:
            current = [row for row in rows if P.num(row, "sigma") == sigma and row["view"] == view]
            fig, axes = plt.subplots(4, 2, figsize=(11.7, 8.3), squeeze=False)
            source = []
            for i, model in enumerate(P.MODELS):
                for j, method in enumerate(METHODS):
                    ax = axes[i, j]
                    selected = [row for row in current if P.model(row) == model and row["method"] == method]
                    largest_count = 0.
                    for size in sorted({P.num(row, "C_label") for row in selected}):
                        by_size = P.subset(selected, size=size); points = []
                        for strength in sorted({P.num(row, "f") for row in by_size}):
                            cases = P.subset(by_size, strength=strength)
                            points.append(dict(f=strength, observed_tied_steps=sum(P.num(row, "antipodal_steps") for row in cases)))
                        largest_count = max(largest_count, max((row["observed_tied_steps"] for row in points), default=0.))
                        P.curve(ax, points, "f", ("observed_tied_steps",), color=P.size_color(model, size), label=f"C {size:g}")
                        ax.scatter([row["f"] for row in points], [row["observed_tied_steps"] for row in points], color=P.size_color(model, size), s=9)
                    ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1); ax.set_ylim(0., max(1., largest_count * 1.1))
                    ax.set_title(P.LABELS[model] + (" · native M1 heading" if j == 0 else " · profile peak"), fontsize=8)
                    P.polish(ax, xlabel="Mean input fraction f" if i == 3 else None,
                             ylabel="Antipodal\nstep count" if j == 0 else None)
                    if ax.get_legend_handles_labels()[0]: ax.legend(frameon=False, fontsize=5.5, ncol=3)
                    source.extend(selected)
            fig.subplots_adjust(left=.08, right=.94, top=.89, bottom=.23, hspace=.7, wspace=.3)
            book.save(fig, f"S5_ambiguity_sigma{sigma:g}_{view}", f"Native antipodal ambiguity · {view.replace('_', ' ')} · σ={sigma:g}", source,
                "Counts sum unresolved sampled antipodal edges across the existing whole-seed ensemble; stage contexts can overlap and do not add independent seeds. All seven declared views are retained. Primary literal coordinates split track labels at each tie, preserving both endpoint frames; eligible pair loss is shown in the paired coverage pages and retained at every lag in the source tables. Plus/minus-π routes are paired deterministic sensitivities. Exact native cadence, qualified/total frames, adjacent steps, and track counts before/after splitting remain in each case's page CSV. Zero means no observed tied step at this sampling, not continuous-path identifiability.")


def select_literal_observations(metrics, by_stage, ambiguities_by_stage):
    """Use exactly the frozen native observations, never overlapping cohorts."""
    references = {S.condition_key(row): row for row in metrics}
    if len(references) != len(metrics):
        raise ValueError("Canonical literal observations contain duplicate condition/view identities")
    curves = S.select_observation_curves(P.main_rows(metrics), by_stage)
    ambiguities = S.select_observation_curves(metrics, ambiguities_by_stage)
    for row in curves + ambiguities:
        reference = S.annotation(references[S.condition_key(row)])
        for key in ("display_forced", "display_forced_definition", "selected_observation_arm", "common_observation_arm"):
            if key in reference:
                row[key] = reference[key]
    return curves, ambiguities


def load_canonical(observation_root, digest, input_spec, literal_spec):
    root = Path(observation_root).resolve(); input_spec = Path(input_spec).resolve(); literal_spec = Path(literal_spec).resolve()
    status_path = root / "status.json"
    Q.verified_file(dict(path=str(status_path), sha256=digest), root)
    status = json.loads(status_path.read_text())
    if (status.get("passed") is not True or status.get("complete") is not True or
            status.get("analysis_arm") not in ("common", "longest") or status.get("sigma") != 1. or status.get("forcing_sigma") != 1. or
            status.get("input_spec_sha256") != P.sha(input_spec)):
        raise ValueError("Canonical literal displays require the exact completed common/longest observation gate")
    metrics_path = root / "condition_metrics.csv"
    Q.verified_file(dict(path=str(metrics_path), sha256=status["outputs"][metrics_path.name]), root)
    metrics = P.read_csv(metrics_path)
    if not metrics or any(row["analysis_arm"] != status["analysis_arm"] for row in metrics):
        raise ValueError("Canonical literal displays cannot pool observation arms")
    for records in status.get("source_stage_gates", {}).values():
        for receipt in records:
            Q.verified_file(receipt, root)
    observations, declared = json.loads(input_spec.read_text()), json.loads(literal_spec.read_text())
    stages = {"E4" if row["selected_observation_arm"] == "preselected_E4" else
              "E3" if row.get("common_observation_arm", row["selected_observation_arm"]) == "common_new_E3_size" else "E2" for row in metrics}
    if set(declared["stages"]) != stages:
        raise ValueError("The literal source closure must contain exactly the selected native stages")
    curves, ambiguities = {}, {}
    sources = dict(observation_status=dict(path=str(status_path), sha256=digest),
                   observation_metrics=dict(path=str(metrics_path), sha256=P.sha(metrics_path)),
                   selected_observation_inputs=dict(path=str(input_spec), sha256=P.sha(input_spec)),
                   literal_source_spec=dict(path=str(literal_spec), sha256=P.sha(literal_spec)))
    for stage in sorted(stages):
        claim = declared["stages"][stage]; native = observations["stages"][stage]
        literal_root = Path(claim["root"]).resolve(); manifest = json.loads((literal_root / "manifest.json").read_text())
        parent = Path(native["root"]).resolve(); parent_status = parent / "results/analysis/status.json"
        if (Path(manifest["parent"]).resolve() != parent or
                manifest["input_receipts"].get(str(parent_status)) != native["status_sha256"]):
            raise ValueError("A literal measurement belongs to different selected native observations")
        curves[stage], ambiguities[stage], sources[stage] = load_literal(literal_root, claim["status_sha256"])
    selected, selected_ambiguities = select_literal_observations(metrics, curves, ambiguities)
    return selected, selected_ambiguities, sources, status["analysis_arm"] + (" checkpoint" if status.get("checkpoint_only") else "")


def render_rows(rows, ambiguities, sources, output, name, context):
    sources["literal_display_source"] = dict(path=str(Path(__file__).resolve()), sha256=P.sha(__file__))
    book = S.Book(output, name, context + " literal material-origin physical-site measurements", sources, section="S5")
    try:
        literal_pages(book, rows)
        ambiguity_pages(book, ambiguities)
        return book.close()
    except Exception:
        book.pdf.close()
        raise


def render(root, digest, output, name):
    rows, ambiguities, sources = load_literal(root, digest)
    return render_rows(rows, ambiguities, sources, output, name, "/".join(sorted({row["source_stage"] for row in rows})))


def render_canonical(observation_root, digest, input_spec, literal_spec, output, name):
    rows, ambiguities, sources, arm = load_canonical(observation_root, digest, input_spec, literal_spec)
    return render_rows(rows, ambiguities, sources, output, name, arm.capitalize() + " selected observations")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--status-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--name", required=True)
    args = parser.parse_args()
    print(json.dumps(render(args.run_root, args.status_sha256, args.output, args.name)["pdf"], indent=2))
