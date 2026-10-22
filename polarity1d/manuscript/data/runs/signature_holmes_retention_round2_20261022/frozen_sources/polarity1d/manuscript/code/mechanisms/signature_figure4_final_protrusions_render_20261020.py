"""Render canonical common/longest A5 tables with exact native context.

The accepted E0 display closure is separate. All region estimators and seed
bootstrap outputs are read from completed independently audited A5 collectors.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import signature_figure4_interim_render_20261020 as P
import signature_figure4_package_20261020 as Q
import signature_figure4_protrusions_render_20261020 as A
import signature_figure4_stage_render_20261020 as S


def key(row):
    return (P.model(row), P.num(row, "C_label"), P.num(row, "f"), row["kind"],
            P.num(row, "sigma"), P.num(row, "native_horizon"))


def variant(row):
    return P.num(row, "threshold"), int(P.num(row, "downsample_factor")), P.boolean(row["smoothed"])


def validate_metrics(rows, expected_conditions):
    lookup = {}
    for row in rows:
        identity = key(row)
        if identity not in expected_conditions:
            raise ValueError("A5 source tables contain an observation outside the exact canonical inventory")
        native = P.num(row, "native_profile_cadence")
        if (not math.isfinite(native) or native <= 0 or
                not math.isclose(P.num(row, "observation_cadence"), native * P.num(row, "downsample_factor"))):
            raise ValueError("Native factors cannot be relabelled as a different physical observation cadence")
        item = identity, variant(row), row["statistic"]
        if item in lookup:
            raise ValueError("A5 canonical metric identities are duplicated")
        lookup[item] = row
        if row["kind"] != "noise" and (P.num(row, "independent_seeds") != 0 or P.num(row, "bootstrap_resamples") != 0):
            raise ValueError("Deterministic controls must have zero independent seeds and bootstrap draws")
    factorial = {(threshold, factor, smooth) for threshold in (.3, .5, .7)
                 for factor in (1, 2, 4, 8) for smooth in (False, True)}
    for identity in expected_conditions:
        selected = [row for (condition, _, _), row in lookup.items() if condition == identity]
        if not selected:
            raise ValueError("A5 canonical condition has no metrics")
        primary = [row for row in selected if P.boolean(row["primary"]) is True]
        if len(primary) != len(A.NAMES) or {row["statistic"] for row in primary} != set(A.NAMES):
            raise ValueError("Every A5 condition needs one complete physical-cadence primary")
        target = 8. if identity[0] == "spring_ring" else 80.
        if any(P.num(row, "threshold") != .5 or P.boolean(row["smoothed"]) is not True or
               not math.isclose(P.num(row, "observation_cadence"), target) for row in primary):
            raise ValueError("A5 primary must be joined by physical cadence, never a fixed native factor")
        declared = factorial | {variant(primary[0])}
        observed = {(variant(row), row["statistic"]) for row in selected}
        if observed != {(v, statistic) for v in declared for statistic in A.NAMES}:
            raise ValueError("All 24 native-factor variants and a distinct common-cadence primary must remain explicit")


def load_tables(run_root, complete_digest, native_digest, observation_root, observation_digest):
    root, observation_root = Path(run_root).resolve(), Path(observation_root).resolve()
    complete_path = root / "results/analysis/complete.json"
    native_path = root / "results/independent_all_views_audit.json"
    for path, digest in ((complete_path, complete_digest), (native_path, native_digest)):
        Q.verified_file(dict(path=str(path), sha256=digest), root)
    complete, native = (json.loads(path.read_text()) for path in (complete_path, native_path))
    inventory_path = root / "inventory.json"
    inventory = json.loads(inventory_path.read_text())
    cases = inventory["cases"]; ids = {case["id"] for case in cases}
    if inventory.get("analysis_arm") not in ("common", "longest") or len(cases) != len(ids):
        raise ValueError("Use a declared complete common or longest canonical A5 inventory")
    if any(P.num(case, "sigma") != (1. if case["kind"] == "noise" else 0.) or
           not math.isfinite(P.num(case, "native_horizon")) or P.num(case, "native_horizon") <= 0 for case in cases):
        raise ValueError("Canonical common/longest observations must retain sigma-one noise and exact native horizons")
    if (complete.get("passed") is not True or native.get("passed") is not True or native.get("complete") is not True or
            complete.get("inventory_sha256") != P.sha(inventory_path) or native.get("inventory_sha256") != P.sha(inventory_path) or
            complete.get("condition_identity_checked_against_inventory") is not True or
            complete.get("cases") != len(cases) or native.get("native_cases") != len(cases) or
            {row["case_id"] for row in native["records"]} != ids or
            any(row.get("passed") is not True or row.get("views_checked", 0) < 24 for row in native["records"])):
        raise ValueError("Final A5 canonical native/statistical coverage is incomplete")
    sources = dict(A5_complete=dict(path=str(complete_path), sha256=complete_digest),
                   A5_all_native_views=dict(path=str(native_path), sha256=native_digest),
                   canonical_inventory=dict(path=str(inventory_path), sha256=P.sha(inventory_path)))
    for row in native["source_shards"]:
        Q.verified_file(row, root)
    primary_path = root / "results/independent_primary_native_frames_audit.json"
    primary = json.loads(primary_path.read_text())
    if (primary.get("passed") is not True or primary.get("complete") is not True or
            primary.get("inventory_sha256") != P.sha(inventory_path) or primary.get("subset_of_all_views") is not True or
            primary.get("all_views_audit_sha256") != native_digest):
        raise ValueError("Primary-frame projection must bind the same accepted all-view native audit")
    sources["A5_primary_projection"] = dict(path=str(primary_path), sha256=P.sha(primary_path))
    worker_path = root / "inputs/worker_freeze.json"; worker = json.loads(worker_path.read_text())
    if worker["inventory_sha256"] != P.sha(inventory_path):
        raise ValueError("A5 worker closure covers a different canonical inventory")
    for name, digest in worker["source_sha256"].items():
        Q.verified_file(dict(path=str(root / "inputs/frozen_sources" / name), sha256=digest), root)
    sources["A5_worker_closure"] = dict(path=str(worker_path), sha256=P.sha(worker_path))
    Q.verified_file(complete["seed_receipts"], root)
    seed_receipts = json.loads(Path(complete["seed_receipts"]["path"]).read_text())
    if {row["id"] for row in seed_receipts} != ids:
        raise ValueError("A5 collector seed receipts omit canonical cases")
    tables = {}
    for name in ("metrics", "distributions", "survival"):
        source = complete["outputs"][name]; path = Q.verified_file(source, root)
        sources[name] = source
        selected, row_count = [], 0
        with path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                row_count += 1
                if row["experiment"] != inventory["analysis_arm"]:
                    raise ValueError("The A5 table declared observation arm changed")
                if name == "metrics" or P.boolean(row["primary"]) is True:
                    selected.append(row)
        if row_count != source["rows"]:
            raise ValueError("The A5 table census or declared observation arm changed")
        tables[name] = selected
    expected_conditions = {key(case) for case in cases}
    if complete["groups"] != len(expected_conditions):
        raise ValueError("A5 canonical group census changed")
    validate_metrics(tables["metrics"], expected_conditions)
    status_path = observation_root / "status.json"
    Q.verified_file(dict(path=str(status_path), sha256=observation_digest), observation_root)
    status = json.loads(status_path.read_text())
    if not status.get("passed") or not status.get("complete") or status.get("analysis_arm") != inventory["analysis_arm"]:
        raise ValueError("A5 display qualifications require the same frozen observation arm")
    heading_path = observation_root / "condition_metrics.csv"
    Q.verified_file(dict(path=str(heading_path), sha256=status["outputs"][heading_path.name]), observation_root)
    headings = {key(dict(row, kind="noise", sigma=1.)): row for row in P.main_rows(P.read_csv(heading_path))}
    for rows in tables.values():
        for row in rows:
            horizon=P.num(row,"native_horizon")
            row.update(display_analysis_start=horizon/2,display_analysis_end=horizon,
                       display_analysis_duration=horizon/2)
            if row["kind"] == "noise":
                if key(row) not in headings:
                    raise ValueError("A5 condition lacks the exact same-horizon heading/forcing qualification")
                heading = headings[key(row)]
                row["display_forced"] = heading.get("forced_combined")
                row["display_forced_definition"] = "forced_combined"
    sources["frozen_observation_status"] = dict(path=str(status_path), sha256=observation_digest)
    sources["heading_forcing_qualifications"] = dict(path=str(heading_path), sha256=P.sha(heading_path))
    return tables, sources, inventory["analysis_arm"]


class Book(S.Book):
    def save(self, fig, page, title, source_rows, note):
        fig.subplots_adjust(bottom=max(fig.subplotpars.bottom,.24),top=min(fig.subplotpars.top,.88))
        note = note.replace("Solid/dashed/dotted also distinguish small/intermediate/large sizes.",
                            "Size colours and line patterns distinguish all declared sizes.")
        note = note.replace("Solid/dashed/dotted distinguish small/intermediate/large;", "Size colours and line patterns distinguish sizes;")
        note = note.replace("Native coarsely saved profiles limit birth/turn observability.",
                            "Finite native and imaging cadences limit birth/turn observability; exact cadences remain in CSV.")
        horizons={model:sorted({P.num(row,"native_horizon") for row in source_rows if P.model(row)==model}) for model in P.MODELS}
        horizon_text="; ".join(P.LABELS[model]+" "+"/".join(f"{value:g}" for value in values)
                              for model,values in horizons.items() if values)
        note+=" Native H: "+horizon_text+"; analysis window [H/2,H]. Event counts per trajectory cover that window; no cross-model rate ranking. Actual independent seed counts are condition-specific; some weakest-input groups have two, and controls/companions add no seeds."
        super().save(fig, "S9_" + page, title, source_rows, note)


def metric_series_by_horizon(ax, rows, *, ls="-", label=None, show_mean=False, count=False):
    """Never connect event counts collected over different native durations."""
    horizons=sorted({P.num(row,"native_horizon") for row in rows if row["kind"]=="noise"})
    for i,horizon in enumerate(horizons):
        selected=[row for row in rows if P.num(row,"native_horizon")==horizon]
        A.metric_series(ax,selected,{},ls=ls,label=label,show_mean=show_mean,
                        point_marker=("o","s","D","^")[i%4],
                        label_suffix=f", H={horizon:g}" if count else "")


def metric_pages(book, tables):
    rows = tables["metrics"]
    for statistic in sorted(A.NAMES):
        all_rows = [row for row in rows if row["statistic"] == statistic]
        primary = [row for row in all_rows if P.boolean(row["primary"]) is True]
        fig, axes = plt.subplots(4, 1, figsize=(11.7, 8.3))
        for ax, model in zip(axes, P.MODELS):
            metric_series_by_horizon(ax, P.subset(primary, which_model=model), show_mean=True,
                                     count=statistic in ("observed_births","observed_takeovers","observed_pairwise_takeovers","unscorable_takeover_tracks"))
            P.polish(ax, xlabel="Input fraction f" if model == P.MODELS[-1] else None,
                     ylabel=P.LABELS[model], logx=True)
            ax.legend(fontsize=6, frameon=False, ncol=5)
        fig.subplots_adjust(left=.09, right=.97, top=.91, bottom=.17, hspace=.5)
        book.save(fig, statistic, A.NAMES[statistic], primary,
                  "Imaging primary: half contour, fixed physical Gaussian 1.5 (Spring 0.15), physical cadence 80 (Spring 8), joined independently of native save factor. Grey marks retain operationally forced conditions; unknown flags stay explicit in CSV. Dashed matched-mean and dotted zero controls have no seed bootstrap. Complete native peak basins define positive active-signal excess, not molecule counts. Observed dominant-newborn takeover requires the newborn to be largest among qualified regions and pairwise D>0.5 against the previous dominant for two frames; this is not 75% of whole-cell activity with three or more regions. Lifetimes retain censoring; RMST restriction remains 16000 (Spring 1600) even when long survival curves extend farther.")
        for smooth in (False, True):
            fig, axes = plt.subplots(4, 3, figsize=(12.5, 9.2), squeeze=False)
            source = []
            for i, model in enumerate(P.MODELS):
                for j, threshold in enumerate((.3, .5, .7)):
                    ax = axes[i, j]
                    selected = [row for row in P.subset(all_rows, which_model=model)
                                if P.num(row, "threshold") == threshold and P.boolean(row["smoothed"]) is smooth
                                and int(P.num(row, "downsample_factor")) in (1, 2, 4, 8)]
                    for factor, style in zip((1, 2, 4, 8), ("-", "--", ":", "-.")):
                        metric_series_by_horizon(ax, [row for row in selected if int(P.num(row, "downsample_factor")) == factor],
                                                ls=style, label=f"native ×{factor}")
                    if i == 0:
                        ax.set_title(f"Contour {threshold:g} · " + ("physical Gaussian" if smooth else "raw optics"), fontsize=8)
                    P.polish(ax, xlabel="Input fraction f" if i == 3 else None,
                             ylabel=P.LABELS[model] if j == 0 else None, logx=True)
                    source.extend(selected)
            fig.subplots_adjust(left=.07, right=.97, top=.91, bottom=.17, hspace=.6, wspace=.3)
            book.save(fig, statistic + ("_smoothed_factorial" if smooth else "_raw_factorial"), A.NAMES[statistic] + " · full native-factor sensitivity", source,
                      "All three contours and both optics views are shown across these paired pages. Solid/dashed/dotted/dash-dot curves are native save factors 1/2/4/8; size colours encode actual C. Native profile cadences differ across source arms, so equal factors are not equal physical cadences. Every physical cadence is retained in the CSV. The separately displayed common-cadence primary may be a 25th variant, outside this 24-way factorial; no missing factor16 contour/optics sensitivity is invented. CIs resample whole independent seeds, not frames.")


def render(root, complete_digest, native_digest, observations, observation_digest, output, name):
    tables, sources, arm = load_tables(root, complete_digest, native_digest, observations, observation_digest)
    sources["final_A5_display_source"] = dict(path=str(Path(__file__).resolve()), sha256=P.sha(__file__))
    sources["shared_A5_plot_source"] = dict(path=str(Path(A.__file__).resolve()), sha256=P.sha(A.__file__))
    book = Book(output, name, arm.title() + " canonical A5 observations", sources, section="S9")
    try:
        metric_pages(book, tables)
        A.distribution_pages(book, tables, {})
        A.survival_pages(book, tables, {})
        return book.close()
    except Exception:
        book.pdf.close()
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--complete-sha256", required=True)
    parser.add_argument("--native-audit-sha256", required=True)
    parser.add_argument("--observation-root", type=Path, required=True)
    parser.add_argument("--observation-status-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--name", required=True)
    args = parser.parse_args()
    report = render(args.run_root, args.complete_sha256, args.native_audit_sha256,
                    args.observation_root, args.observation_status_sha256, args.output, args.name)
    print(json.dumps(report["pdf"], indent=2))
