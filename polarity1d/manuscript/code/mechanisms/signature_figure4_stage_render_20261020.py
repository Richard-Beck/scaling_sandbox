"""Separate S7/S8 and common/long Figure 4 displays from audited tidy exports.

No native integrations, new estimators, framing choice, or implicit horizon
pooling occurs here. Existing accepted E0 exports are never overwritten.
"""
from __future__ import annotations

import argparse
import csv
import gzip
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import LogLocator, NullFormatter

import signature_figure4_interim_render_20261020 as P
import signature_figure4_package_20261020 as Q
verified_file=Q.verified_file


def condition_key(row):
    return P.model(row), P.num(row, "C_label"), P.num(row, "f"), row.get("view", "main")


def annotation(row, flag="forced_combined"):
    result = dict(row)
    # Keep the original scientific forced column and all component flags.
    # The display-only column makes the chosen operational flag explicit.
    if flag in row:
        result["display_forced"] = row[flag]
        result["display_forced_definition"] = flag
    else:
        result["display_forced"] = row.get("forced")
        result["display_forced_definition"] = "matched_mean_only"
    return result


def stream_rows(path, *, view=None):
    with Path(path).open(newline="") as stream:
        for row in csv.DictReader(stream):
            if view is None or row.get("view", "main") == view:
                yield row


def stage_tables(run_root, expected_status_sha256, required_tables=None):
    """Require completed child independent audits and exact collected CSVs."""
    run_root = Path(run_root).resolve()
    status_path = run_root / "results/analysis/status.json"
    verified_file(dict(path=str(status_path), sha256=expected_status_sha256), run_root)
    status = json.loads(status_path.read_text())
    arms_path = run_root / "results/arms.json"
    arms = json.loads(arms_path.read_text())
    if not status.get("passed") or not status.get("complete") or not arms.get("verified") or arms["stage"] != status["stage"]:
        raise ValueError("Stage displays require complete verified stage exports")
    frozen_path = run_root / "freeze.json"
    frozen = json.loads(frozen_path.read_text())
    for path, digest in {**frozen["source_sha256"], **frozen.get("authority_sha256", {})}.items():
        verified_file(dict(path=path, sha256=digest), run_root)
    if P.sha(run_root / "inventory.json") != frozen["inventory_sha256"]:
        raise ValueError("The native stage inventory changed")
    sources = dict(stage_status=dict(path=str(status_path), sha256=P.sha(status_path)),
                   stage_freeze=dict(path=str(frozen_path), sha256=P.sha(frozen_path)),
                   verified_sigma_arms=dict(path=str(arms_path), sha256=P.sha(arms_path)))
    claims = {float(row["sigma"]):row for row in status["arms"]}
    if set(claims) != {float(row["sigma"]) for row in arms["arms"]}:
        raise ValueError("Collected and verified sigma arms differ")
    for arm in arms["arms"]:
        child = Path(arm["path"])
        manifest_path = child / "analysis_manifest.json"
        verified_file(dict(path=str(manifest_path), sha256=arm["analysis_manifest_sha256"]), child)
        claim = claims[float(arm["sigma"])]
        for filename, claim_key in (("status.json", "status_sha256"),
                                    ("independent_estimand_audit.json", "audit_sha256")):
            path = child / "results/analysis" / filename
            verified_file(dict(path=str(path), sha256=claim[claim_key]), child)
            payload = json.loads(path.read_text())
            if (not payload.get("passed") or not payload.get("complete") or
                    payload.get("condition_groups") != arm["condition_groups"] or
                    payload.get("analysis_manifest_sha256") != arm["analysis_manifest_sha256"]):
                raise ValueError("A sigma arm is missing complete independent estimand coverage")
            sources[f"sigma{arm['sigma']}_{filename}"] = dict(path=str(path), sha256=P.sha(path))
    tables = {}
    for name, digest in status["outputs"].items():
        path = run_root / "results/analysis" / name
        verified_file(dict(path=str(path), sha256=digest), run_root)
        sources[name] = dict(path=str(path), sha256=digest)
        if required_tables is not None and name not in required_tables:
            continue
        # Curves retain all lags of the main view. Window and cadence
        # comparisons use their separate metrics, not multi-million-row copies.
        view = "main" if name in ("heading_curves.csv", "site_curves.csv") else None
        tables[name] = P.display_rows(list(stream_rows(path, view=view)))
    return tables, sources, status


def select_observation_curves(metrics, by_stage):
    """Use the frozen observation choice; refuse response-dependent fallback."""
    selections = {}
    for row in metrics:
        key = condition_key(row)
        origin = row.get("common_observation_arm", row["selected_observation_arm"])
        stage = ("E4" if row["selected_observation_arm"] == "preselected_E4" else
                 "E3" if origin == "common_new_E3_size" else "E2")
        if key in selections:
            raise ValueError("An observation arm contains duplicate condition/view identities")
        selections[key] = (stage, P.num(row, "native_horizon"))
    result = []
    seen = set()
    for stage, rows in by_stage.items():
        for row in rows:
            key = condition_key(row)
            selected = selections.get(key)
            if selected is None or selected[0] != stage:
                continue
            sigma = P.num(row, "analysis_sigma", "sigma", default=1.)
            if sigma != 1.:
                continue
            if not math.isclose(P.num(row, "native_horizon"), selected[1]):
                raise ValueError("A selected curve belongs to a different native horizon")
            result.append(dict(row, analysis_arm=metrics[0]["analysis_arm"]))
            seen.add(key)
    if seen != set(selections):
        raise ValueError("A selected observation has no exact curve; common fallback is prohibited")
    return result


def presentation_bounds(fig):
    """Check visible tick labels, legends and document text against the page."""
    fig.canvas.draw(); renderer = fig.canvas.get_renderer()
    width, height = fig.canvas.get_width_height(); records = []
    artists = [("document_text", item) for item in fig.texts]
    artists.extend(("legend",item) for item in fig.legends)
    plot_boxes = []
    for index, ax in enumerate(fig.axes):
        if not ax.get_visible():
            continue
        plot_boxes.append((index, ax.get_window_extent(renderer)))
        artists.extend(("annotation", item) for item in ax.texts)
        if not ax.axison:
            continue
        # Locators can retain visible=True labels beyond the view interval;
        # Matplotlib does not draw those ticks. Check only rendered positions.
        for axis in (ax.xaxis, ax.yaxis):
            lower, upper = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if lower <= tick.get_loc() <= upper:
                    artists.extend(("tick", item) for item in (tick.label1, tick.label2))
        artists.extend(("axis_text", item) for item in (ax.xaxis.label, ax.yaxis.label, ax.title,
                                                       ax.xaxis.get_offset_text(), ax.yaxis.get_offset_text()))
        if ax.get_legend() is not None:
            artists.append(("legend", ax.get_legend()))
    for kind, artist in artists:
        if not artist.get_visible() or (hasattr(artist, "get_text") and not artist.get_text()):
            continue
        box = artist.get_window_extent(renderer)
        inside = bool(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= width and box.y1 <= height)
        records.append(dict(kind=kind, text=artist.get_text() if hasattr(artist, "get_text") else None,
                            document_role=artist.get_gid(),
                            bbox_pixels=[float(value) for value in box.extents], within_page=inside))
    collisions=[]
    def intersects(first,second):
        return min(first[2],second[2])>max(first[0],second[0]) and min(first[3],second[3])>max(first[1],second[1])
    for document in records:
        if document["document_role"] not in ("document_header","document_footer"):
            continue
        for other in records:
            if other is document:
                continue
            if other["kind"] in ("tick","axis_text","legend","annotation","document_text") and intersects(document["bbox_pixels"],other["bbox_pixels"]):
                collisions.append(dict(document_role=document["document_role"],overlap_with=other))
        for index,box in plot_boxes:
            if intersects(document["bbox_pixels"],list(box.extents)):
                collisions.append(dict(document_role=document["document_role"],overlap_with_plot_axis=index))
    return dict(passed=all(row["within_page"] for row in records) and not collisions,
                within_page=all(row["within_page"] for row in records),
                header_footer_clear_of_axes=not collisions,collisions=collisions,
                canvas_pixels=[int(width), int(height)], items=records)


class Book:
    def __init__(self, output, name, context, sources, section=None):
        if name.startswith("figure4_E0") or name == "figure4_S9_E0_protrusion_statistics":
            raise ValueError("Accepted E0 display names are immutable")
        self.output, self.name, self.context = Path(output).resolve(), name, context
        self.sources, self.section, self.pages = sources, section, []
        for folder in ("data", "vectors", "metadata"):
            (self.output / "supporting" / folder).mkdir(parents=True, exist_ok=True)
        self.path = self.output / f"supporting/vectors/{name}.pdf"
        if self.path.exists():
            raise ValueError("Use a new named display closure rather than replacing existing context exports")
        self.pdf = PdfPages(self.path)

    def save(self, fig, page, title, source_rows, note):
        if not source_rows:
            plt.close(fig)
            raise ValueError("Every displayed page needs complete source rows")
        compressed=len(source_rows)>50000
        csv_path = self.output / (f"supporting/data/{self.name}_{len(self.pages)+1:03d}_{page}.csv"+(".gz" if compressed else ""))
        fields = list(dict.fromkeys(key for row in source_rows for key in row))
        opener=gzip.open if compressed else open
        with opener(csv_path,"wt",newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(source_rows)
        section = self.section or page.split("_")[0]
        title_text=f"{section} · {self.context} · {title}"
        header=fig.suptitle(textwrap.fill(title_text,width=int(fig.get_size_inches()[0]*9)), x=.055, ha="left", y=.97, fontsize=12)
        header.set_gid("document_header")
        text = self.context + ". " + note + " × at the axis floor means unidentified, not zero. Original scientific fields and display annotations remain in CSV."
        if any(key.startswith('display_') and key.endswith('_CI_marker') for row in source_rows for key in row):
            text+=' CI∞ next to a finite point marks an upper-unbounded bootstrap CI; finite lower endpoints are retained. CI[0,∞) spans the full nonnegative range.'
        footer=fig.text(.055, .015, "\n".join(textwrap.wrap(text, 155)) + "\nCSV: " + csv_path.name, fontsize=6, va="bottom")
        footer.set_gid("document_footer")
        bounds = presentation_bounds(fig)
        if not bounds["passed"]:
            raise ValueError("Page text is outside the physical page or header/footer intersects axes: " +
                             str(dict(outside=[row for row in bounds["items"] if not row["within_page"]],collisions=bounds["collisions"])))
        self.pdf.savefig(fig); plt.close(fig)
        self.pages.append(dict(section=section, label=page, title=title, pdf_page=len(self.pages)+1,
                               csv=str(csv_path), csv_sha256=P.sha(csv_path), rows=len(source_rows),
                               presentation_bounds=bounds))

    def close(self):
        self.pdf.close()
        report = dict(passed=True, complete=True, scope="Scientific display closure only; full Figure 4 export/framing is separate",
            context=self.context, producer_state="stopped", generated_utc=datetime.now(timezone.utc).isoformat(),
            native_integrations_run=0, final_main_framing_selected=False,
            pdf=dict(path=str(self.path), sha256=P.sha(self.path), pages=len(self.pages)), pages=self.pages,
            source_files=self.sources, renderer_source=dict(path=str(Path(__file__).resolve()), sha256=P.sha(__file__)),
            shared_plot_source=dict(path=str(Path(P.__file__).resolve()), sha256=P.sha(P.__file__)),
            shared_receipt_source=dict(path=str(Path(Q.__file__).resolve()),sha256=P.sha(Q.__file__)))
        P.write_json(self.output / f"supporting/metadata/{self.name}_display_audit.json", report)
        return report


def sigma_metric_confidence(row,key):
    raw_key={"single_fraction":"single_front_fraction","multiple_fraction":"multiple_fronts_fraction"}.get(key,key)
    stem=raw_key.removesuffix('_fraction')
    return (P.num(row,key+'_ci_low',stem+'_ci_lower'),P.num(row,key+'_ci_high',stem+'_ci_upper'))


def sigma_pages(book, tables):
    rows = P.main_rows(tables["condition_metrics.csv"])
    controls = [row for row in tables["control_metrics.csv"] if row.get("view") == "main" and row.get("kind") == "mean"]
    sigmas = {P.num(row, "sigma") for row in rows}
    if sigmas != {.25, .5, 1., 2.}:
        raise ValueError("S7 requires all four observed sigma arms")
    for row in controls:
        if P.num(row, "bootstrap_resamples") != 0 or P.num(row, "n_seeds") != 0:
            raise ValueError("A deterministic matched control must not receive seed bootstrap uncertainty")
    metrics = (("tau_p", r"Empirical $\tau_p$", True), ("Dtheta", r"Short-lag $D_\theta$", False),
               ("multiple_fraction", "Multiple-region\ntime fraction", False),
               ("low_contrast_fraction", "Low-contrast\ntime fraction", False),
               ("single_fraction", "Single-region\ntime fraction", False))
    for key, label, logy in metrics:
        fig, axes = plt.subplots(4, 2, figsize=(11.7, 8.3), squeeze=False); plotted=[]
        for i, which_model in enumerate(P.MODELS):
            strengths = sorted({P.num(row, "f") for row in P.subset(rows, which_model=which_model)})
            if len(strengths) != 2: raise ValueError("Every model needs its two reviewed fixed means")
            for j, strength in enumerate(strengths):
                ax = axes[i, j];bound_entries=[]
                selected = P.subset(rows, which_model=which_model, strength=strength)
                reference = P.subset(controls, which_model=which_model, strength=strength)
                for size in sorted({P.num(row, "C_label") for row in selected}):
                    series = sorted(P.subset(selected, size=size), key=lambda row:P.num(row, "sigma"))
                    mean = P.subset(reference, size=size)
                    if len(mean) != 1: raise ValueError("S7 requires one horizon-matched deterministic mean")
                    qualifier=next(row for row in series if P.num(row,"sigma")==1.)
                    mean_row=dict(mean[0],display_sigma=0.,row_type="deterministic_mean",condition_qualifier_source_sigma=1.)
                    for field in ("forced_mean","forced_E1","forced_combined"):
                        if field in qualifier:mean_row[field]=qualifier[field]
                    combined = [mean_row] + [dict(row, display_sigma=row["sigma"], row_type="noisy") for row in series]
                    color = P.size_color(which_model, size)
                    if key != "tau_p":
                        ax.plot([P.num(row, "display_sigma") for row in combined], [P.num(row, key) for row in combined], color=color, lw=.8, label=f"C {size:g}")
                    for row in combined:
                        current = P.size_color(which_model, size, P.forced(annotation(row)) is True)
                        x = P.num(row, "display_sigma")
                        if key == "tau_p":
                            first=len(ax.lines)
                            state=P.interval_mark(ax, x, row, key, color=current, marker="s" if x==0 else "o", label=f"C {size:g}" if x==0 else None)
                            if state in ('bounded','unresolved'):bound_entries.append((row,state,ax.lines[first]))
                        else:
                            value=P.num(row,key)
                            if math.isfinite(value):
                                ax.plot(x,value,"s" if x==0 else "o",color=current,ms=3)
                                low,high=sigma_metric_confidence(row,key)
                                if x and math.isfinite(low) and math.isfinite(high):ax.vlines(x,low,high,color=current,lw=.7)
                            else:ax.plot(x,.025,"x",color=current,transform=ax.get_xaxis_transform(),ms=3)
                    plotted.extend(combined)
                if key=='tau_p':P.coincident_bound_markers(ax,bound_entries,key,x_field='display_sigma')
                ax.set_title(f"{P.LABELS[which_model]} · fixed mean f={strength:g}", fontsize=8)
                ax.set_xscale("symlog",linthresh=.25);ax.set_xticks([0,.25,.5,1,2],labels=["0",".25",".5","1","2"]);ax.set_xticks([],minor=True)
                P.polish(ax,xlabel="Lognormal fluctuation parameter σ" if i==3 else None,ylabel=label if j==0 else None,logy=logy)
                if logy:
                    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=4))
                    ax.yaxis.set_minor_formatter(NullFormatter())
                    ax.margins(y=.35)
                if i==0:ax.legend(frameon=False,fontsize=6)
        fig.subplots_adjust(left=.10,right=.97,top=.86,bottom=.23,hspace=.7,wspace=.3)
        book.save(fig,"S7_"+key,"Fixed-mean fluctuation sweep",plotted,
            "At σ=0, resolved deterministic means use squares; triangles retain bounds and × marks unidentified values. Nested triangles and ×n retain coincident sizes at unchanged coordinates. Chemical C20/40/60 correspond to Spring resting C₀1/2/3. Controls have no ensemble CI. Other arms use twelve common latent Gaussian seed blocks, verified by knot hashes, and whole-seed paired bootstrap. Grey uses combined mean/E1 qualification; σ is not mean input. Native time units differ; bounds remain bounds.")
    sigma_comparison_pages(book,tables)


def sigma_comparison_pages(book,tables):
    qualifiers={condition_key(row):row for row in P.main_rows(tables["condition_metrics.csv"]) if P.num(row,"sigma")==1.}
    rows=[]
    for row in P.main_rows(tables["sigma_comparisons.csv"]):
        qualifier=qualifiers[condition_key(row)]
        merged=dict(row,condition_qualifier_source_sigma=1.)
        for field in ("forced_mean","forced_E1","forced_combined"):
            if field in qualifier:merged[field]=qualifier[field]
        rows.append(annotation(merged))
    for key,label in (("tau_ratio","Persistence ratio\nto σ=1"),("Dtheta_ratio","Angular mobility ratio\nto σ=1"),
                      ("multiple_fronts_difference","Multiple-region fraction\nminus σ=1"),
                      ("low_contrast_difference","Low-contrast fraction\nminus σ=1")):
        fig,axes=plt.subplots(4,2,figsize=(11.7,8.3),squeeze=False);plotted=[]
        for i,which_model in enumerate(P.MODELS):
            strengths=sorted({P.num(row,"f") for row in P.subset(rows,which_model=which_model)})
            if len(strengths)!=2:raise ValueError("Paired sigma comparisons omit a reviewed fixed mean")
            for j,strength in enumerate(strengths):
                ax=axes[i,j];selected=P.subset(rows,which_model=which_model,strength=strength)
                for size in sorted({P.num(row,"C_label") for row in selected}):
                    series=sorted(P.subset(selected,size=size),key=lambda row:P.num(row,"sigma"))
                    color=P.size_color(which_model,size,P.forced(series[0]) is True)
                    for row in series:
                        x,value=P.num(row,"sigma"),P.num(row,key)
                        low,high=P.num(row,key+"_ci_lower"),P.num(row,key+"_ci_upper")
                        if math.isfinite(value):ax.plot(x,value,"o",color=color,ms=3)
                        else:ax.plot(x,.025,"x",color=color,ms=3,transform=ax.get_xaxis_transform())
                        if math.isfinite(low) and math.isfinite(high):
                            ax.vlines(x,low,high,color=color,lw=.8)
                        elif math.isfinite(low):ax.plot(x,low,"^",mfc="none",color=color,ms=4)
                    ax.plot([],[],color=color,label=f"C {size:g}")
                ax.axhline(1 if "ratio" in key else 0,color=".7",ls=":",lw=.7)
                ax.set_title(f"{P.LABELS[which_model]} · fixed mean f={strength:g}",fontsize=8)
                P.polish(ax,xlabel="Lognormal fluctuation parameter σ" if i==3 else None,ylabel=label if j==0 else None,logx=True)
                ax.set_xticks([.25,.5,1,2],labels=[".25",".5","1","2"])
                ax.set_xticks([],minor=True)
                if i==0:ax.legend(frameon=False,fontsize=6)
                plotted.extend(selected)
        fig.subplots_adjust(left=.10,right=.97,top=.86,bottom=.23,hspace=.7,wspace=.3)
        book.save(fig,"S7_paired_"+key,label,plotted,
            "Direct paired comparisons use twelve identical latent Gaussian seed blocks across σ. Grey marks retain combined-qualified forced conditions. Vertical ranges are bootstrap confidence intervals; open upward triangles mark finite lower confidence endpoints when upper endpoints are unbounded. × marks an unidentified point; no midpoint is invented. The σ=1 reference line is an identity, not an additional ensemble. Every unavailable value and valid-bootstrap fraction is retained in CSV.")
    flags=P.main_rows(tables["E1_forced_flags.csv"])
    fig,axes=plt.subplots(4,1,figsize=(11.7,8.3))
    for ax,which_model in zip(axes,P.MODELS):
        ax.axis("off");selected=P.subset(flags,which_model=which_model)
        ax.text(0,.95,P.LABELS[which_model],va="top",fontsize=11,weight="bold")
        lines=[]
        for row in sorted(selected,key=lambda r:(P.num(r,"f"),P.num(r,"C_label"))):
            lines.append(f"f={P.num(row,'f'):g}, C={P.num(row,'C_label'):g}: E1 qualification {P.yesno(row.get('forced_E1'))}; lower-σ abnormal regime persists {P.yesno(row.get('low_sigma_abnormal_persistence'))}; matched mean agrees {P.yesno(row.get('mean_agrees_with_sigma1_abnormal'))}.")
        ax.text(0,.72,"\n".join(lines),va="top",fontsize=8)
    fig.subplots_adjust(left=.07,right=.96,top=.86,bottom=.21,hspace=.6)
    book.save(fig,"S7_forced_flags","Observed low-σ and matched-mean qualification",flags,
        "The approved operational rule uses the same abnormal majority at σ=1,0.5,0.25 and the observed matched σ=0 mean. This is finite-control evidence, not an analytic σ→0 limit. The original matched-mean forced flag remains a separate field.")


def display_sigma_one(book, tables, *, include_legacy_site_sensitivities=False):
    selected = {}
    for name, rows in tables.items():
        if name.endswith(".csv"):
            selected[name] = [annotation(row) for row in rows if P.num(row, "analysis_sigma", "sigma", default=1.) == 1.]
    for function in (P.regime_pages, P.heading_pages, P.size_pages, P.ratio_page,
                     P.site_pages if include_legacy_site_sensitivities else P.physical_mobility_page, P.noise_tolerance_pages):
        function(book, selected)


def observation_tables(observation_root,digest,input_spec,*,include_legacy_site_sensitivities=False):
    observation_root,input_spec=Path(observation_root).resolve(),Path(input_spec).resolve()
    status_path=observation_root/"status.json"
    verified_file(dict(path=str(status_path),sha256=digest),observation_root)
    status=json.loads(status_path.read_text())
    if (not status.get("passed") or not status.get("complete") or status.get("analysis_arm") not in ("common","longest")
            or float(status.get("sigma",math.nan))!=1. or float(status.get("forcing_sigma",math.nan))!=1.
            or status["input_spec_sha256"]!=P.sha(input_spec)):
        raise ValueError("Observation displays require exact frozen common/longest sigma-one inputs")
    sources=dict(observation_status=dict(path=str(status_path),sha256=P.sha(status_path)),
                 selected_observation_inputs=dict(path=str(input_spec),sha256=P.sha(input_spec)))
    tables={}
    for name,digest in status["outputs"].items():
        path=observation_root/name;verified_file(dict(path=str(path),sha256=digest),observation_root)
        tables[name]=P.display_rows([annotation(row) for row in stream_rows(path)])
        sources[name]=dict(path=str(path),sha256=digest)
    required={"condition_metrics.csv","size_ratios.csv","noise_tolerance.csv"}
    if not required<=set(tables):raise ValueError("The complete frozen observation tables are required")
    metrics=tables["condition_metrics.csv"]
    if any(row["analysis_arm"]!=status["analysis_arm"] for row in metrics):
        raise ValueError("Different observation arms cannot be pooled")
    for records in status.get("source_stage_gates",{}).values():
        for receipt in records:verified_file(receipt,observation_root)
    spec=json.loads(input_spec.read_text())
    stages={}
    selected_stages={"E4" if row["selected_observation_arm"]=="preselected_E4" else
                     "E3" if row.get("common_observation_arm",row["selected_observation_arm"])=="common_new_E3_size" else "E2" for row in metrics}
    curve_names={"heading_curves.csv"}|({"site_curves.csv"} if include_legacy_site_sensitivities else set())
    for stage in sorted(selected_stages):
        declared=spec["stages"][stage]
        selected,proof,_=stage_tables(declared["root"],declared["status_sha256"],
                                      required_tables=curve_names)
        stages[stage]=selected;sources[stage]=proof
    primary=P.main_rows(metrics)
    references={condition_key(row):row for row in primary}
    for name in sorted(curve_names):
        selected=select_observation_curves(primary,{stage:table[name] for stage,table in stages.items()})
        for row in selected:
            reference=references[condition_key(row)]
            for flag in ("forced_mean","forced_E1","forced_combined"):
                if flag in reference:row[flag]=reference[flag]
        tables[name]=P.display_rows(selected)
    return tables,sources,status


def render_observation(observation_root,digest,input_spec,output,name):
    tables,sources,status=observation_tables(observation_root,digest,input_spec)
    context="Common observations" if status["analysis_arm"]=="common" else "Longest qualified observations"
    if status.get("checkpoint_only"):context+=" checkpoint"
    book=Book(output,name,context,sources)
    try:
        display_sigma_one(book,tables)
        return book.close()
    except Exception:
        book.pdf.close();raise


def render_stage(run_root, digest, output, name):
    tables, sources, status = stage_tables(run_root, digest,required_tables={"condition_metrics.csv","control_metrics.csv",
        "sigma_comparisons.csv","E1_forced_flags.csv","heading_curves.csv","site_curves.csv","size_ratios.csv","noise_tolerance.csv"})
    stage = status["stage"]
    book = Book(output, name, stage + (" fixed-mean σ sweep" if stage=="E1" else " completed native observations"), sources,
                section="S7" if stage=="E1" else "S8" if stage=="E2" else None)
    try:
        if stage=="E1": sigma_pages(book,tables)
        else: display_sigma_one(book,tables)
        return book.close()
    except Exception:
        book.pdf.close();raise


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest="command",required=True)
    stage=commands.add_parser("stage")
    stage.add_argument("--run-root",required=True,type=Path)
    observation=commands.add_parser("observation")
    observation.add_argument("--observation-root",required=True,type=Path)
    observation.add_argument("--input-spec",required=True,type=Path)
    for command in (stage,observation):
        command.add_argument("--status-sha256",required=True)
        command.add_argument("--output",required=True,type=Path)
        command.add_argument("--name",required=True)
    arguments=parser.parse_args()
    if arguments.command=="stage":
        report=render_stage(arguments.run_root,arguments.status_sha256,arguments.output,arguments.name)
    else:
        report=render_observation(arguments.observation_root,arguments.status_sha256,arguments.input_spec,arguments.output,arguments.name)
    print(json.dumps(report["pdf"],indent=2))
