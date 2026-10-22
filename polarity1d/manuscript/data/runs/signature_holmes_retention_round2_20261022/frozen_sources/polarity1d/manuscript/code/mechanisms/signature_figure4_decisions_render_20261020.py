"""Display the audited final Figure 4 criteria and observation qualifications.

This renderer does not choose a main framing, recompute criteria, fit models or
integrate native trajectories. Raw criterion truth and window/horizon sensitivity
are kept on separate pages, with every displayed record in its page CSV.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

import signature_figure4_interim_render_20261020 as P
import signature_figure4_package_20261020 as Q
import signature_figure4_stage_render_20261020 as S

VIEWS = ("main", "analysis_first_half", "analysis_second_half", "late_quarter",
         "cadence_x2", "cadence_x4", "cadence_x8")
VIEW_LABELS = ("Main", "First half", "Second half", "Late quarter", "Cadence ×2", "Cadence ×4", "Cadence ×8")
ARMS = ("common", "longest")
STATUS = {"pass": 2, "inconclusive": 1, "fail": 0}
STATUS_TEXT = {"pass": "criterion met", "inconclusive": "inconclusive", "fail": "criterion not met"}
COLORMAP = ListedColormap(("#efc1b8", "#ededed", "#b8d8c0"))
REQUIRED = ("criteria.csv", "strong_input_candidates.csv", "condition_metrics.csv", "size_ratios.csv",
            "noise_tolerance.csv", "window_qualification.csv", "horizon_qualification.csv",
            "original_seed_sensitivity.csv", "E3_reference_checks.csv", "fingerprint_effects.csv", "E1_condition_metrics.csv")


def nested(value):
    if isinstance(value, str):
        try:
            return json.loads(value)
        except (ValueError, TypeError):
            return value
    return value


def words(value):
    """Plain display text; absence never becomes zero or an infinite endpoint."""
    if value is None or value == "":
        return "unavailable"
    if value is True or value == "True":
        return "yes"
    if value is False or value == "False":
        return "no"
    if isinstance(value, float):
        if math.isinf(value):
            return "∞" if value > 0 else "−∞"
        if math.isnan(value):
            return "unavailable"
        return f"{value:.4g}"
    if isinstance(value, list):
        return ", ".join(words(item) for item in value) if value else "none observed"
    if isinstance(value, dict):
        return "; ".join(f"{words(key)}: {words(item)}" for key, item in value.items())
    value = str(value)
    if value in STATUS_TEXT:
        return STATUS_TEXT[value]
    return value.replace("_", " ").replace("None", "unavailable")


def load_final(bundle, producer_sha256, audit_sha256):
    """Bind every plotted table to the active frozen producer and its audit."""
    bundle = Path(bundle).resolve()
    directory = bundle / "results/final_analysis"
    producer_path = directory / "complete.json"
    audit_path = directory / "independent_final_decision_audit.json"
    Q.verified_file(dict(path=str(producer_path), sha256=producer_sha256), bundle)
    Q.verified_file(dict(path=str(audit_path), sha256=audit_sha256), bundle)
    producer, audit = (json.loads(path.read_text()) for path in (producer_path, audit_path))
    if not all(payload.get("passed") is True and payload.get("complete") is True for payload in (producer, audit)):
        raise ValueError("Final displays require complete producer and independent decision audit")
    if audit.get("producer_complete_sha256") != producer_sha256:
        raise ValueError("The independent audit covers a different final producer")
    if audit.get("final_outputs_sha256") != producer["outputs"]:
        raise ValueError("The independent audit covers different final scientific tables")
    active_path = bundle / "ACTIVE_FINAL_ANALYSIS.json"
    active = json.loads(active_path.read_text())
    freeze_path = Path(active["active_source_freeze"])
    Q.verified_file(dict(path=str(freeze_path), sha256=active["active_source_freeze_sha256"]), bundle)
    if producer.get("source_freeze_sha256") != active["active_source_freeze_sha256"]:
        raise ValueError("The final producer does not use the active reviewed source closure")
    freeze = json.loads(freeze_path.read_text())
    if producer["source_sha256"] != freeze["source_sha256"]:
        raise ValueError("The frozen final producer source inventory differs")
    for path, digest in freeze["source_sha256"].items():
        Q.verified_file(dict(path=path, sha256=digest), bundle)
    spec_path = Path(producer["input_spec_path"])
    Q.verified_file(dict(path=str(spec_path), sha256=producer["input_spec_sha256"]), bundle)
    sources = dict(final_producer=dict(path=str(producer_path), sha256=producer_sha256),
                   independent_final_decision_audit=dict(path=str(audit_path), sha256=audit_sha256),
                   active_final_analysis=dict(path=str(active_path), sha256=P.sha(active_path)),
                   final_source_freeze=dict(path=str(freeze_path), sha256=P.sha(freeze_path)),
                   exact_final_inputs=dict(path=str(spec_path), sha256=P.sha(spec_path)))
    tables = {}
    for name, digest in producer["outputs"].items():
        path = directory / name
        Q.verified_file(dict(path=str(path), sha256=digest), bundle)
        sources[name] = dict(path=str(path), sha256=digest)
        if name.endswith(".csv"):
            tables[name] = P.read_csv(path)
    if not set(REQUIRED) <= set(tables):
        raise ValueError("Final scientific table coverage is incomplete")
    validate_criteria(tables["criteria.csv"])
    return tables, sources


def criterion_order():
    return [(model, criterion) for model in P.MODELS for criterion in ("H1_model", "H2_model")] + [
        ("all", criterion) for criterion in ("H1_overall", "H2_overall", "H3_overall")]


def validate_criteria(rows):
    expected = {(arm, view, inclusive, model, criterion) for arm in ARMS for view in VIEWS
                for inclusive in (False, True) for model, criterion in criterion_order()}
    actual = [(row.get("arm"), row.get("view"), P.boolean(row.get("forced_inclusive")),
               row.get("model"), row.get("criterion")) for row in rows]
    if len(actual) != len(set(actual)) or set(actual) != expected:
        raise ValueError("All 308 raw criteria must cover both arms, forcing views and seven windows exactly")
    if any(row.get("status") not in STATUS for row in rows):
        raise ValueError("An unknown criterion status cannot become inconclusive silently")


def status_grid(ax, cells, xlabels, ylabels, *, title):
    values = np.full((len(ylabels), len(xlabels)), np.nan)
    for (i, j), value in cells.items():
        if value not in STATUS:
            raise ValueError("A declared qualification must have an explicit three-way status")
        values[i, j] = STATUS[value]
    ax.imshow(np.ma.masked_invalid(values), cmap=COLORMAP, vmin=0, vmax=2, aspect="auto", interpolation="none")
    for (i, j), value in cells.items():
        ax.text(j, i, {"pass": "met", "fail": "not met", "inconclusive": "?"}[value],
                ha="center", va="center", fontsize=6.5)
    ax.set_xticks(range(len(xlabels)), xlabels, rotation=40, ha="right", fontsize=7)
    ax.set_yticks(range(len(ylabels)), ylabels, fontsize=7)
    ax.set_title(title, loc="left", fontsize=9)
    ax.set_xticks(np.arange(-.5, len(xlabels)), minor=True)
    ax.set_yticks(np.arange(-.5, len(ylabels)), minor=True)
    ax.grid(which="minor", color="white", linewidth=.8)
    ax.tick_params(which="minor", length=0)


def forcing_label(inclusive):
    return "all conditions including operationally forced" if inclusive else "operationally unforced conditions"


def summary_pages(book, tables):
    rows = tables["criteria.csv"]
    order = criterion_order()
    for arm in ARMS:
        for inclusive in (False, True):
            selected = [row for row in rows if row["arm"] == arm and P.boolean(row["forced_inclusive"]) is inclusive]
            lookup = {(row["model"], row["criterion"], row["view"]): row for row in selected}
            cells = {(i, j): lookup[(model, criterion, view)]["status"]
                     for i, (model, criterion) in enumerate(order) for j, view in enumerate(VIEWS)}
            fig, ax = plt.subplots(figsize=(11.7, 8.3))
            status_grid(ax, cells, VIEW_LABELS,
                        [(P.LABELS.get(model, "All models") + " · " + criterion.replace("_model", "").replace("_overall", " overall"))
                         for model, criterion in order], title=arm.title() + " · " + forcing_label(inclusive))
            fig.subplots_adjust(left=.24, right=.96, bottom=.28, top=.86)
            fig.text(.24, .12, "Green: criterion met. Red: declared criterion not met. ?: inconclusive.\n"
                     "Raw criterion truth is shown here. A7 window, horizon and seed-cohort qualifications are separate S12 evidence.\n"
                     "An unmet criterion does not establish its biological opposite; no main framing is selected by this renderer.", fontsize=8)
            book.save(fig, f"S0_{arm}_{'inclusive' if inclusive else 'primary'}", "Declared hypothesis criteria across every view", selected,
                      "Both forcing views are retained. All seven declared windows/cadences are displayed; missing work is refused.")
            main = [lookup[(model, criterion, "main")] for model, criterion in order]
            detail_pages(book, main, f"S0_{arm}_{'inclusive' if inclusive else 'primary'}_details",
                         "Main-view criterion components · " + arm + " · " + forcing_label(inclusive),
                         "Criterion components are copied from the audited final tables, without recomputing truth.", rows_per_page=6)


def detail_lines(row):
    header = " · ".join(filter(None, (P.LABELS.get(row.get("model"), row.get("model", "")),
        words(row.get("criterion")), words(row.get("status")))))
    if row.get("f") not in (None, ""):
        header += " · f " + str(row["f"])
    lines = [header]
    for key in ("reason", "components", "reversal_inputs", "peak_drop_ranges", "evaluable_inputs", "passes", "inconclusive", "required", "total"):
        if key in row and row[key] not in (None, ""):
            lines.extend(textwrap.wrap(words(key) + ": " + words(nested(row[key])), 150))
    if row.get("criterion") == "H3_overall":
        pairs=nested(row.get("pairs",[]))
        lines.append(f"All {len(pairs)} preset distinct-model/mode pairs are displayed separately on H3 pair pages, including every inconclusive or unmet component.")
    return lines


def encoded(value):
    return json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(",",":"))


def three_way(value):
    value=nested(value)
    if not isinstance(value,dict) or value.get("status") not in STATUS:
        raise ValueError("An H3 component lacks its existing audited three-way status")
    return STATUS_TEXT[value["status"]]


def bounds_text(row,stem="tau"):
    if row is None:return "unavailable"
    lo,hi=P.num(row,stem+"_lower"),P.num(row,stem+"_upper")
    if P.boolean(row.get(stem+"_upper_unbounded")) is True:hi=math.inf
    return "["+words(lo)+", "+words(hi)+"]"


def fraction_text(row,stem):
    if row is None:return "unavailable"
    ci_stem=stem[:-len("_fraction")]+"_ci" if stem.endswith("_fraction") else stem+"_ci"
    return words(P.num(row,stem))+" "+bounds_text(row,ci_stem)


def flag_triplet(row):
    if row is None:return "?/?/?"
    return "/".join({True:"Y",False:"N",None:"?"}[P.boolean(row.get(name))]
                    for name in ("forced_mean","forced_E1","forced_combined"))


def effect_text(row):
    if row is None:return "unavailable"
    text=words(P.num(row,"estimate"))+" "+bounds_text(row,"ci")
    return text+"\nN "+words(P.num(row,"independent_small_seeds"))+"/"+words(P.num(row,"independent_large_seeds"))+"; valid "+words(P.num(row,"valid_bootstrap_fraction"))


def native_horizon_text(row):
    value=P.num(row or {},"native_horizon")
    return f"{value:g}" if math.isfinite(value) else "unavailable"


def sigma_text(check,source):
    mode=(words(source.get("regime"))+"; H "+native_horizon_text(source)) if source else "unavailable"
    native=("multi "+words(P.num(source,"multiple_fronts_fraction"))+", low "+words(P.num(source,"low_contrast_fraction"))+
            "; forced mean "+P.yesno(source.get("forced_mean",source.get("forced")))) if source else "source unavailable"
    if check is None:
        return "not declared" if source is None else "no eligible check\n"+mode+"\n"+native
    if P.boolean(check.get("present")) is not True:return "not present"
    return (mode+"\nqualified "+P.yesno(check.get("qualified"))+
            "; same "+P.yesno(check.get("same_abnormal_majority"))+"\n"+
            native)


def h3_projection(tables):
    """Join existing decisions/measurements without constructing an H3 gate."""
    metrics=tables["condition_metrics.csv"]
    metric_map={(row["arm"],P.model(row),P.num(row,"C_label"),P.num(row,"f"),row["view"]):row for row in metrics}
    if len(metric_map)!=len(metrics):raise ValueError("Duplicate H3 native observation identity")
    sizes={(arm,model):sorted({P.num(row,"C_label") for row in metrics if row["arm"]==arm and P.model(row)==model})
           for arm in ARMS for model in P.MODELS}
    if any(len(values)<2 for values in sizes.values()):raise ValueError("Every H3 arm/model needs its declared size endpoints")
    candidates=tables["strong_input_candidates.csv"]
    expected={(arm,view,inclusive,model,P.num(row,"f"))
              for arm in ARMS for model in P.MODELS for view in VIEWS for inclusive in (False,True)
              for row in metrics if row["arm"]==arm and P.model(row)==model and row["view"]==view
              and P.num(row,"C_label")==sizes[(arm,model)][-1] and P.num(row,"f")>.0001}
    actual=[(row["arm"],row["view"],P.boolean(row["forced_inclusive"]),P.model(row),P.num(row,"f")) for row in candidates]
    if len(actual)!=len(set(actual)) or set(actual)!=expected:
        raise ValueError("Every predeclared H3 model/input/view/forcing candidate must be displayed")
    effects={(row["arm"],P.model(row),P.num(row,"f"),row["view"],row["statistic"]):row for row in tables["fingerprint_effects.csv"]}
    if len(effects)!=len(tables["fingerprint_effects.csv"]):raise ValueError("Duplicate H3 fingerprint quotient identity")
    e1={(P.model(row),P.num(row,"C_label"),P.num(row,"f"),row["view"],P.num(row,"sigma")):row for row in tables["E1_condition_metrics.csv"]}
    if len(e1)!=len(tables["E1_condition_metrics.csv"]):raise ValueError("Duplicate preset H3 E1 observation")
    joined=[]
    for original in candidates:
        arm,view,model,f=original["arm"],original["view"],P.model(original),P.num(original,"f")
        small,large=sizes[(arm,model)][0],sizes[(arm,model)][-1]
        strong=metric_map[(arm,model,large,f,view)];weak=metric_map.get((arm,model,large,.0001,view))
        component,corroboration=nested(original["strong"]),nested(original["E1"])
        three_way(component);three_way(corroboration)
        checks=corroboration.get("sigma_checks",[])
        if len({P.num(check,"sigma") for check in checks})!=len(checks) or any(P.num(check,"sigma") not in (.25,.5,2.) for check in checks):
            raise ValueError("H3 preset sigma checks were duplicated or redefined")
        raw={name:effects.get((arm,model,f,view,name)) for name in ("region_count","normalized_contrast")}
        record=dict(original)
        record.update(display_large_strong_condition=encoded(strong),display_large_weak_condition=encoded(weak),
                      display_small_strong_condition=encoded(metric_map.get((arm,model,small,f,view))),
                      display_small_weak_condition=encoded(metric_map.get((arm,model,small,.0001,view))),
                      display_raw_fingerprint_effects=encoded(raw),
                      display_preset_E1_conditions=encoded({str(sigma):e1.get((model,large,f,view,sigma)) for sigma in (.25,.5,1.,2.)}))
        joined.append(record)
    return joined,effects


def scientific_table_pages(book,records,cells,headers,widths,label,title,note,*,rows_per_page=12):
    if len(records)!=len(cells):raise ValueError("H3 displayed records and table rows differ")
    def wrapped(value,column):
        limit=max(10,int(11.7*72*.925*widths[column]/(6.2*.65)))
        return "\n".join(textwrap.fill(line,limit) for line in str(value).split("\n"))
    all_values=[[wrapped(value,j) for j,value in enumerate(row)] for row in cells]
    display_headers=[wrapped(value,j) for j,value in enumerate(headers)]
    max_lines=max((value.count("\n")+1 for row in [display_headers]+all_values for value in row),default=1)
    # This only chooses pagination; the exact renderer below still refuses
    # any actual cell overflow. No evidence is truncated or omitted.
    capacity=max(1,int(8.3*72*.62/(max_lines*6.2*1.6+1))-1)
    rows_per_page=min(rows_per_page,capacity)
    for start in range(0,len(records),rows_per_page):
        source=records[start:start+rows_per_page];values=all_values[start:start+rows_per_page]
        fig,ax=plt.subplots(figsize=(11.7,8.3));ax.axis("off")
        fraction=(len(values)+1)/(rows_per_page+1)
        table=ax.table(cellText=values,colLabels=display_headers,colWidths=widths,cellLoc="left",bbox=[0,1-fraction,1,fraction])
        table.auto_set_font_size(False);table.set_fontsize(6.2)
        for (i,j),cell in table.get_celld().items():
            cell.PAD=.035
            if i==0:cell.set_facecolor("#edf1f4");cell.set_text_props(weight="bold")
        fig.subplots_adjust(left=.045,right=.97,top=.86,bottom=.24)
        fig.canvas.draw();renderer=fig.canvas.get_renderer();cell_bounds=[]
        for identity,cell in table.get_celld().items():
            outer=cell.get_window_extent(renderer);inner=cell.get_text().get_window_extent(renderer)
            inside=bool(inner.x0>=outer.x0 and inner.x1<=outer.x1 and inner.y0>=outer.y0 and inner.y1<=outer.y1)
            cell_bounds.append(dict(cell=list(identity),text=cell.get_text().get_text(),inside_cell=inside,
                                    cell_bbox_pixels=[float(value) for value in outer.extents],
                                    text_bbox_pixels=[float(value) for value in inner.extents]))
        if not all(item["inside_cell"] for item in cell_bounds):
            raise ValueError("An H3 table cell clips its evidence: "+str([row for row in cell_bounds if not row["inside_cell"]]))
        book.save(fig,f"{label}_{start//rows_per_page+1:03d}",title,source,note)
        book.pages[-1]["table_cell_bounds"]=cell_bounds


def h3_pages(book,tables):
    joined,effects=h3_projection(tables)
    common_note=("All values project existing audited truth. S/W=strong/weak. Native A4/E1 raw mode and fractions are separate from A5 optical fingerprints. Strong/weak τ intervals are observation bounds; weak input is the actual f=0.0001. The existing strong-input test uses strong upper < 0.5×weak lower, subject to source eligibility/forcing. M/E1/C flags mean matched-mean/preset-E1/combined. The 95% bound check is the existing confidence-bound qualifier. Y=yes, N=no, ?=unavailable. Raw quotient CIs do not override excluded candidate fingerprints; upper∞ follows the source unbounded flag. Whole-seed sizes/valid-bootstrap fractions remain explicit. No component or H3 gate is redefined.")
    for arm in ARMS:
        for view in VIEWS:
            for inclusive in (False,True):
                records=sorted([row for row in joined if row["arm"]==arm and row["view"]==view and P.boolean(row["forced_inclusive"]) is inclusive],
                               key=lambda row:(P.MODELS.index(P.model(row)),P.num(row,"f")))
                suffix=f"{arm}_{view}_{'inclusive' if inclusive else 'primary'}"
                title=arm+" · "+view.replace("_"," ")+" · "+forcing_label(inclusive)
                basic=[];fingerprints=[]
                for row in records:
                    strong=nested(row["display_large_strong_condition"]);weak=nested(row["display_large_weak_condition"])
                    e1=nested(row["E1"]);checks={P.num(check,"sigma"):check for check in e1.get("sigma_checks",[])}
                    e1_sources=nested(row["display_preset_E1_conditions"]);raw=nested(row["display_raw_fingerprint_effects"])
                    identity=P.LABELS[P.model(row)]+"\nf="+words(P.num(row,"f"))
                    basic.append([identity,words(row["mode"]),three_way(row["strong"]),three_way(row["E1"]),
                        bounds_text(strong)+"\n"+words(strong.get("tau_status")),bounds_text(weak)+"\n"+words(weak.get("tau_status") if weak else None),
                        "S "+words(P.num(strong,"tau_coverage"))+" / W "+words(P.num(weak or {},"tau_coverage"))+"\nH "+native_horizon_text(strong)+" / "+native_horizon_text(weak)+"\n95% bound check: "+{True:"Y",False:"N",None:"?"}[P.boolean(nested(row["strong"]).get("confidence_bound_qualifier"))],
                        "S multi "+fraction_text(strong,"multiple_fronts_fraction")+"\nS low "+fraction_text(strong,"low_contrast_fraction")+"\nW multi "+words(P.num(weak or {},"multiple_fronts_fraction"))+"; low "+words(P.num(weak or {},"low_contrast_fraction")),
                        "S "+flag_triplet(strong)+"\nW "+flag_triplet(weak)])
                    qualified=nested(row["fingerprints"])
                    fingerprints.append([identity]+[sigma_text(checks.get(sigma),e1_sources.get(str(sigma))) for sigma in (.25,.5,2.)]+
                        [effect_text(raw.get("region_count")),effect_text(raw.get("normalized_contrast")),
                         "count "+P.yesno(qualified.get("region_count",{}).get("available"))+"\ncontrast "+P.yesno(qualified.get("normalized_contrast",{}).get("available"))])
                scientific_table_pages(book,records,basic,["Model / input","Native mode","Strong-input\ncomponent","E1 component","Strong τ bounds\n(native units)","Weak τ bounds\n(native units)","Coverage S/W\nNative H S/W","Native largest-size fractions\nS=strong; W=weak","Forcing M/E1/C"],
                    [.09,.09,.09,.09,.12,.12,.10,.21,.09],"S0_H3_candidates_"+suffix,"H3 candidate components · "+title,common_note)
                scientific_table_pages(book,records,fingerprints,["Model / input","Preset σ=0.25","Preset σ=0.5","Preset σ=2","Raw optical count\nlarge/small quotient CI","Raw normalized contrast\nlarge/small quotient CI","Candidate fingerprint\navailable after qualification"],
                    [.10,.14,.14,.14,.18,.18,.12],"S0_H3_fingerprints_"+suffix,"H3 preset-σ checks and optical quotients · "+title,
                    common_note+" Empty E1 checks outside the preset two-input grid mean no declared corroboration at that input, not missing executed native work.",rows_per_page=10)
                overall=[row for row in tables["criteria.csv"] if row["arm"]==arm and row["view"]==view and P.boolean(row["forced_inclusive"]) is inclusive and row["criterion"]=="H3_overall"]
                if len(overall)!=1:raise ValueError("Exactly one audited H3 overall result is required per context")
                pair_records=[];pair_cells=[]
                for pair in nested(overall[0].get("pairs",[])):
                    three_way(pair)
                    components=pair.get("components",[]);separation=[item for item in components if item.get("criterion")=="H3_fingerprint_separation"]
                    if (len(components)!=5 or len(separation)!=1 or
                            [component.get("criterion") for component in components[:4]]!=["H3_strong_input","H3_strong_input","H3_E1_corroboration","H3_E1_corroboration"]):
                        raise ValueError("An H3 pair lost existing strong/E1/fingerprint components")
                    fingerprint_checks={item["criterion"]:item for item in separation[0]["components"]}
                    if set(fingerprint_checks)!={"region_count","normalized_contrast"}:raise ValueError("Both audited H3 fingerprint comparisons must be shown")
                    context=dict(arm=arm,view=view,forced_inclusive=inclusive)
                    record=dict(context,**{k:encoded(v) if isinstance(v,(dict,list)) else v for k,v in pair.items()},
                                display_pair_source=encoded(pair),display_overall_source=encoded(overall[0]))
                    pairs_effects=[]
                    for which in ("first","second"):
                        key=(arm,pair[which+"_model"],float(pair[which+"_f"]),view)
                        source={stat:effects.get(key+(stat,)) for stat in ("region_count","normalized_contrast")}
                        record["display_"+which+"_raw_effects"]=encoded(source);pairs_effects.append(source)
                    pair_records.append(record)
                    pair_cells.append([P.LABELS[pair['first_model']]+"\nf="+words(pair['first_f'])+"\n"+words(pair['first_mode']),
                        P.LABELS[pair['second_model']]+"\nf="+words(pair['second_f'])+"\n"+words(pair['second_mode']),three_way(pair),
                        "Strong A/B: "+three_way(components[0])+" / "+three_way(components[1])+"\nE1 A/B: "+three_way(components[2])+" / "+three_way(components[3]),
                        "A "+effect_text(pairs_effects[0]["region_count"])+"\nB "+effect_text(pairs_effects[1]["region_count"]),
                        "A "+effect_text(pairs_effects[0]["normalized_contrast"])+"\nB "+effect_text(pairs_effects[1]["normalized_contrast"]),
                        three_way(separation[0]),"count "+three_way(fingerprint_checks['region_count'])+"\ncontrast "+three_way(fingerprint_checks['normalized_contrast'])])
                if not pair_records:
                    pair_records=[dict(overall[0],display_planned_pairs=0)];pair_cells=[["No preset distinct-model / distinct-mode pairs"]+[""]*7]
                scientific_table_pages(book,pair_records,pair_cells,["First model/input\n(A)","Second model/input\n(B)","Audited pair\nstatus","Strong / E1\ncomponents","Optical count quotient\nA / B (each within model)","Normalized contrast quotient\nA / B (each within model)","Fingerprint pair\nseparation status","Separate quotient\ncomparison statuses"],
                    [.11,.11,.10,.16,.17,.17,.09,.09],"S0_H3_pairs_"+suffix,"Every preset H3 model pair · "+title,
                    common_note+" Pair status and every nested comparison are copied unchanged; only within-model dimensionless effects are compared across models, never raw activity/counts. Empty context pages state zero planned pairs explicitly.",rows_per_page=10)


def detail_pages(book, rows, label, title, note, rows_per_page=6):
    # Each source record stays whole in the page CSV. Long component text is
    # explicitly continued across pages, rather than clipped or omitted.
    pages, current, entries = [], [], []
    for row in rows:
        lines = detail_lines(row)
        while lines:
            room = 43 - len(current)
            if room < 4:
                pages.append((current, entries)); current, entries = [], []; room = 43
            take = min(len(lines), room - 1)
            current.extend(lines[:take] + [""]); lines = lines[take:]
            if row not in entries:
                entries.append(row)
            if lines:
                pages.append((current, entries)); current, entries = ["Continued criterion components"], [row]
    if current:
        pages.append((current, entries))
    for page, (lines, source) in enumerate(pages, 1):
        fig = plt.figure(figsize=(11.7, 8.3))
        fig.text(.06, .91, "\n".join(lines), va="top", fontsize=7.4, linespacing=1.3)
        book.save(fig, f"{label}_{page:02d}", title, source, note)


def qualification_pages(book, tables):
    windows = tables["window_qualification.csv"]
    component_names = ("A7_window_qualification", "analysis_first_half_versus_analysis_second_half",
                       "main_versus_late_quarter", "main_versus_cadence_x2", "main_versus_cadence_x4", "main_versus_cadence_x8")
    component_labels = ("Combined", "First/second half", "Main/late quarter", "Cadence ×2", "Cadence ×4", "Cadence ×8")
    for arm in ARMS:
        for model in P.MODELS:
            selected = [row for row in windows if row["arm"] == arm and row["model"] == model]
            if not selected:
                raise ValueError("A model/arm lacks required A7 qualifications")
            sizes, strengths = (sorted({P.num(row, key) for row in selected}) for key in ("C_label", "f"))
            lookup = {(P.num(row, "C_label"), P.num(row, "f")): row for row in selected}
            cells = {}
            for i, (size, name) in enumerate((size, name) for size in sizes for name in component_names):
                for j, strength in enumerate(strengths):
                    row = lookup.get((size, strength))
                    if row is None:
                        continue
                    components = {component["criterion"]: component["status"] for component in nested(row["components"])}
                    cells[i, j] = row["status"] if name == component_names[0] else components[name]
            labels = [f"C {size:g} · {text}" for size in sizes for text in component_labels]
            fig, ax = plt.subplots(figsize=(11.7, 8.3))
            status_grid(ax, cells, [f"{f:g}" for f in strengths], labels,
                        title=P.LABELS[model] + " · " + arm + " · finite-window/cadence agreement")
            ax.set_xlabel("Input strength f"); fig.subplots_adjust(left=.29, right=.97, top=.84, bottom=.22)
            book.save(fig, f"S12_windows_{arm}_{model}", "A7 qualification by condition and comparison", selected,
                      "These are interval-and-native-regime agreement checks. Angular-mobility disagreement remains separately recorded in components. Blank cells are outside the declared grid.")
    for name, context in (("horizon_qualification.csv", "Common versus selected longest observations"),
                          ("original_seed_sensitivity.csv", "Original retained Gaussian paths versus all 12 long paths")):
        rows = tables[name]
        for model in P.MODELS:
            selected = [row for row in rows if row["model"] == model]
            if not selected:
                raise ValueError("A model lacks required horizon/cohort sensitivity evidence")
            sizes, strengths = (sorted({P.num(row, key) for row in selected}) for key in ("C_label", "f"))
            lookup = {(P.num(row, "C_label"), P.num(row, "f"), row["view"]): row for row in selected}
            cells = {(i, j): lookup[(size, strength, view)]["status"]
                     for i, (size, view) in enumerate((size, view) for size in sizes for view in VIEWS)
                     for j, strength in enumerate(strengths) if (size, strength, view) in lookup}
            fig, ax = plt.subplots(figsize=(11.7, 8.3))
            status_grid(ax, cells, [f"{f:g}" for f in strengths],
                        [f"C {size:g} · {label}" for size in sizes for label in VIEW_LABELS], title=P.LABELS[model] + " · " + context)
            ax.set_xlabel("Input strength f"); fig.subplots_adjust(left=.28, right=.97, top=.84, bottom=.22)
            book.save(fig, f"S12_{name[:-4]}_{model}", context, selected,
                      "All seven views remain distinct. Exact original paths are a same-horizon subset, not additional independent replicates. A non-selected common observation is retained explicitly; no missing long observation is replaced.")
    references = tables["E3_reference_checks.csv"]
    if not references:
        raise ValueError("The retained E3 reference checks cannot be omitted")
    for start in range(0, len(references), 35):
        selected = references[start:start + 35]
        lines = ["Retained E3 observations remain references; they do not replace canonical common observations.", ""]
        for row in selected:
            identity = words(nested(row["key"]))
            lines.append(f"{identity}: present in common arm {P.yesno(row.get('matched_in_common'))}; {words(row.get('source'))}.")
        fig = plt.figure(figsize=(11.7, 8.3))
        fig.text(.06, .91, "\n".join(lines), va="top", fontsize=7.4, linespacing=1.3)
        book.save(fig, f"S12_E3_references_{start//35+1:02d}", "E3 retained-reference inventory", selected,
                  "All retained reference identities are shown, including unmatched designed inputs. A missing designed comparison is not a missing executed run; native execution and intake remain separately audited.")


def render(bundle, producer_digest, audit_digest, output, name):
    tables, sources = load_final(bundle, producer_digest, audit_digest)
    book = S.Book(output, name, "Completed audited Figure 4 criteria and qualifications", sources)
    # The standard Book receipt also binds this actual additional producer.
    book.sources["decision_display_source"] = dict(path=str(Path(__file__).resolve()), sha256=P.sha(__file__))
    try:
        summary_pages(book, tables)
        h3_pages(book,tables)
        qualification_pages(book, tables)
        return book.close()
    except Exception:
        book.pdf.close()
        raise


def reviewed_memo_rows(memo,tables,producer_digest,audit_digest):
    """Bind the root's interpretation to unchanged actual criterion records."""
    if (memo.get("approved_for_candidate_framing") is not True or type(memo.get('no_clean_mapping')) is not bool or
            memo.get("producer_complete_sha256")!=producer_digest or
            memo.get("independent_final_decision_audit_sha256")!=audit_digest or not memo.get("framing")):
        raise ValueError("The interpretation memo must bind the independently completed actual criteria")
    summaries=memo.get("model_summaries",[]);identities=[P.model(row) for row in summaries]
    if len(identities)!=4 or set(identities)!=set(P.MODELS):
        raise ValueError("S0 needs an explicit pattern and limitations for every model")
    alternatives=memo.get("alternative_framings",[])
    if memo.get("no_clean_mapping") is True and not 2<=len(alternatives)<=3:
        raise ValueError("No clean mapping requires two or three explicit alternative framings")
    lookup={(row['arm'],row['view'],P.boolean(row['forced_inclusive']),P.model(row),row['criterion']):row
            for row in tables['criteria.csv']}
    result=[]
    for row in summaries:
        if not all(row.get(key) for key in ('pattern','qualifications','numeric_shortfalls','criterion_refs')):
            raise ValueError("A model memo must expose its pattern, qualifications and deciding numerical shortfalls")
        references=[]
        for ref in row['criterion_refs']:
            key=(ref['arm'],ref['view'],P.boolean(ref['forced_inclusive']),P.model(ref),ref['criterion'])
            if key not in lookup or ref.get('status')!=lookup[key]['status']:
                raise ValueError("A memo criterion reference differs from actual audited truth")
            references.append(lookup[key])
        result.append(dict(section='model',model=P.model(row),heading=P.LABELS[P.model(row)],
            pattern=row['pattern'],qualifications=row['qualifications'],numeric_shortfalls=row['numeric_shortfalls'],
            actual_criterion_records=encoded(references)))
    for alternative in alternatives:
        if not all(alternative.get(key) for key in ('framing','supporting_panels','weaknesses','further_runs')):
            raise ValueError("Every alternative needs support, weaknesses and resolving runs")
        result.append(dict(section='alternative',heading=alternative['framing'],
            supporting_panels=alternative['supporting_panels'],weaknesses=alternative['weaknesses'],
            further_runs=alternative['further_runs']))
    if not memo.get('overall_interpretation') or not memo.get('retained_qualifications'):
        raise ValueError("S0 must state the reviewed overall interpretation and retained limits")
    result.insert(0,dict(section='overall',heading=memo['framing'],
        overall_interpretation=memo['overall_interpretation'],retained_qualifications=memo['retained_qualifications'],
        no_clean_mapping=memo.get('no_clean_mapping')))
    return result


def render_review_memo(bundle,producer_digest,audit_digest,memo_record,output,name):
    tables,sources=load_final(bundle,producer_digest,audit_digest)
    memo_path=Q.verified_file(memo_record,bundle);memo=json.loads(memo_path.read_text())
    rows=reviewed_memo_rows(memo,tables,producer_digest,audit_digest)
    sources['reviewed_interpretation_memo']=dict(path=str(memo_path),sha256=P.sha(memo_path))
    sources['memo_display_source']=dict(path=str(Path(__file__).resolve()),sha256=P.sha(__file__))
    book=S.Book(output,name,'Reviewed final Figure 4 interpretation',sources,section='S0')
    try:
        for sequence,row in enumerate(rows,1):
            paragraphs=[row['heading']]
            for key in ('overall_interpretation','retained_qualifications','pattern','numeric_shortfalls',
                        'qualifications','supporting_panels','weaknesses','further_runs'):
                if key in row:paragraphs.extend(['',words(key)+': '+words(row[key])])
            if row['section']=='model':
                for ref in nested(row['actual_criterion_records']):
                    paragraphs.extend(['','Audited deciding record:']+detail_lines(ref))
            lines=[]
            for paragraph in paragraphs:lines.extend(textwrap.wrap(paragraph,142) if paragraph else [''])
            for start in range(0,len(lines),30):
                fig=plt.figure(figsize=(11.7,8.3))
                fig.text(.06,.86,'\n'.join(lines[start:start+30]),va='top',fontsize=8.4,linespacing=1.3)
                book.save(fig,f'S0_memo_{sequence:02d}_{start//30+1:02d}',
                    'Interpretation and alternative framings · '+row['heading'],[row],
                    'Reviewed interpretation remains distinct from raw hypothesis truth and its window/horizon qualifications. '
                    'All deciding criterion records and complete memo text remain in the page CSV. '
                    'Operational angular mobility and geometric conversion are not resolved long-lag persistence or measured migration.')
        return book.close()
    except Exception:
        book.pdf.close();raise


def prefix_recovery_rows(summary,rows):
    """Check the engineering census; no biological outcome enters selection."""
    total=summary['selected_native_trajectories'];repairs=summary['actual_prefix_repairs']
    if (summary.get('passed') is not True or summary.get('numerical_integrity_only') is not True or
            summary.get('scientific_long_arm_qualification') is not False or
            summary.get('response_outcomes_used_for_selection') is not False or
            summary.get('all_failed_originals_retained') is not True or
            summary.get('no_passing_or_live_originals_rerun') is not True or
            summary.get('selected_native_failures')!=0 or summary.get('missing_native_trajectories')!=0):
        raise ValueError('E4 recovery display requires the accepted numerical-only, failure-only source selection')
    if (summary['original_attempted_native_trajectories']!=total or
            summary['original_native_integrations_reaching_full_horizon']!=total or
            summary['original_prefix_validation_passed']+summary['original_prefix_validation_failed']!=total or
            summary['original_prefix_validation_failed']!=repairs or
            summary['qualified_old_cap_restarts']!=repairs or
            summary['passing_originals_preferred']!=summary['original_prefix_validation_passed'] or
            summary['total_native_trajectories_including_failed_attempts']!=total+repairs):
        raise ValueError('Original full-horizon attempts and selected repair census disagree')
    identities=[int(row['original_slot']) for row in rows]
    counts=Counter(P.model(row) for row in rows)
    if (len(rows)!=repairs or len(identities)!=len(set(identities)) or dict(counts)!=summary['repair_model_counts'] or
            sum(int(row['prefix_checks']) for row in rows)!=summary['prefix_checks'] or
            sum(int(row['bitwise_equal_checks']) for row in rows)!=summary['bitwise_equal_prefix_checks'] or
            not math.isclose(max(float(row['maximum_absolute_difference']) for row in rows),
                             summary['maximum_absolute_prefix_difference'],rel_tol=1e-14,abs_tol=0)):
        raise ValueError('Every actual prefix repair and all original numeric checks must remain in the display source')
    rule=summary['unchanged_native_prefix_rule']
    if rule['relative_tolerance']!=2e-7 or rule['absolute_tolerance']!=2e-9:
        raise ValueError('Recovery does not authorize a changed native prefix tolerance')
    context={key:summary[key] for key in ('original_attempted_native_trajectories','original_native_integrations_reaching_full_horizon',
             'original_prefix_validation_passed','original_prefix_validation_failed','actual_prefix_repairs',
             'total_native_trajectories_including_failed_attempts','selected_native_trajectories','prefix_checks',
             'bitwise_equal_prefix_checks','maximum_absolute_prefix_difference')}
    return [dict(row,display_summary_counts=encoded(context),display_prefix_rule=encoded(rule),
                 display_scope='Numerical engineering repair only; not stationarity, persistence, or scalar-observation convergence')
            for row in rows]


def render_prefix_recovery(summary_record,output,name):
    summary_path=Q.verified_file(summary_record,Path.cwd());summary=json.loads(summary_path.read_text())
    rows=P.read_csv(Q.verified_file(summary['cases_csv'],summary_path.parent))
    rows=prefix_recovery_rows(summary,rows)
    sources=dict(prefix_recovery_summary=dict(path=str(summary_path),sha256=P.sha(summary_path)),
                 exact_prefix_repair_cases=summary['cases_csv'])
    for key in ('selected_native_manifest','selection_authority','recovery_recipe','original_attempt_census',
                'original_independent_attempt_audit','selected_native_master','selected_numerical_control_prefix_gate'):
        path=Q.verified_file(summary[key],summary_path.parent)
        sources[key]=dict(path=str(path),sha256=P.sha(path))
    native=json.loads(Path(sources['selected_native_master']['path']).read_text())
    if (native.get('passed') is not True or native.get('full_native_run_covered') is not True or
            native.get('selected_workers')!=summary['selected_native_trajectories'] or
            native.get('native_manifest_sha256')!=sources['selected_native_manifest']['sha256']):
        raise ValueError('Recovery display must bind the full independent selected native-source audit')
    sources['display_source']=dict(path=str(Path(__file__).resolve()),sha256=P.sha(__file__))
    book=S.Book(output,name,'Accepted E4 numerical prefix repair; scientific outcomes assessed separately',sources,section='S12')
    try:
        fig=plt.figure(figsize=(11.7,8.3));ax=fig.add_axes([.08,.36,.33,.46])
        counts=[summary['repair_model_counts'].get(model,0) for model in P.MODELS]
        ax.bar(range(4),counts,color='#638a9b')
        ax.set_xticks(range(4),[P.LABELS[model] for model in P.MODELS],rotation=20,ha='right',fontsize=9)
        ax.set_ylabel('Qualified exact-cap prefix repairs',fontsize=10)
        ax.set_ylim(0,max(counts)*1.2)
        for i,count in enumerate(counts):ax.text(i,count+.2,str(count),ha='center',fontsize=10)
        lines=[f"All {summary['original_attempted_native_trajectories']} original integrations reached the full horizon.",
               f"{summary['original_prefix_validation_passed']} passed the original prefix gate and were retained unchanged.",
               f"{summary['original_prefix_validation_failed']} failed the prefix gate; precisely these were replayed to the exact old cap and continued.",
               f"All {summary['selected_native_trajectories']} selected histories pass the independent native/source/forcing/prefix gates.",
               '',f"{summary['prefix_checks']} native numeric prefix checks; {summary['bitwise_equal_prefix_checks']} bitwise equal.",
               f"Maximum absolute prefix difference: {summary['maximum_absolute_prefix_difference']:.5g}.",
               'Unchanged gate: relative tolerance 2×10⁻⁷; absolute tolerance 2×10⁻⁹.',
               '', 'One original normalization and global-time forcing closure; native kinetics, mesh, fields, observation clocks and seeds are unchanged.',
               'Failed attempts remain archived. No biological response was used to select replacements.',
               '', 'This is engineering equivalence, not a stationarity or persistence result. The separate historical Spring scalar-cadence q90 sensitivity remains unresolved.']
        body=[]
        for line in lines:body.extend(textwrap.wrap(line,74) if line else [''])
        fig.text(.48,.84,'\n'.join(body),va='top',fontsize=9.1,linespacing=1.25)
        book.save(fig,'S12_E4_prefix_recovery','E4 exact-cap repair retains the declared native prefix rule',rows,
            'All 51 repair records, original failure reasons and exact source/proof identities remain in the page CSV. '
            'Counts are engineering integration attempts, not additional independent scientific seeds. '
            'Successful replay neither resolves the historical Spring observation-cadence limitation nor qualifies long-horizon biological claims.')
        return book.close()
    except Exception:
        book.pdf.close();raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--producer-sha256")
    parser.add_argument("--audit-sha256")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--review-memo",type=Path)
    parser.add_argument("--review-memo-sha256")
    parser.add_argument("--prefix-recovery-summary",type=Path)
    parser.add_argument("--prefix-recovery-summary-sha256")
    args = parser.parse_args()
    if (args.review_memo is None)!=(args.review_memo_sha256 is None):parser.error('Memo path and SHA256 must be supplied together')
    if (args.prefix_recovery_summary is None)!=(args.prefix_recovery_summary_sha256 is None):parser.error('Recovery summary path and SHA256 must be supplied together')
    if args.prefix_recovery_summary is not None:
        if args.review_memo is not None:parser.error('Recovery engineering display and final interpretation memo are separate exports')
        report=render_prefix_recovery(dict(path=str(args.prefix_recovery_summary.resolve()),sha256=args.prefix_recovery_summary_sha256),args.output,args.name)
    elif any(value is None for value in (args.bundle,args.producer_sha256,args.audit_sha256)):
        parser.error('Scientific decision/memo displays require the complete final producer and independent audit')
    elif args.review_memo is not None:
        report=render_review_memo(args.bundle,args.producer_sha256,args.audit_sha256,
            dict(path=str(args.review_memo.resolve()),sha256=args.review_memo_sha256),args.output,args.name)
    else:report = render(args.bundle, args.producer_sha256, args.audit_sha256, args.output, args.name)
    print(json.dumps(report["pdf"], indent=2))
