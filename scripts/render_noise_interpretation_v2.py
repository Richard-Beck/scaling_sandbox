#!/usr/bin/env python3
"""Render an interpretation supplement from the completed, audited noise campaign.

No integration, filtering of seeds, estimator change, or native source mutation.
Run with the existing scaling_sandbox Python and -B. Large source CSVs are read
one row at a time; compact source rows retain every lag and all three windows.
"""
from __future__ import annotations

import csv
import argparse
import gzip
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_noise_20261017"
DATA = SOURCE / "supporting/data"
OUT = ROOT / "manuscript_outline/assets/noise_v2"
MODELS = ("wave_pinning", "goryachev", "holmes_model3")
LABELS = {"wave_pinning": "Wave-pinning", "goryachev": "Goryachev", "holmes_model3": "Holmes 3"}
STRENGTHS = (.001, .01)
COLORS = { .001: "#2673b8", .01: "#c75032" }
WINDOWS = ("full", "first_half", "last_half")
WINLABELS = {"full": "Full 0–32,000", "first_half": "First 0–16,000", "last_half": "Last 16,000–32,000"}
MARKERS = {"full": "o", "first_half": "<", "last_half": ">"}
SOURCE_INFO = {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rel(path):
    return str(Path(path).resolve().relative_to(ROOT))


def read_small(name):
    path = DATA / name
    with path.open(newline="") as f:
        rows = [dict(row,source_csv_row=line) for line,row in enumerate(csv.DictReader(f),2)]
    SOURCE_INFO[name] = dict(path=rel(path), sha256=sha(path), rows=len(rows))
    return rows


def stream_export(name, predicate, target):
    """Keep exact source strings and source row numbers, without bulk CSV loads."""
    path = DATA / name
    selected = []
    count = 0
    with path.open(newline="") as f, gzip.open(OUT / target, "wt", newline="", compresslevel=6) as g:
        reader = csv.DictReader(f)
        writer = csv.DictWriter(g, fieldnames=reader.fieldnames + ["source_csv_row"])
        writer.writeheader()
        for line, row in enumerate(reader, 2):
            count += 1
            if predicate(row):
                row["source_csv_row"] = line
                writer.writerow(row)
                selected.append(row)
    SOURCE_INFO[name] = dict(path=rel(path), sha256=sha(path), source_rows=count,
                             exported_rows=len(selected), export=target)
    return selected


def key(r):
    return r["model"], float(r["circumference"]), float(r["fraction"])


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def json_clean(value):
    if isinstance(value, dict):
        return {k:json_clean(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_clean(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def f(r, field):
    return float(r[field])


def native_reconciliation(seed_rows, direction, sites, traces):
    """Independently recompute all plotted lag-8000 and final-time point values.

    Intervals remain the already audited 1000 whole-seed bootstrap intervals;
    their exact strings are copied with source hashes rather than recomputed
    under a different random stream. Native descriptors were independently
    checked at all saved profiles by the completed original analysis audit.
    """
    groups = defaultdict(list)
    raw_checks = []
    for row in seed_rows:
        path = Path(row["archive"])
        assert sha(path) == row["archive_sha256"]
        with np.load(path, allow_pickle=False) as z:
            t = z["time"].copy()
            p = z["P_real"] + 1j*z["P_imag"]
            magnitude = abs(p)
            gate = ((z["current_half_contrast_regions"] == 1)
                    & (z["peaks_prominence_ge_005"] <= 1)
                    & (z["contrast"] >= .05) & (magnitude >= .1))
            site = z["dominant_region_position"].copy()
            total = z["M0"].copy()
        assert len(t) == 6401 and np.array_equal(t, np.arange(0., 32001., 5.))
        np.testing.assert_allclose(magnitude[-1], float(row["final_R"]), rtol=1e-13, atol=1e-13)
        record = dict(row=row, t=t, p=p, magnitude=magnitude, gate=gate, site=site, total=total)
        groups[(*key(row), row["condition"])].append(record)
        raw_checks.append(dict(case_id=row["case_id"], condition=row["condition"],
                               archive=rel(path), sha256=row["archive_sha256"],
                               scalar_samples=len(t), final_R=float(magnitude[-1])))
    checks = []
    direct_by_key = {}
    site_by_key = {}
    for ident, group in groups.items():
        if ident[-1] != "noise":
            continue
        assert len(group) == 12
        C = ident[1]
        for window in WINDOWS:
            start, end = {"full": (0,32000), "first_half": (0,16000), "last_half": (16000,32000)}[window]
            unit_product = site_square = site_cosine = 0.
            dn = sn = available = 0
            for case in group:
                mask = (case["t"] >= start) & (case["t"] <= end)
                p, magnitude, eligible, position = (case[k][mask] for k in ("p", "magnitude", "gate", "site"))
                nlag = 1600
                good = (magnitude[nlag:] >= .1) & (magnitude[:-nlag] >= .1)
                product = (p[nlag:] / magnitude[nlag:]) * np.conj(p[:-nlag] / magnitude[:-nlag])
                unit_product += float(product.real[good].sum())
                dn += int(good.sum())
                pair = eligible[nlag:] & eligible[:-nlag]
                delta = (position[nlag:] - position[:-nlag] + C/2) % C - C/2
                site_square += float(np.square(delta[pair]).sum())
                site_cosine += float(np.cos(2*np.pi*delta[pair]/C).sum())
                sn += int(pair.sum())
                available += len(pair)
            wanted = (*ident[:3], window)
            direct_by_key[wanted] = dict(value=unit_product/dn if dn else math.nan,
                                          eligible=dn/available, pairs=dn, available=available)
            site_by_key[wanted] = dict(value=site_square/sn if sn else math.nan,
                                       cosine=site_cosine/sn if sn else math.nan,
                                       eligible=sn/available, pairs=sn, available=available)
    maxerr = defaultdict(float)
    for row in direction:
        if float(row["lag"]) != 8000:
            continue
        native = direct_by_key[(*key(row), row["window"])]
        for field, rawfield in (("direction_correlation", "value"), ("eligible_pair_fraction", "eligible")):
            err = abs(float(row[field])-native[rawfield])
            assert err < 5e-12, (row, field, err)
            maxerr[field] = max(maxerr[field], err)
        checks.append(dict(kind="direction", model=row["model"], circumference=float(row["circumference"]),
                           fraction=float(row["fraction"]), window=row["window"], **native))
    for row in sites:
        if float(row["lag"]) != 8000:
            continue
        native = site_by_key[(*key(row), row["window"])]
        want = native["value" if row["metric"] == "shortest_physical_displacement_squared" else "cosine"]
        np.testing.assert_allclose(float(row["value"]), want, rtol=2e-13, atol=2e-13, equal_nan=True)
        assert abs(float(row["eligible_pair_fraction"])-native["eligible"]) < 1e-14
        if math.isfinite(want):
            maxerr["site_"+row["metric"]] = max(maxerr["site_"+row["metric"]], abs(float(row["value"])-want))
        checks.append(dict(kind="site_"+row["metric"], model=row["model"], circumference=float(row["circumference"]),
                           fraction=float(row["fraction"]), window=row["window"], **native))
    for row in traces:
        if row["metric"] != "R":
            continue
        ident = (*key(row),row["condition"])
        value = np.mean([g["magnitude"][-1] for g in groups[ident]])
        np.testing.assert_allclose(float(row["mean"]), value, rtol=1e-13, atol=1e-13)
        maxerr["final_mean_R"] = max(maxerr["final_mean_R"], abs(float(row["mean"])-value))
    return dict(cases=len(seed_rows), noisy_cases=sum(r["condition"] == "noise" for r in seed_rows),
                native_case_checks=raw_checks, summary_checks=checks, maximum_errors=dict(maxerr)), groups, site_by_key


def common_axes(ax, ylabel=None, ylim=None):
    ax.set_xlim(17,63)
    ax.set_xticks([20,40,60])
    ax.grid(axis="y", alpha=.18, lw=.6)
    ax.spines[["top","right"]].set_visible(False)
    if ylabel:
        ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)


def error_curve(ax, rows, value, lower="ci_lower", upper="ci_upper", color="k", marker="o", offset=0, style="-", transform=None, label=None):
    rows = sorted(rows,key=lambda r:float(r["circumference"]))
    x = np.array([float(r["circumference"])+offset for r in rows])
    y,lo,hi = [np.array([float(r[k]) for r in rows]) for k in (value,lower,upper)]
    if transform:
        y,lo,hi = transform(y),transform(lo),transform(hi)
    ax.plot(x,y,linestyle=style,marker=marker,ms=4,color=color,lw=1.25,label=label)
    finite = np.isfinite(y)&np.isfinite(lo)&np.isfinite(hi)
    ax.errorbar(x[finite],y[finite],yerr=np.maximum(0,np.stack([y[finite]-lo[finite],hi[finite]-y[finite]])),
                fmt="none",color=color,lw=.85,capsize=2)


def select(rows, model, fraction, window="full", **extra):
    return [r for r in rows if r["model"]==model and float(r["fraction"])==fraction and
            (not r.get("window") or r["window"]==window) and all(r.get(k)==v for k,v in extra.items())]


def savefig(fig, stem):
    fig.savefig(OUT/(stem+".png"), dpi=200, facecolor="white")
    fig.savefig(OUT/(stem+".pdf"), facecolor="white")
    plt.close(fig)


def render_main(persistence,direction,sites,traces,site_counts):
    fig, axes = plt.subplots(5,3,figsize=(12.4,14.8))
    fig.subplots_adjust(left=.085,right=.97,top=.87,bottom=.12,hspace=.43,wspace=.28)
    fig.suptitle("Angular memory, physical site motion and observation coverage",fontsize=17,x=.5,y=.975)
    fig.text(.5,.947,"Completed 32,000-unit ring histories • 12 independent seeds per condition • both fixed input strengths",ha="center",fontsize=10)
    fig.legend(handles=[Line2D([0],[0],color=COLORS[s],marker="o",label=f"Mean fraction f = {s:g}") for s in STRENGTHS],
               loc="upper center",bbox_to_anchor=(.5,.933),ncol=2,frameon=False,fontsize=10)
    rownames = ["A  Restricted vector correlation area", "B  Conditional angular direction memory", "C  Resultant magnitude at final time", "D  Conditional physical site displacement", "E  Eligible endpoint pairs at lag 8,000"]
    ylabels=["Signed area to lag 8,000\n(native time units)","Direction cosine\nat lag 8,000","Final relative magnitude R","RMS shortest displacement\n(native length units)","Fraction of all origin pairs"]
    for col,model in enumerate(MODELS):
        axes[0,col].set_title(LABELS[model],fontsize=13,pad=19,fontweight="bold")
        for row in range(5):
            ax=axes[row,col]
            common_axes(ax,ylabels[row] if col==0 else None,
                        [(-200,8500),(-.15,1.08),(-.08,1.16),(-.25,8.5),(-.07,1.08)][row])
            ax.text(0,1.04,rownames[row].split("  ",1)[0]+"  "+(rownames[row].split("  ",1)[1] if col==0 else ""),transform=ax.transAxes,fontsize=9,fontweight="bold")
            if row==4:
                ax.set_xlabel("Circumference (native length units)")
        for strength in STRENGTHS:
            color=COLORS[strength]
            offset=-.7 if strength==.001 else .7
            error_curve(axes[0,col],select(persistence,model,strength),"restricted_integral",color=color,offset=offset)
            error_curve(axes[1,col],select(direction,model,strength),"direction_correlation",color=color,offset=offset)
            error_curve(axes[2,col],select(traces,model,strength,condition="noise",metric="R"),"mean",color=color,offset=offset)
            mean=select(traces,model,strength,condition="mean",metric="R")
            error_curve(axes[2,col],mean,"mean",color=color,marker="^",offset=offset,style=":")
            physical=select(sites,model,strength,metric="shortest_physical_displacement_squared")
            if model != "goryachev":
                error_curve(axes[3,col],physical,"value",color=color,offset=offset,transform=np.sqrt)
            direction_rows=select(direction,model,strength)
            error_curve(axes[4,col],direction_rows,"eligible_pair_fraction",lower="eligible_ci_lower",upper="eligible_ci_upper",color=color,offset=offset)
            physical=sorted(physical,key=lambda r:float(r["circumference"]))
            axes[4,col].plot([float(r["circumference"])+offset for r in physical], [float(r["eligible_pair_fraction"]) for r in physical],
                             color=color,ls="--",marker="s",ms=3.5,lw=1.25)
        axes[0,col].axhline(8000,color=".65",ls=":",lw=.8)
        axes[1,col].axhline(0,color=".65",lw=.6)
        axes[2,col].axhline(1,color=".65",lw=.6)
        if model=="goryachev":
            axes[0,col].annotate("C20 high mean changes\nthe forced-state regime",xy=(20.7,220),xytext=(29,2000),fontsize=8,color=COLORS[.01],arrowprops=dict(arrowstyle="->",color=COLORS[.01],lw=.8))
            axes[2,col].annotate("Mean-input polarity loss",xy=(20.7,0),xytext=(26,.20),fontsize=8,color=COLORS[.01],arrowprops=dict(arrowstyle="->",color=COLORS[.01],lw=.8))
            axes[3,col].text(.5,.85,"No population drift curve",ha="center",transform=axes[3,col].transAxes,fontweight="bold",fontsize=10)
            axes[3,col].text(.5,.45,"High f: 1 / 0 / 0 eligible pairs\n(out of 57,612 at each size)\nLow f: 8.43 / 1.00 / 0.23% coverage",ha="center",transform=axes[3,col].transAxes,fontsize=9,linespacing=1.6)
            axes[3,col].text(.5,.12,"Computed fragments remain in the data.\nUnavailability does not mean zero drift.",ha="center",transform=axes[3,col].transAxes,fontsize=8,color=".3")
    axes[2,2].text(.04,.09,"● noisy ensemble   ▲ uniform mean",transform=axes[2,2].transAxes,fontsize=8)
    axes[4,2].text(.04,.10,"— direction gate   - - single-site gate",transform=axes[4,2].transAxes,fontsize=8)
    fig.text(.085,.073,"Angular geometry: Δθ = 2π Δs / C. Equal physical displacement subtends a smaller angle on a larger boundary.",fontsize=10,fontweight="bold")
    fig.text(.085,.048,"Wave-pinning: high-input physical displacement is smaller at C60 than C20. Holmes: higher angular memory coexists with\nnonmonotone physical displacement and decreasing single-site coverage. Goryachev requires a forced-state qualification.",fontsize=10,linespacing=1.5)
    fig.text(.085,.015,"Bars: source 95% whole-seed bootstrap intervals. Signed finite-window areas are not fitted or extrapolated persistence times.\nSite endpoints must pass the complete local-profile and R ≥ 0.1 gate; continued site identity is unproven. All windows are retained in the companion sheet.",fontsize=8,color=".25",linespacing=1.4)
    savefig(fig,"noise_memory_geometry_v2")


def render_adaptation(persistence,direction,sites):
    fig,axes=plt.subplots(3,3,figsize=(12.4,10.5))
    fig.subplots_adjust(left=.085,right=.97,top=.84,bottom=.12,hspace=.44,wspace=.29)
    fig.suptitle("Full and half-window diagnostics retain ongoing adaptation",fontsize=16,y=.975)
    fig.text(.5,.94,"Same maximum lag 8,000 in every window; no selected burn-in or selected seeds",ha="center",fontsize=11)
    handles=[Line2D([0],[0],color=COLORS[s],lw=2,label=f"f = {s:g}") for s in STRENGTHS]
    handles += [Line2D([0],[0],color=".3",marker=MARKERS[w],ls="-" if w=="full" else "none",markerfacecolor=".3",label=WINLABELS[w]) for w in WINDOWS]
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.5,.916),ncol=3,frameon=False,fontsize=9)
    for row,model in enumerate(MODELS):
        for col in range(3):
            common_axes(axes[row,col],None,[(-200,8500),(-.15,1.08),(-.25,8.5)][col])
            axes[row,col].set_title(["Restricted vector area","Direction cosine, lag 8,000","Site RMS displacement, lag 8,000"][col],fontsize=10)
            if col==0: axes[row,col].set_ylabel(LABELS[model]+"\nNative time units")
            if col==1: axes[row,col].set_ylabel("Conditional direction cosine")
            if col==2: axes[row,col].set_ylabel("Conditional native length")
            if row==2: axes[row,col].set_xlabel("Circumference")
        for strength in STRENGTHS:
            for wi,window in enumerate(WINDOWS):
                offset=(-1.5 if strength==.001 else 1.5)+(wi-1)*.8
                for col,source,value,extra,transform in [(0,persistence,"restricted_integral",{},None),(1,direction,"direction_correlation",{},None),(2,sites,"value",{"metric":"shortest_physical_displacement_squared"},np.sqrt)]:
                    if model=="goryachev" and col==2:continue
                    error_curve(axes[row,col],select(source,model,strength,window,**extra),value,color=COLORS[strength],marker=MARKERS[window],offset=offset,style="-" if window=="full" else "none",transform=transform)
        if model=="goryachev":
            axes[row,2].text(.5,.55,"All windows remain poorly eligible.\nNo population curve is assigned.\nComplete values and counts are exported.",ha="center",transform=axes[row,2].transAxes,fontsize=10,linespacing=1.8)
    fig.text(.085,.067,"The major restricted-area ordering survives both halves. This agreement does not establish stationarity.",fontsize=10,fontweight="bold")
    fig.text(.085,.035,"Physical motion needs a narrower claim: high-input Wave-pinning last-half RMS is 3.19 / 3.49 / 1.93 at C20 / C40 / C60.\nHolmes C60 high-input direction cosine shifts from 0.798 in the first half to 0.670 in the last; large-ring site coverage remains conditional.",fontsize=9,linespacing=1.5)
    savefig(fig,"noise_adaptation_v2")


def native_examples(seed_rows,reps,groups):
    selected=[r for r in reps if r["role"]=="within_condition_descriptor_medoid" and float(r["fraction"])==.01 and float(r["circumference"]) in (20,60)]
    arrays={}; medoid_audit=[]
    fig,axes=plt.subplots(3,2,figsize=(12.4,10.4))
    fig.subplots_adjust(left=.07,right=.93,top=.88,bottom=.12,hspace=.43,wspace=.16)
    fig.suptitle("Matched high-input native activity histories",fontsize=17,y=.975)
    fig.text(.5,.939,"C20 and C60 within-condition medoids • f = 0.01 • individual histories, never ensemble-average fields",ha="center",fontsize=10)
    for mi,model in enumerate(MODELS):
        native=[]
        for C in (20.,60.):
            choice=next(r for r in selected if r["model"]==model and float(r["circumference"])==C)
            group=sorted([r for r in seed_rows if key(r)==(model,C,.01) and r["condition"]=="noise"],key=lambda r:r["case_id"])
            descriptors=np.array([[float(r[k]) for k in ("final_memory","minimum_R","time_mean_R","multi_region_sample_fraction")] for r in group])
            center=np.median(descriptors,axis=0)
            scale=np.subtract(*np.percentile(descriptors,(75,25),axis=0))
            scale=np.where(scale>1e-12,scale,1.)
            distance=np.sum(abs((descriptors-center)/scale),axis=1)
            winner=group[int(np.argmin(distance))]
            assert winner["case_id"]==choice["case_id"]
            np.testing.assert_allclose(distance.min(),float(choice["medoid_distance"]),rtol=1e-13)
            path=Path(choice["profile_archive"])
            assert sha(path)==choice["profile_archive_sha256"]
            with np.load(path) as z:
                primary=z["primary"].copy();x=z["x"].copy();t=z["time"].copy()
            # Shift the static plotting cut so the initial maximum sits at C/2.
            shift=len(x)//2-int(np.argmax(primary[0]))
            display=np.roll(primary,shift,axis=1)
            prefix=f"{model}_C{int(C)}"
            arrays[prefix+"_primary"]=primary
            arrays[prefix+"_x"]=x
            arrays[prefix+"_time"]=t
            arrays[prefix+"_display_roll"]=np.array(shift)
            trace=next(g for g in groups[(model,C,.01,"noise")] if g["row"]["case_id"]==choice["case_id"])
            for k in ("p","magnitude","total","site","gate","t"):
                arrays[prefix+"_"+k]=trace[k]
            native.append((choice,display,x,t,prefix))
            medoid_audit.append(dict(model=model,circumference=C,fraction=.01,case_id=choice["case_id"],seed=choice["seed"],
                                     medoid_distance=float(distance.min()),profile_archive=rel(path),profile_sha256=sha(path),
                                     trace_archive=rel(Path(choice["archive"])),trace_sha256=choice["archive_sha256"],
                                     compact_prefix=prefix,display_roll=shift,descriptor_median=center.tolist(),descriptor_iqr_scale=scale.tolist(),
                                     all_candidates=[dict(case_id=r["case_id"],distance=float(d),descriptors=v.tolist()) for r,d,v in zip(group,distance,descriptors)]))
        vmax=max(float(d[1].max()) for d in native)
        for ci,(choice,display,x,t,prefix) in enumerate(native):
            ax=axes[mi,ci]
            image=ax.imshow(display,origin="lower",aspect="auto",extent=(0,float(choice["circumference"]),0,32000),vmin=0,vmax=vmax,cmap="viridis",rasterized=True)
            ax.set_title(f"{LABELS[model]}   C{int(float(choice['circumference']))}   {choice['case_id']}",fontsize=10)
            ax.set_xlabel("Physical boundary arclength (shifted fixed cut)")
            ax.set_ylabel("Native time" if ci==0 else "")
            ax.set_yticks([0,16000,32000])
            trace=arrays[prefix+"_magnitude"]
            ax.text(.02,.97,f"Final R = {trace[-1]:.3f}",transform=ax.transAxes,va="top",color="white",fontsize=9,bbox=dict(facecolor="black",alpha=.4,edgecolor="none",pad=2))
        box=axes[mi,1].get_position()
        cax=fig.add_axes([.947,box.y0,.012,box.height])
        fig.colorbar(image,cax=cax,label="Native activity")
    fig.text(.07,.073,"C20/C60 share one zero-based activity scale within each model; physical widths differ. The fixed cut keeps the initial peak intact.",fontsize=9)
    fig.text(.07,.039,"Medoids minimize median/IQR-standardized distance using final initial-axis memory, minimum R, time-mean R and multisite fraction.\nThese examples illustrate native local organization; the complete seed population, magnitude gates and pair coverage determine quantitative claims.",fontsize=9,linespacing=1.5)
    savefig(fig,"noise_native_examples_v2")
    np.savez_compressed(OUT/"noise_native_examples_v2.npz",**arrays)
    (OUT/"noise_native_examples_v2.selection.json").write_text(json.dumps(medoid_audit,indent=2)+"\n")
    return medoid_audit


def audit_native():
    OUT.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"pdf.fonttype":42,"ps.fonttype":42,"axes.labelsize":9,"xtick.labelsize":9,"ytick.labelsize":9})
    seeds=read_small("vision_noise_seed_summary.csv")
    reps=read_small("vision_noise_representative_histories.csv")
    persistence=read_small("vision_noise_persistence.csv")
    assert len(seeds)==243 and sum(r["condition"]=="noise" for r in seeds)==216
    write_csv(OUT/"noise_seed_summary_v2.csv",[dict(r,archive=rel(Path(r["archive"]))) for r in seeds])
    write_csv(OUT/"noise_persistence_v2.csv",persistence)
    correlations=stream_export("vision_noise_lag_correlations.csv",lambda r:r["condition"]=="noise","noise_vector_lags_v2.csv.gz")
    directions=stream_export("vision_noise_direction_correlations.csv",lambda r:r["condition"]=="noise" and float(r["gate"])==.1,"noise_direction_lags_v2.csv.gz")
    sites=stream_export("vision_noise_site_displacement.csv",lambda r:r["condition"]=="noise","noise_site_lags_v2.csv.gz")
    traces=stream_export("vision_noise_ensemble_traces.csv",lambda r:float(r["time"])==32000 and r["metric"] in ("R","total_active_ratio","weak_fraction_gate0p1"),"noise_final_magnitude_controls_v2.csv.gz")
    # All source-lag values remain exported; only a compact terminal lag is plotted.
    integral_check=[]
    bylag=defaultdict(list)
    for r in correlations:bylag[(*key(r),r["window"])].append(r)
    for r in persistence:
        if r["condition"]!="noise":continue
        lagrows=sorted(bylag[(*key(r),r["window"])],key=lambda q:float(q["lag"]))
        lag=np.array([float(q["lag"]) for q in lagrows]);c=np.array([float(q["correlation"]) for q in lagrows])
        assert np.array_equal(lag,np.arange(0,8001,5))
        value=float(np.trapezoid(c,lag))
        np.testing.assert_allclose(value,float(r["restricted_integral"]),rtol=2e-14,atol=2e-10)
        integral_check.append(dict(model=r["model"],circumference=float(r["circumference"]),fraction=float(r["fraction"]),window=r["window"],integral=value,source_integral=float(r["restricted_integral"]),absolute_error=abs(value-float(r["restricted_integral"]))))
    del correlations,bylag
    direction_terminal=[r for r in directions if float(r["lag"])==8000]
    site_terminal=[r for r in sites if float(r["lag"])==8000]
    del directions,sites
    reconciliation,groups,site_counts=native_reconciliation(seeds,direction_terminal,site_terminal,traces)
    reconciliation["restricted_integral_reintegration"]=integral_check
    # Exact source summary rows with source-file identity are the plotting input.
    rows=[]
    for name,collection in [("vision_noise_persistence.csv",[r for r in persistence if r["condition"]=="noise"]),("vision_noise_direction_correlations.csv",direction_terminal),("vision_noise_site_displacement.csv",site_terminal),("vision_noise_ensemble_traces.csv",traces)]:
        for row in collection:rows.append(dict(source_file=name,source_sha256=SOURCE_INFO[name]["sha256"],**row))
    write_csv(OUT/"noise_plotted_source_rows_v2.csv",rows)
    noisy_persistence=[r for r in persistence if r["condition"]=="noise"]
    render_main(noisy_persistence,direction_terminal,site_terminal,traces,site_counts)
    render_adaptation(noisy_persistence,direction_terminal,site_terminal)
    medoids=native_examples(seeds,reps,groups)
    # All condition/window counts accompany the conditional physical point values.
    countrows=[]
    for ident,counts in sorted(site_counts.items()):
        countrows.append(dict(model=ident[0],circumference=ident[1],fraction=ident[2],window=ident[3],**counts,
                              population_curve_shown=ident[0]!="goryachev",rms=math.sqrt(counts["value"]) if math.isfinite(counts["value"]) else "unavailable"))
    write_csv(OUT/"noise_native_pair_counts_v2.csv",countrows)
    (OUT/"noise_reconciliation_v2.json").write_text(json.dumps(json_clean(reconciliation),indent=2,allow_nan=False)+"\n")
    content=make_content()
    (OUT/"noise_interpretation_v2.content.json").write_text(json.dumps(content,indent=2,ensure_ascii=False)+"\n")
    manifest=dict(generated_at=datetime.now(timezone.utc).isoformat(),analysis_only=True,new_integrations=0,
                  native_cases=243,noisy_realizations=216,seeds_per_condition=12,circumferences=[20,40,60],strengths=list(STRENGTHS),
                  source_tables=SOURCE_INFO,generator=rel(Path(__file__)),generator_sha256=sha(__file__),
                  definitions=dict(direction="Both endpoint relative first-moment magnitudes R >= 0.1; conditional cosine; coverage denominator is all within-window lag-8000 endpoint pairs across the twelve seeds.",
                    site="Both endpoints have exactly one current half-contrast region, at most one native prominence>=0.05 peak, contrast>=0.05, and R>=0.1. Shortest physical endpoint displacement modulo C; no continued ancestry assumed.",
                    intervals="Exact stored 95% percentile intervals from 1000 complete-seed resamples; correlation ratios and conditional physical means are recomputed per bootstrap draw. No time-origin pseudo-replication.",
                    physical_display="Entire Goryachev physical population curve withheld at BOTH strengths because the prespecified gate yields sparse support. All actual values/counts remain exported; no numeric zero substitutes for unavailability. No newly chosen coverage threshold changes the estimator.",
                    windows="Full[0,32000], first[0,16000], last[16000,32000], inclusive endpoints, maxlag8000, cadence5; no burn-in selected.",
                    medoids="Existing frozen within-condition four-descriptor medoid; six C20/C60 high-input histories; exact full arrays included; no selected seed used for population inference."),
                  audit_scope="Every plotted point value checked against exact source rows. All 54 noisy restricted integrals reintegrated from all1601source lags; all54lag-8000direction/site summaries independently recomputed from all216native traces; final noisy/control magnitude checked against all243native traces. Bootstrap intervals preserve audited source values and hashes; the bootstrap random stream is not rerun here.",
                  outputs={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!="noise_interpretation_v2.manifest.json"})
    (OUT/"noise_interpretation_v2.manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(dict(output=rel(OUT),native_cases=243,plotted_source_rows=len(rows),integral_checks=len(integral_check),medoids=[r["case_id"] for r in medoids],maximum_errors=reconciliation["maximum_errors"]),indent=2))


def render_retained_native():
    """Render exact medoid arrays retained in Git, with no native worker access."""
    selected=json.loads((OUT/"noise_native_examples_v2.selection.json").read_text())
    with np.load(OUT/"noise_native_examples_v2.npz",allow_pickle=False) as z:
        arrays={k:z[k].copy() for k in z.files}
    fig,axes=plt.subplots(3,2,figsize=(12.4,10.4))
    fig.subplots_adjust(left=.07,right=.90,top=.88,bottom=.12,hspace=.43,wspace=.16)
    fig.suptitle("Matched high-input native activity histories",fontsize=17,y=.975)
    fig.text(.5,.939,"C20 and C60 within-condition medoids • f = 0.01 • individual histories, never ensemble-average fields",ha="center",fontsize=10)
    for mi,model in enumerate(MODELS):
        choices=[next(r for r in selected if r["model"]==model and r["circumference"]==C) for C in (20.,60.)]
        vmax=max(float(arrays[r["compact_prefix"]+"_primary"].max()) for r in choices)
        for ci,choice in enumerate(choices):
            prefix=choice["compact_prefix"]
            primary=arrays[prefix+"_primary"]
            display=np.roll(primary,int(arrays[prefix+"_display_roll"]),axis=1)
            t=arrays[prefix+"_time"]
            C=choice["circumference"]
            ax=axes[mi,ci]
            im=ax.imshow(display,origin="lower",aspect="auto",extent=(0,C,0,32000),vmin=0,vmax=vmax,cmap="viridis",rasterized=True)
            ax.set_title(f"{LABELS[model]}   C{int(C)}   {choice['case_id']}",fontsize=10)
            ax.set_xlabel("Physical boundary arclength (shifted fixed cut)")
            ax.set_ylabel("Native time" if ci==0 else "")
            ax.set_yticks([0,16000,32000])
            ax.text(.02,.97,f"Final R = {arrays[prefix+'_magnitude'][-1]:.3f}",transform=ax.transAxes,va="top",color="white",fontsize=9,bbox=dict(facecolor="black",alpha=.4,edgecolor="none",pad=2))
        box=axes[mi,1].get_position()
        cax=fig.add_axes([.918,box.y0,.012,box.height])
        cb=fig.colorbar(im,cax=cax)
        cb.set_label("Native activity",fontsize=9,labelpad=5)
        cb.ax.tick_params(labelsize=8)
    fig.text(.07,.073,"C20/C60 share one zero-based activity scale within each model; physical widths differ. The fixed cut keeps the initial peak intact.",fontsize=9)
    fig.text(.07,.039,"Medoids minimize median/IQR-standardized distance using final initial-axis memory, minimum R, time-mean R and multisite fraction.\nThese examples illustrate native local organization; the complete seed population, magnitude gates and pair coverage determine quantitative claims.",fontsize=9,linespacing=1.5)
    savefig(fig,"noise_native_examples_v2")


def render_retained():
    """Portable default: render exclusively from the compact committed inputs."""
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"pdf.fonttype":42,"ps.fonttype":42,"axes.labelsize":9,"xtick.labelsize":9,"ytick.labelsize":9})
    with (OUT/"noise_plotted_source_rows_v2.csv").open(newline="") as h:
        rows=list(csv.DictReader(h))
    collections={name:[r for r in rows if r["source_file"]==name] for name in (
        "vision_noise_persistence.csv","vision_noise_direction_correlations.csv",
        "vision_noise_site_displacement.csv","vision_noise_ensemble_traces.csv")}
    p=collections["vision_noise_persistence.csv"]
    d=collections["vision_noise_direction_correlations.csv"]
    s=collections["vision_noise_site_displacement.csv"]
    t=collections["vision_noise_ensemble_traces.csv"]
    assert len(p)==54 and len(d)==54 and len(s)==108 and len(t)==135
    # Exact source row strings, including bootstrap intervals, are reconciled
    # independently of their union-schema display CSV. No bulk source is read.
    for filename,plotted,source_name in [
        ("noise_direction_lags_v2.csv.gz",d,"direction"),
        ("noise_site_lags_v2.csv.gz",s,"site"),
        ("noise_final_magnitude_controls_v2.csv.gz",t,"magnitude")]:
        lookup={str(r["source_csv_row"]):r for r in plotted}
        found=set()
        with gzip.open(OUT/filename,"rt",newline="") as h:
            for row in csv.DictReader(h):
                rowid=row["source_csv_row"]
                if rowid in lookup:
                    for field,value in row.items():
                        assert lookup[rowid][field]==value,(source_name,rowid,field)
                    found.add(rowid)
        assert found==set(lookup),(source_name,len(found),len(lookup))
    # Revalidate every plotted area against its retained full-resolution lag curve.
    lags=defaultdict(list)
    with gzip.open(OUT/"noise_vector_lags_v2.csv.gz","rt",newline="") as h:
        for row in csv.DictReader(h):
            lags[(*key(row),row["window"])].append((float(row["lag"]),float(row["correlation"])))
    for row in p:
        pairs=sorted(lags[(*key(row),row["window"])])
        grid,values=np.array(pairs).T
        assert np.array_equal(grid,np.arange(0,8001,5))
        np.testing.assert_allclose(np.trapezoid(values,grid),float(row["restricted_integral"]),rtol=2e-14,atol=2e-10)
    del lags
    render_main(p,d,s,t,None)
    render_adaptation(p,d,s)
    render_retained_native()
    (OUT/"noise_interpretation_v2.content.json").write_text(json.dumps(make_content(),indent=2,ensure_ascii=False)+"\n")
    manifest_path=OUT/"noise_interpretation_v2.manifest.json"
    manifest=json.loads(manifest_path.read_text())
    manifest["rendered_at"]=datetime.now(timezone.utc).isoformat()
    manifest["default_render_native_reads"]=0
    manifest["default_render_bulk_source_reads"]=0
    manifest["native_audit_mode"]="--audit-native performs optional one-time native verification/export; default uses retained compact inputs only"
    manifest["generator_sha256"]=sha(__file__)
    manifest["outputs"]={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!=manifest_path.name}
    manifest_path.write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(dict(output=rel(OUT),mode="portable retained-data render",native_reads=0,bulk_source_reads=0,restricted_area_checks=len(p)),indent=2))


def make_content():
    return dict(
        title="Activity memory increases with ring size, while physical site motion and forced-state eligibility qualify its meaning",
        assets=[dict(path="manuscript_outline/assets/noise_v2/noise_memory_geometry_v2.png",caption="Continued Figure 4: finite memory, angular geometry, native magnitude and observation coverage; both input strengths."),
                dict(path="manuscript_outline/assets/noise_v2/noise_adaptation_v2.png",caption="Supporting/continued sheet: full/first/last windows with the same lag horizon; no selected burn-in."),
                dict(path="manuscript_outline/assets/noise_v2/noise_native_examples_v2.png",caption="Supporting/continued sheet: exact high-input C20/C60 within-condition medoid histories, matched zero-based scale within each model.")],
        legend=[
            "A–E, Completed noisy-ring histories from three native periodic chemical models at circumferences C20/C40/C60. Every model/size/strength contains 12 independent realizations, each observed for 32,000 native time units. Blue/red denote the two prospectively frozen mean input fractions f=0.001/0.01. The driving law is positive lognormal activation with fixed local spatial covariance parameter 1.5 and analytically variance-normalized interpolation of OU knots (correlation time 50; spacing 5). No per-realization mean removal or circumference renormalization is applied. Local input statistics are fixed within model, while total input and angular forcing modes vary with circumference. Native amplitudes and times are not cross-model biological calibrations.",
            "A, Actual signed full-window vector-correlation area from lag 0 through 8000 on the 5-unit grid. The uncentered vector ratio retains magnitude fluctuations and weak resultants, subtracts no average vector and normalizes by both endpoint energies. The dotted 8000 line is the maximum area for correlation identically one; proximity marks unresolved decay. Intervals are stored 95% whole-seed bootstrap percentile bounds on the integrated recomputed correlation ratio. Sixteen of 18 conditions do not cross 1/e within this lag horizon. No fitted exponential time, extrapolated persistence time or cell motility time is assigned.",
            "B, Direction-only cosine at the same lag 8000 and over the same full time-origin window. Both endpoints require relative moment magnitude R>=0.1. Bars are stored pointwise 95% whole-seed bootstrap intervals. The eligible fraction of all possible endpoint pairs is shown in E. Conditional angular memory and raw vector area answer distinct questions.",
            "C, Mean relative first-moment magnitude at the final time 32000 (circles), with the matched deterministic uniform-mean input controls (dotted triangles). This final-time magnitude is a forced-state diagnostic, not the lag-correlation denominator. Initial R = 1; R can exceed one. In Goryachev high-input C20 the uniform-mean control ends with essentially zero resultant; larger rings retain substantial polarity. Its small vector area cannot be attributed solely to noise-driven wandering around an otherwise unchanged polarized state.",
            "D, Square root of the conditional mean squared shortest physical dominant-region center displacement at lag 8000; confidence bounds are monotone square roots of the source squared-displacement intervals. Both endpoints must have exactly one current half-contrast connected region, at most one peak with native prominence>=0.05, contrast>=0.05, and R>=0.1. Site descriptors are recorded every 5 native units and were independently validated at every saved 40-unit primary profile in the source audit. Endpoints do not establish uninterrupted site identity. All actual Goryachev values remain in the exported tables; the population curve is withheld at both strengths because eligibility is sparse. High input supplies only 1/0/0 eligible pair at C20/C40/C60 among 57,612 possible pairs per condition, and low-input coverage is 8.43/1.00/0.23%. The isolated one-pair result is not assigned a population precision interval or represented as zero drift.",
            "E, Fraction of all available lag 8000 endpoint pairs satisfying the direction gate (solid circles; source seed-bootstrap intervals) or the complete physical-site gate (dashed squares; exact counts exported). At high input, Holmes physical-site coverage decreases 92.1%→66.4%→37.3%, although direction coverage remains above 99%. The declared absolute-prominence gate is deliberately conservative: small secondary peaks may fail it despite a dominant incumbent. Eligibility is an observation restriction, not a functional competitor count or an instability metric.",
            "Companion adaptation sheet, full [0, 32000], first [0, 16000] and last [16000, 32000] windows retain identical maximum lag 8000 and inclusive endpoints. Symbols distinguish windows; colors retain strengths. The source estimator and all 12 seeds are unchanged, with no chosen burn-in. All windows and coverage counts remain exported. Agreement of major area ordering does not establish stationarity; physical displacement is less uniformly ordered.",
            "Companion native sheet, exact high-input C20/C60 within-condition medoids for all three models. Medoids minimize the sum of absolute deviations from within-condition descriptor medians, scaled by IQR (unit fallback when IQR ≤ 10⁻¹²), for final initial-axis memory, minimum R, time-mean R and multisite-frame fraction. Ties use case identifier. Each model uses one shared zero-based native activity scale across sizes; physical arclength axes retain circumference. The static plotting cut is shifted so the initial maximum lies at mid-boundary. The exact full-resolution arrays, shift indices, descriptor distances, native source hashes and trace arrays are included. These individual histories illustrate the quantitative population evidence without supplying independent uncertainty."
        ],
        narrative=[
            "At high input, the measured full-window restricted vector areas rise C20 → C40 → C60 in Wave-pinning 6257 → 7627 → 7885 and Holmes 4304 → 6332 → 7103 native time units. The direction-only comparisons retain the size ordering, while large-ring direction remains measurable. This supplies a direct, declared-noise activity-memory observation alongside deterministic replacement experiments; it does not close the causal migration-persistence hypothesis.",
            "The physical interpretation differs by model. High-input Wave-pinning full-window physical RMS is 3.25 [2.71,3.74],2.83 [2.40,3.29],2.25 [2.03,2.48] at C20/C40/C60, so the larger/smaller endpoint contrast includes reduced physical displacement under the declared gate. Holmes instead gives 4.57 [3.88,5.17],6.33 [5.22,7.28],5.54 [4.25,6.70] with decreasing coverage: angular memory alone cannot establish slower physical-site motion. The geometric identity Δθ=2πΔs/C explains why equal physical displacement is less angular change on a larger ring, without assigning the observed effect to geometry causally.",
            "Goryachev high-input C20 is a distinct forced-state regime, not a clean noise-specific reorientation comparison around the same maintained polarity. Conservative site eligibility is sparse even at low input and effectively absent at high input. This leaves physical comparison unavailable rather than proving stationary local activity; the raw histories can still reveal a strong dominant incumbent or locally structured activity.",
            "The lower input mostly leaves decay unresolved over lag 8000. It is retained as an informative weak-response regime; apparent areas near 8000 are bounded observation descriptors rather than arbitrarily long estimated memory times.",
            "Full/first/last restricted areas retain the major size ordering, but physical motion needs a bounded endpoint statement. Wave-pinning high-input last-half RMS is 3.19/3.49/1.93 at C20/C40/C60, so a universal monotone size rule is unsupported. Holmes C60 high-input direction changes 0.798 first-half → 0.670 last-half. These windows are shown rather than silently selecting a favorable stationary segment."
        ],
        measurements=[
            "In independent cell recordings, jointly measure whole-boundary activity vectors, their magnitudes, total reporter activity and local region/peak histories over a declared observation horizon and cadence. Report weak-magnitude fractions, pair coverage and finite-window restricted correlations; compare early and late windows before interpreting a stationary time.",
            "Register physical boundary arclength and track local activity-site displacement separately from angular reorientation. Retain multisite, low-contrast and weak-resultant intervals in the coverage denominator and report unavailable single-site comparisons explicitly. Calibrated signal is needed to interpret integrated fluorescence as molecular amount.",
            "Measure local protrusion or extension/retraction alongside reporter activity to test whether persistent or renewed activity sites correspond to protrusive behaviour. Then associate these independently measured quantities with matched cell trajectories; do not impose a velocity law on the polarity model.",
            "Relate the noisy-ring distinction to the primary cue comparison: quantify early challenger activation, incumbent loss/rebound, replacement/return/coexistence and local site lifetime after controlled cues, alongside spontaneous activity fluctuations. A local recovery response, finite cue resistance and spontaneous activity memory remain separate observables."
        ],
        limitations=[
            "The input process is a declared perturbation ensemble, not a measured cell-noise law. Changed circumference also changes total forcing and angular forcing modes, material, basal architecture and the geometric embedding, so no causal allocation among these changes is demonstrated.",
            "The restricted correlation is uncentered and may retain a nonzero average axis or adaptation. All windows remain finite; absence of a 1/e crossing gives no resolved exponential time. Confidence intervals quantify independent-seed variation at the fixed native model/input condition, not model discrepancy or biological uncertainty.",
            "Site displacement is shortest endpoint distance modulo C, not unwrapped travel, diffusion, cell movement or proof of uninterrupted site ancestry. Gate failure can be caused by small secondary peaks under the conservative prominence rule; unavailable comparisons are never zeros."
        ],
        sources=[dict(label="Audited source results",path=rel(SOURCE/"RESULTS.md")),dict(label="Frozen statistical design",path=rel(SOURCE/"DESIGN.md")),dict(label="Source/native reconciliation and exact pair counts",path="manuscript_outline/assets/noise_v2/noise_reconciliation_v2.json"),dict(label="All plotted source rows and hashes",path="manuscript_outline/assets/noise_v2/noise_plotted_source_rows_v2.csv"),dict(label="Compact native histories and selection",path="manuscript_outline/assets/noise_v2/noise_native_examples_v2.selection.json"),dict(label="Complete compact lag exports and provenance",path="manuscript_outline/assets/noise_v2/noise_interpretation_v2.manifest.json")])


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-native",action="store_true",help="Optional one-time full native/source reconciliation and compact export; needs excluded local worker archives")
    args=parser.parse_args()
    if args.audit_native:
        audit_native()
    render_retained()
