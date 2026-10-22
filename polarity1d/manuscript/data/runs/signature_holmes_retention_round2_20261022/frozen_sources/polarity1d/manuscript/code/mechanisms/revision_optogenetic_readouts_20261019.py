"""Small finite-response discussion sketch from saved native fields; no solver."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[4]
ROOT = REPO / "polarity1d/manuscript"
F2 = ROOT / "figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2"
OUT = F2 / "optogenetic_revision/discussion_sheet"
REC = ROOT / "data/runs/vision_linear_20261017/recovery"
MODEL_NAMES = {"wave_pinning": "Wave-pinning", "goryachev": "Goryachev"}
COLORS = {19: "#d27828", 39: "#347ba7"}
EPS = (.001, .005, .01)
SOURCES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source(path):
    path = Path(path)
    SOURCES[str(path.relative_to(REPO))] = sha(path)
    return path


def rows(path):
    with source(path).open() as stream:
        return list(csv.DictReader(stream))


def csvwrite(name, values):
    path = OUT / "supporting/data" / name
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)
    return path


def selected_activation():
    records = rows(F2 / "supporting/data/all127_native_activation_and_fate.csv")
    result, exported = {}, []
    for model in MODEL_NAMES:
        for size in (19, 39):
            r, = [r for r in records if r["model"] == model and float(r["size_native"]) == size and r["primary"] == "True"]
            p = source(REPO / r["archive"])
            assert sha(p) == r["archive_sha256"]
            meta = json.loads(source(p.with_suffix(".json")).read_text())
            c = meta["job"]["case"]
            assert c["width"] == 3. and meta["job"]["duration"] == 1000.
            assert c["response_amplitude"] == (.1 if model == "wave_pinning" else 6.4)
            with np.load(p) as z:
                t = z["cue_time"].copy()
                u0 = z["originalstate"][0].copy()
                q = z["x"] / size
                threshold = u0.min() + .5*np.ptp(u0)
                j = np.flatnonzero(u0 < threshold)[0]
                edge = q[j-1] + (threshold-u0[j-1])/(u0[j]-u0[j-1])*(q[j]-q[j-1])
                newmask = q >= 1-edge
                peak = z["cue_field_active"][:, newmask].max(axis=1)
                np.testing.assert_allclose(peak, z["cue_regional_challenger_peak"], rtol=0, atol=1e-12)
                delta = peak-peak[0]
                assert t[0] == 0. and t[1] == .1 and t[-1] == 1000.
                result[model, size] = dict(time=t, peak=peak, delta=delta, initial=u0,
                                          archive=str(p.relative_to(REPO)), amplitude=c["response_amplitude"])
                for i, time in enumerate(t):
                    exported.append(dict(model=model, size_native=size, time_native=float(time),
                                         challenger_peak_native=float(peak[i]), initial_peak_in_same_region=float(peak[0]),
                                         peak_increment_native=float(delta[i]), native_input_peak=c["response_amplitude"],
                                         physical_support_native=3., cue_duration_native=1000.,
                                         region="fixed original basal half-contrast core mirrored to challenger end",
                                         native_channel="stimulus*inactive transferred to active",
                                         source_archive=str(p.relative_to(REPO))))
    return result, exported


def selected_recovery(activation):
    summaries = rows(REC / "results/recovery_summary.csv")
    gates = rows(REC / "results/recovery_candidate_gates.csv")
    result, exported, errors = {}, [], []
    for model in MODEL_NAMES:
        for size in (19, 39):
            s, = [r for r in summaries if r["model"] == model and float(r["size_native"]) == size]
            p = source(REC / s["archive"])
            assert sha(p) == s["archive_sha256"]
            status = json.loads(source(p.with_name("status.json")).read_text())
            assert status["complete"] and status["zero_kick_full_field_identity"]
            with np.load(p) as z:
                t = z["tight_sham_time"].copy()
                a = z["tight_sham_core_amount"].copy()
                ar = z["tighter_sham_core_amount"]
                np.testing.assert_array_equal(z["tight_sham_field_active"][0], activation[model, size]["initial"])
                recomputed_sham = z["tight_sham_field_active"] @ z["core_weights"]
                errors.append(float(abs(recomputed_sham-a).max()))
                np.testing.assert_allclose(recomputed_sham, a, rtol=0, atol=2e-12)
                series = {}
                for eps in EPS:
                    label = f"eps{eps:g}"
                    probe = z["tight_"+label+"_core_amount"].copy()
                    pr = z["tighter_"+label+"_core_amount"]
                    recomputed_probe = z["tight_"+label+"_field_active"] @ z["core_weights"]
                    errors.append(float(abs(recomputed_probe-probe).max()))
                    np.testing.assert_allclose(recomputed_probe, probe, rtol=0, atol=2e-12)
                    d, dr = a-probe, ar-pr
                    f, fr = 1-d/d[0], 1-dr/dr[0]
                    roundoff = 64*np.finfo(float).eps*(abs(a)+abs(probe)+abs(a[0])+abs(probe[0]))/d[0]
                    e = np.maximum.reduce((2*abs(f-fr), 2*abs(a-z["tight_eps0_core_amount"])/d[0], roundoff))
                    series[eps] = dict(probe=probe, deficit=d, restored=f, empirical_E=e)
                    for i, time in enumerate(t):
                        exported.append(dict(model=model, size_native=size, epsilon=eps, time_native=float(time),
                                             sham_core_amount_native=float(a[i]), probe_core_amount_native=float(probe[i]),
                                             initial_damage_native=float(d[0]), initial_damage_fraction_of_initial_core=float(d[0]/a[0]),
                                             core_deficit_native=float(d[i]), restored_fraction=float(f[i]),
                                             sham_change_native=float(a[i]-a[0]), probe_change_native=float(probe[i]-probe[0]),
                                             empirical_numerical_sensitivity=float(e[i]),
                                             inside_qualified_early_window=bool(time <= float(s["horizon_native"])),
                                             source_archive=str(p.relative_to(REPO))))
                clock = float(s["horizon_native"])
                at = int(np.argmin(abs(t-clock)))
                assert t[0] == 0. and t[-1] == 1. and abs(t[at]-clock) < 1e-14
                assert abs(series[.005]["restored"][at]-float(s["restored_fraction"])) < 1e-12
                selected_gate, = [r for r in gates if r["model"] == model and float(r["size_native"]) == size and abs(float(r["horizon_native"])-clock) < 1e-14]
                assert selected_gate["passes"] == "True"
                late_gate, = [r for r in gates if r["model"] == model and float(r["size_native"]) == size and float(r["horizon_native"]) == 1.]
                assert late_gate["linearity_pass"] == "True" and late_gate["saturation_pass"] == "False"
                result[model, size] = dict(time=t, sham=a, series=series, clock=clock, core_samples=int(s["core_sample_centers"]),
                                          primary_initial_fraction=float(s["initial_damage_fraction_of_core_active"]),
                                          late_gate=late_gate, selected_gate=selected_gate,
                                          archive=str(p.relative_to(REPO)))
    return result, exported, max(errors)


def audit_windows(a, b):
    activation, recovery, amplitudes, raw_crossings = [], [], [], []
    for model in MODEL_NAMES:
        time = a[model, 19]["time"]
        np.testing.assert_array_equal(time, a[model, 39]["time"])
        ratio = a[model, 39]["delta"][1:]/a[model, 19]["delta"][1:]
        for lo, hi in ((.1, 1.), (.1, 10.), (.1, 1000.)):
            mask = (time[1:] >= lo) & (time[1:] <= hi)
            activation.append(dict(model=model, window_start_native=lo, window_end_native=hi,
                                   observed_large_small_peak_increment_ratio_min=float(ratio[mask].min()),
                                   observed_large_small_peak_increment_ratio_max=float(ratio[mask].max()),
                                   ratio_at_window_end=float(ratio[np.flatnonzero(time[1:] == hi)[0]])))
        t = b[model, 19]["time"]
        np.testing.assert_array_equal(t, b[model, 39]["time"])
        for h in sorted(set((.001,.01,.1,1.,b[model,19]["clock"]))):
            at = int(np.argmin(abs(t-h)))
            assert abs(t[at]-h) < 1e-14
            s, l = [b[model,size]["series"][.005] for size in (19,39)]
            recovery.append(dict(model=model, time_native=float(t[at]),
                                 previous_common_qualified_clock=bool(abs(t[at]-b[model,19]["clock"]) < 1e-14),
                                 restored_fraction_small=float(s["restored"][at]), restored_fraction_large=float(l["restored"][at]),
                                 restored_fraction_large_over_small=float(l["restored"][at]/s["restored"][at]),
                                 raw_deficit_small_native=float(s["deficit"][at]), raw_deficit_large_native=float(l["deficit"][at])))
        difference = b[model,39]["series"][.005]["deficit"]-b[model,19]["series"][.005]["deficit"]
        cross = np.flatnonzero(difference <= 0.)
        raw_crossings.append(dict(model=model, observed_crossing=bool(len(cross)),
                                  last_sample_before_crossing_native=float(t[cross[0]-1]) if len(cross) else None,
                                  first_sample_after_crossing_native=float(t[cross[0]]) if len(cross) else None,
                                  interpolation_used=False, followup_native=float(t[-1])))
        for size in (19,39):
            c = b[model,size]
            f = np.stack([c["series"][eps]["restored"] for eps in EPS])
            at = int(np.argmin(abs(t-c["clock"])))
            amplitudes.append(dict(model=model,size_native=size,epsilon_values="0.001/0.005/0.01",
                                   initial_primary_damage_native=float(c["series"][.005]["deficit"][0]),
                                   initial_primary_fraction_of_core=c["primary_initial_fraction"],
                                   qualified_clock_native=c["clock"],F_at_clock_min=float(f[:,at].min()),F_at_clock_max=float(f[:,at].max()),
                                   F_at_1_min=float(f[:,-1].min()),F_at_1_max=float(f[:,-1].max()),
                                   maximum_absolute_F_span_over_time_0p0001_to_1=float(np.ptp(f[:,t>=1e-4],axis=0).max()),
                                   paired_amplitude_linearity_gate_pass_to_1=c["late_gate"]["linearity_pass"] == "True",
                                   early_restoration_under20pct_gate_pass_to_1=False))
    return activation, recovery, amplitudes, raw_crossings


def sketch(a, b):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,
                         "axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
    fig, axes = plt.subplots(2,2,figsize=(10,8))
    fig.subplots_adjust(left=.095,right=.975,bottom=.135,top=.82,hspace=.61,wspace=.28)
    fig.text(.5,.97,"Two measured-response options for discussion",ha="center",fontsize=17,weight="bold")
    fig.text(.5,.929,"Retained starting contrast: WP mainly widens; Goryachev mainly grows taller. Original all-five artwork is linked unchanged.",ha="center",fontsize=9.2)
    fig.text(.065,.877,"A   Local excitation → challenger activity",fontsize=12,weight="bold")
    fig.text(.065,.851,"Native source = stimulus × inactive; support 3; peaks WP 0.1 / Goryachev 6.4. One fixed dose per model.",fontsize=9)
    fig.text(.065,.478,"B   Incumbent deactivation → original-core deficit and recovery",fontsize=12,weight="bold")
    fig.text(.065,.452,"Instantaneous 0.5% tapered basal-excess direction moves active into inactive; no countercue; paired sham.",fontsize=9)
    for col, model in enumerate(MODEL_NAMES):
        ax = axes[0,col]
        for size in (19,39):
            c=a[model,size]
            ax.plot(c["time"][1:],c["delta"][1:],color=COLORS[size],ls="-" if size == 19 else "--",lw=1.7,label=f"L={size}")
            ax.scatter(c["time"][1],c["delta"][1],s=23,color=COLORS[size],zorder=4)
        ax.set_xscale("log");ax.set_yscale("log")
        ax.set_xlim(.1,1000);ax.set_xticks((.1,1,10,100,1000),labels=("0.1","1","10","100","1000"))
        ax.set_ylim(.01,2. if model == "wave_pinning" else 60.)
        ax.set_ylabel("Peak increment in fixed challenger region\n(native concentration; log scale)",fontsize=8.4)
        ax.set_xlabel("Time after cue onset (native units; log scale)",fontsize=8.4)
        ax.set_title(MODEL_NAMES[model],fontsize=10.5,weight="bold",loc="left")
        ax.legend(frameon=False,loc="lower right",fontsize=8)
        r=a[model,39]["delta"][1]/a[model,19]["delta"][1]
        ax.text(.025,.965,f"First saved interval 0–0.1:\nΔpeak {a[model,19]['delta'][1]:.5f} / {a[model,39]['delta'][1]:.5f}\nL39/L19 = {r:.4f}",transform=ax.transAxes,va="top",fontsize=8,
                bbox=dict(facecolor="white",edgecolor="none",alpha=.9,pad=2))
        ax.grid(alpha=.13)
        ax=axes[1,col]
        clock=b[model,19]["clock"]
        ax.axvspan(0,clock,color="#8b76ac",alpha=.09)
        ax.axvline(clock,color="#77628e",ls=":",lw=1)
        for size in (19,39):
            c=b[model,size];s=c["series"][.005]
            ax.plot(c["time"],s["deficit"],color=COLORS[size],ls="-" if size == 19 else "--",lw=1.7,
                    label=rf"L={size}; $D_0$={s['deficit'][0]:.5f}")
            ax.scatter([0,1],[s["deficit"][0],s["deficit"][-1]],color=COLORS[size],s=17,zorder=4)
        ax.set_xscale("symlog",linthresh=.001,linscale=.5)
        ax.set_xlim(0,1);ax.set_xticks((0,.001,.01,.1,1),labels=("0","0.001","0.01","0.1","1"))
        ax.set_ylim(bottom=0,top=.045 if model == "wave_pinning" else .36)
        ax.set_xlabel("Time after kick (native; linear to 0.001, then log)",fontsize=8.4)
        ax.set_ylabel("Sham − probe amount in original core\n(native concentration × length)",fontsize=8.4)
        ax.legend(frameon=False,loc="upper left",fontsize=8)
        at=np.argmin(abs(b[model,19]["time"]-clock))
        f19=b[model,19]["series"][.005]["restored"];f39=b[model,39]["series"][.005]["restored"]
        ax.text(.98,.98,f"Qualified τ={clock:.6g}\nF(L19, L39) at τ: ({f19[at]:.3f}, {f39[at]:.3f})\nAt t=1: ({f19[-1]:.3f}, {f39[-1]:.3f})\nF₃₉/F₁₉: {f39[at]/f19[at]:.3f} → {f39[-1]/f19[-1]:.3f}",transform=ax.transAxes,va="top",ha="right",fontsize=8,
                bbox=dict(facecolor="white",edgecolor="none",alpha=.92,pad=2))
        ax.grid(alpha=.13)
    fig.text(.5,.018,"B shading marks the previously qualified early estimator window (≤20% restoration); later lines are valid finite-time responses.\n"
             "Native doses and clocks have no shared optical or reporter calibration. These are conditional readout options, not a replacement Figure 2 layout.",ha="center",va="bottom",fontsize=8.5,linespacing=1.5)
    paths=[]
    for ext in ("png","pdf","svg"):
        p=(OUT if ext == "png" else OUT/"supporting/vectors")/("finite_response_options."+ext)
        fig.savefig(p,dpi=240,facecolor="white");paths.append(p)
    plt.close(fig)
    return paths


def docs(audit, recovery_windows, amplitudes, cross):
    original = source(REPO/"current manuscript figs/initial_polarity_length.png")
    canonical = source(ROOT/"figures/initial_polarity_length/initial_polarity_length.png")
    assert sha(original) == sha(canonical) == "b94d05c1caabeb136a2cb744782ad56ac4e366bf2a00ca266caebd6af534f9ba"
    organization = "<"+os.path.relpath(original,OUT)+">"
    readme = f"""# Native finite-response options for author discussion

This is a small evidence/readout discussion sheet, not a replacement Figure2 or an optical-actuator choice. The reviewed prototype, earlier options note and all original artwork remain unchanged. The central question remains whether similar reversal-size trends arise through mechanisms predicting different challenge responses.

## Retained organization: unchanged original evidence

[Original all-five organization figure]({organization}). Its exact97 basal measurements distinguish broader poles from the taller Goryachev peak. Over L19→39, WP amount/width/peak folds are2.124/2.042/1.000; Goryachev folds are2.152/1.119/1.903. This passive comparison motivates the conditional challenges below; it does not by itself establish a reversal mechanism.

![Unchanged original incumbent organization artwork]({organization})

## Two finite-response options

![Observed activation and paired-core recovery](finite_response_options.png)

[PDF](supporting/vectors/finite_response_options.pdf) · [SVG](supporting/vectors/finite_response_options.svg) · [Caption](supporting/metadata/CAPTION.md)

**Excitation:** taller, substrate-depleting organization → a local native coefficient multiplying inactive material → weaker measured challenger activity. Both highlighted models add stimulus×inactive to active and subtract it from inactive. A GEF-recruitment actuator, direct activity forcing and inhibition have not been shown to implement the same channel. Inactive-independent forcing would lack this instantaneous substrate factor. The curves are actual saved field changes, including all native reaction and transport, not the calculated input gain.

**Deactivation:** different organization/coupled kinetics → a declared local active-to-inactive loss → different evolution of the original-core activity deficit relative to a paired sham. The graph uses actual amounts and their absolute initial deficits D₀. It does not plot an instantaneous reaction derivative. This controlled recovery cannot isolate feedback from transport or uniquely allocate reversal resistance.

## What survives a wider observation window?

The activation contrast is not confined to t=0.1: over all saved points from0.1 to10, the L39/L19 peak-increment ratio is0.99718–0.99808 in WP and0.65089–0.65325 in Goryachev. Across the entire1000-unit cue it remains below1:0.96652–0.99808 versus0.25344–0.65325. The late nonlinear response therefore preserves the direction of the contrast, not a fixed fold. The first post-cue sample is0.1; the interval0–0.1 is measured only by its two endpoints. There is no matched excitation-amplitude sensitivity in this selected fixed-support set; computed gain cannot fill that gap.

The recovery contrast also persists through the complete stored1-unit follow-up, but its fold weakens:

| Model | Qualified early clock | F at that clock, L19 / L39 | F at t=0.1, L19 / L39 | F at t=1, L19 / L39 |
| --- | ---: | ---: | ---: | ---: |
| WP |0.177828|0.1339 /0.1373|0.07831 /0.07981|0.5136 /0.5563|
| Goryachev |0.00562341|0.05404 /0.1859|0.3749 /0.7648|0.6844 /0.8764|

Goryachev's F₃₉/F₁₉ fold drops3.439→2.040→1.281 over its qualified clock→0.1→1. WP's fold is1.026 at its clock and1.083 at1. In native amount units the initially larger Goryachev deficit crosses below the smaller-cell deficit between the saved times0.047315 and0.050119. This is an observed crossover bracket, without interpolation.

The0.1%, 0.5%, 1% kick sensitivities retain amplitude linearity through1 for all four cases. Their maximum absolute spread in F over10⁻⁴–1 is0.00165 in WP and0.00302 in Goryachev. Later restoration exceeds the20% qualification gate for the early estimator; that does not invalidate the later native trajectories. The qualified clocks were the largest common tested pilot windows satisfying numerical/linearity/early-restoration gates, not experimentally chosen sampling times.

## Experimental inference remains conditional

These are finite-time model signatures with documented wider-window persistence, rather than only instantaneous algebraic diagnostics. They are not yet validated optical protocols. The primary deactivation removes only0.217–0.223% of WP core activity and0.285% of Goryachev core activity: no reporter precision or response time has been established to resolve those changes. Native time has no shared experimental calibration; a slow reporter or actuator could average away the earliest contrast. The deactivation is an instantaneous state kick, not a simulated finite optical pulse. Goryachev L39 has four native sample centres in its original core; the prior checks do not establish spatial convergence. There is no active-loss follow-up beyond1 in these archives. Zero-kick and tighter-solver comparisons qualify the saved numerical response, not experimental measurability.

Thus the useful discussion choice is which characterized perturbation and reporter could test the distinct conditional responses. The present data support that question; they do not select the actuator or prove a unique molecular explanation.

[All native excitation reads](supporting/data/activation_cue_curves.csv) · [All paired sham/probe amount traces and three amplitudes](supporting/data/recovery_native_curves.csv) · [Window comparisons](supporting/data/recovery_window_contrasts.csv) · [Amplitude sensitivity](supporting/data/recovery_amplitude_sensitivity.csv) · [Sources, checks and missing flags](supporting/metadata/validation.json)
"""
    # Keep prose spacing distinct from compact generated CSV values.
    fixes={"Figure2":"Figure 2","exact97":"exact 97","are2.124":"are 2.124","are2.152":"are 2.152",
           "from0.1":"from 0.1","to10":"to 10","is0.99718":"is 0.99718","and0.65089":"and 0.65089",
           "entire1000":"entire 1000","below1:0.96652":"below 1: 0.96652","versus0.25344":"versus 0.25344",
           "is0.1":"is 0.1","interval0–0.1":"interval 0–0.1","stored1-unit":"stored 1-unit",
           "drops3.439":"drops 3.439","is1.026":"is 1.026","and1.083":"and 1.083","at1.":"at 1.",
           "times0.047315":"times 0.047315","and0.050119":"and 0.050119","The0.1%":"The 0.1%",
           "through1":"through 1","over10⁻⁴–1":"over 10⁻⁴–1","is0.00165":"is 0.00165","and0.00302":"and 0.00302",
           "the20%":"the 20%","only0.217":"only 0.217","and0.285%":"and 0.285%","beyond1":"beyond 1"}
    for x,y in fixes.items():readme=readme.replace(x,y)
    (OUT/"README.md").write_text(readme)
    caption = r"""# Finite-response discussion sketch

Two response options, each faceted by model, from the exact retained L19/L39 corrected basal checkpoints. Original all-five incumbent organization artwork is reused unchanged in the discussion entrypoint. No reversal-duration panel or activation/fate scatter is added.

A: Native fixed-support3 activation, peak0.1 in WP and6.4 in Goryachev, held for1000 model-native units. The input transfers stimulus×inactive into active. Curves give the actual saved peak increment from its initial value in the mirrored original half-contrast core. Both axes are logarithmic; only positive post-cue samples are drawn. Baseline increment is zero at t=0, and the first saved post-cue sample is t=0.1. Increments are native concentrations, without independent peak scaling. Annotations report actual0–0.1 changes, not a derivative or reporter response. Neither native doses nor time units are matched optical calibrations.

B: Uncued paired sham and probe after an instantaneous kick removing0.5% of the smooth frozen basal-excess direction, adding it pointwise to inactive. The actual whole-core loss is about0.22% in WP and0.285% in Goryachev. Curves show the paired core-amount deficit in native integrated-activity units; D₀ is explicit in each legend. All native fields evolve. Time is linear through0.001 then logarithmic; t=0 is retained. Shading/dotted lines mark the previously qualified common early estimator clock for each model. The restored fraction F is reported at that clock and the last stored time1. Later trajectories remain valid finite-time observations even when restoration exceeds the early estimator's20% qualification gate. The0.1%/1% amplitude and tighter-numerics checks are linked separately. Recovery after1 is unavailable in this assay; the finite-cue13,000-unit withdrawal histories are a different experiment and do not fill it.

The measurement definitions are

$$
D(t)=A_{\mathrm{sham}}(t)-A_{\mathrm{probe}}(t),\qquad
F(t)=\frac{D(0)-D(t)}{D(0)}.
\tag{1}
$$

This is a discussion sketch, not an accepted Figure2 layout or optical-actuator choice. Initial gain is an analytic input-source derivative; these plots instead use finite saved fields/amounts. Their experimentally resolvable windows, reporter kinetics and actuator mapping are unknown. Stronger restoration after this one small disturbance does not mean universal stability or causally allocate size-dependent reversal resistance.
"""
    for x,y in {"support3":"support 3","peak0.1":"peak 0.1","and6.4":"and 6.4","for1000":"for 1000","actual0–0.1":"actual 0–0.1","removing0.5%":"removing 0.5%","about0.22%":"about 0.22%","and0.285%":"and 0.285%","through0.001":"through 0.001","time1":"time 1","estimator's20%":"estimator's 20%","The0.1%":"The 0.1%","after1":"after 1","finite-cue13,000":"finite-cue 13,000","Figure2":"Figure 2"}.items():caption=caption.replace(x,y)
    (OUT/"supporting/metadata/CAPTION.md").write_text(caption)


def main():
    for folder in (OUT,OUT/"supporting/data",OUT/"supporting/metadata",OUT/"supporting/vectors"):
        folder.mkdir(parents=True,exist_ok=True)
    unchanged={str(p):sha(p) for p in F2.rglob("*") if p.is_file() and OUT not in p.parents}
    a, ar=selected_activation()
    b, br, core_error=selected_recovery(a)
    aw, rw, amps, cross=audit_windows(a,b)
    artifacts=[csvwrite(name,value) for name,value in (("activation_cue_curves.csv",ar),("recovery_native_curves.csv",br),
                ("activation_window_contrasts.csv",aw),("recovery_window_contrasts.csv",rw),
                ("recovery_amplitude_sensitivity.csv",amps),("raw_deficit_crossing_brackets.csv",cross))]
    artifacts.extend(sketch(a,b))
    docs(aw,rw,amps,cross)
    artifacts.extend((OUT/"README.md",OUT/"supporting/metadata/CAPTION.md"))
    source(Path(__file__))
    source(ROOT/"code/mechanisms/all_five_incumbent_recovery_20261017.py")
    source(ROOT/"data/runs/vision_linear_20261017/recovery/results/verification.json")
    assert all(sha(Path(p)) == h for p,h in unchanged.items())
    missing=dict(recovery_after_native_time_1="not saved/unmeasured for this active-loss assay",
                 activation_before_first_postcue_sample_0p1="only baseline and t0.1 endpoints saved; no interior transient samples",
                 activation_matched_dose_amplitude_sensitivity="unmeasured in selected fixed-support native comparison",
                 optical_actuator_mapping="not selected or validated",reporter_kinetics_and_noise="not modeled or calibrated")
    validation=dict(generated_at_utc=datetime.now(timezone.utc).isoformat(),passed=True,no_integrations=True,
                    scope="standalone readout discussion options; not a replacement Figure2 composition",
                    source_sha256=SOURCES,artifact_sha256={str(p.relative_to(OUT)):sha(p) for p in artifacts},
                    preserved_earlier_figure2_files=len(unchanged),core_amount_recomputed_from_saved_fields_max_abs_error=core_error,
                    activation_curve_rows=len(ar),recovery_curve_rows=len(br),activation_window_comparisons=aw,
                    recovery_window_comparisons=rw,amplitude_sensitivity=amps,raw_deficit_crossing_brackets=cross,
                    wider_window_interpretation="Later restoration >20% fails only the early-estimator qualification gate, not finite native dynamics.",
                    missing_or_unvalidated=missing)
    (OUT/"supporting/metadata/validation.json").write_text(json.dumps(validation,indent=2)+"\n")
    print(json.dumps(dict(output=str(OUT),activation_rows=len(ar),recovery_rows=len(br),
                          maximum_core_amount_recomputation_error=core_error,prior_files_unchanged=len(unchanged)),indent=2))


if __name__ == "__main__":
    main()
