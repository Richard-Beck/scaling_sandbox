"""Separate provisional Figure 2 candidate: retained native finite responses only."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

REPO = Path(__file__).resolve().parents[4]
ROOT = REPO / "polarity1d/manuscript"
F2 = ROOT / "figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2"
OUT = F2 / "finite_response_candidate_20261019"
LABELS = {"wave_pinning":"Wave-pinning", "goryachev":"Goryachev", "debelly":"deBelly", "holmes_model3":"Holmes 3", "spring":"Spring"}
COLORS = dict(zip(LABELS,("#0072b2","#d55e00","#009e73","#7c5ba7","#a6783d")))
REFERENCES = dict(zip(LABELS,(19.,19.,19.,20.,1.)))
PURPLE = "#9760a7"
SOURCES = {}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def source(p):
    p=Path(p); SOURCES[str(p.relative_to(REPO))]=sha(p); return p


def rows(p):
    with source(p).open() as f:return list(csv.DictReader(f))


def csvwrite(name, values):
    p=OUT/"supporting/data"/name
    with p.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(values[0]));w.writeheader();w.writerows(values)
    return p


def end_width(x,u,threshold,length):
    """Boundary-connected fixed-threshold width; linear native contour crossing."""
    if u[0] < threshold:return 0.
    below=np.flatnonzero(u<threshold)
    if not len(below):return length
    j=below[0]
    return float(x[j-1]+(threshold-u[j-1])/(u[j]-u[j-1])*(x[j]-x[j-1]))


def native_cases():
    records=rows(F2/"supporting/data/all127_native_activation_and_fate.csv")
    cases={};traces=[];snapshots=[];summary=[]
    for model in ("wave_pinning","goryachev"):
        for length in (19,39):
            row,=[r for r in records if r["model"]==model and float(r["size_native"])==length and r["primary"]=="True"]
            p=source(REPO/row["archive"]);assert sha(p)==row["archive_sha256"]
            meta=json.loads(source(p.with_suffix(".json")).read_text());c=meta["job"]["case"]
            assert c["width"]==3. and meta["job"]["duration"]==1000.
            with np.load(p) as z:
                x=z["x"].copy();initial=z["originalstate"].copy();u0=initial[0]
                threshold=float(u0.min()+.5*np.ptp(u0))
                phase_times=[z["cue_time"].copy(),z["post_time"][1:].copy()]
                t=np.r_[phase_times[0],1000+phase_times[1]]
                u=np.vstack((z["cue_field_active"],z["post_field_active"][1:]))
                oldpeak=np.r_[z["cue_regional_old_peak"],z["post_regional_old_peak"][1:]]
                regions=np.r_[z["cue_regional_sites_above_initial_halfcontrast"],z["post_regional_sites_above_initial_halfcontrast"][1:]]
                assert t[0]==0. and t[-1]==14000.
                assert np.array_equal(z["cue_field_active"][-1],z["post_field_active"][0])
                withdrawal=len(phase_times[0])-1
                assert t[withdrawal]==1000.
                cue_mask=x>=length-3
                cuepeak=u[:,cue_mask].max(axis=1)
                cell_edges=np.linspace(0,length,len(x)+1)
                weights=np.maximum(0.,np.minimum(cell_edges[1:],length)-np.maximum(cell_edges[:-1],length-3))
                assert abs(weights.sum()-3.)<1e-12
                cue_excess=(u-u0)@weights
                left=np.array([end_width(x,f,threshold,length) for f in u])
                right=np.array([end_width(length-x[::-1],f[::-1],threshold,length) for f in u])
                final=np.stack((z["post_field_active"][-1],z["post_field_inactive"][-1]))
                scales=np.array([c["fixed_scales"]["active"],c["fixed_scales"]["inactive"]])[:,None]
                original_error=float(np.max(abs(final-initial)/scales))
                reflected_error=float(np.max(abs(final-initial[:,::-1])/scales))
                after=z["post_regional_sites_above_initial_halfcontrast"]
                if model=="wave_pinning" and length==39:assert np.all(after==2)
                input_shape=np.where(length-x<=3.,.5*(1+np.cos(np.pi*np.clip((length-x)/3.,0,1))),0.)
                case_id=f"{model}_L{length}"
                cases[model,length]=dict(time=t,u=u,x=x,threshold=threshold,withdrawal=withdrawal,cuepeak=cuepeak,
                                         oldpeak=oldpeak,left_width=left,right_width=right,regions=regions,
                                         initial=initial,cue_excess=cue_excess,case_id=case_id)
                for i,time in enumerate(t):
                    traces.append(dict(case_id=case_id,time_since_cue_onset_native=float(time),cue_on=bool(time<=1000),
                                       time_after_withdrawal_native=float(time-1000) if time>=1000 else None,
                                       peak_in_physical_cue_roi_native=float(cuepeak[i]),old_original_core_peak_native=float(oldpeak[i]),
                                       cue_roi_activity_increment_native=float(cue_excess[i]),
                                       fixed_initial_threshold_native=threshold,initial_side_width_native=float(left[i]),
                                       challenged_side_width_native=float(right[i]),extent_beyond_cue_native=max(0.,float(right[i]-3)),
                                       connected_regions_above_initial_halfcontrast=int(regions[i])))
                for state,idx in (("initial",0),("withdrawal",withdrawal),("final",len(t)-1)):
                    if length==39:
                        for j in range(len(x)):
                            snapshots.append(dict(case_id=case_id,state=state,time_since_cue_onset_native=float(t[idx]),
                                                  physical_position_native=float(x[j]),active_concentration_native=float(u[idx,j])))
                summary.append(dict(case_id=case_id,model=model,length_native=length,native_input_peak=c["response_amplitude"],
                                    requested_physical_support=3.,sampled_input_peak=float(c["response_amplitude"]*input_shape.max()),
                                    sampled_input_physical_integral=float(c["response_amplitude"]*input_shape.sum()*length/len(x)),
                                    cue_roi_sample_centres=int(cue_mask.sum()),cue_duration_native=1000.,withdrawal_followup_native=13000.,
                                    native_input_channel="stimulus*inactive transferred conservatively to active",
                                    native_cells=len(x),initial_threshold=threshold,initial_peak=float(u0.max()),
                                    initial_side_width=float(left[0]),cue_end_cue_roi_peak=float(cuepeak[withdrawal]),
                                    cue_end_old_peak=float(oldpeak[withdrawal]),cue_end_initial_side_width=float(left[withdrawal]),
                                    cue_end_challenged_side_width=float(right[withdrawal]),final_cue_roi_peak=float(cuepeak[-1]),
                                    final_old_peak=float(oldpeak[-1]),final_initial_side_width=float(left[-1]),final_challenged_side_width=float(right[-1]),
                                    final_resolved_regions=int(regions[-1]),native_field_tail_pass=row["field_tail_pass"],
                                    full_two_field_endpoint_error_original=original_error,full_two_field_endpoint_error_reflected=reflected_error,
                                    archive=str(p.relative_to(REPO)),archive_sha256=sha(p),checkpoint=c["checkpoint"],checkpoint_sha256=c["checkpoint_sha256"]))
    return cases,traces,snapshots,summary


def organization():
    table=ROOT/"data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv"
    native=rows(table);assert len(native)==97
    provenance=json.loads(source(table.with_name("report_provenance.json")).read_text())
    assert sha(table)==provenance["summary_csv_sha256"]
    result=[];excluded=[]
    for model in LABELS:
        rr=[r for r in native if r["model"]==model]
        reference,=[r for r in rr if float(r["domain_length"])==REFERENCES[model]]
        for r in rr:
            result.append(dict(model=model,size_native=float(r["domain_length"]),reference_size_native=REFERENCES[model],
                               size_over_reference=float(r["domain_length"])/REFERENCES[model],
                               excess_amount_native=float(r["old_pole_excess_active"]),width_native=float(r["old_pole_width"]),peak_native=float(r["old_pole_peak_active"]),
                               excess_amount_over_reference=float(r["old_pole_excess_active"])/float(reference["old_pole_excess_active"]),
                               width_over_reference=float(r["old_pole_width"])/float(reference["old_pole_width"]),
                               peak_over_reference=float(r["old_pole_peak_active"])/float(reference["old_pole_peak_active"]),
                               checkpoint=r["checkpoint"],checkpoint_sha256=r["checkpoint_sha256"]))
        for size in provenance["excluded_lengths"][model]:excluded.append(dict(model=model,size_native=size,status="preparation_ineligible"))
    assert len(excluded)==13
    return result,excluded


def finite_figure(cases):
    fig,axes=plt.subplots(2,2,figsize=(11.2,8.6))
    fig.subplots_adjust(left=.085,right=.975,bottom=.17,top=.81,hspace=.60,wspace=.23)
    fig.text(.5,.962,"Distinct finite responses to a local polarity challenge",ha="center",fontsize=17,weight="bold")
    fig.text(.5,.913,"Requested support 3 and native peak fixed across sizes within each model · 1000 cue units + 13,000 cue-free units",ha="center",fontsize=9.3)
    fig.text(.5,.876,"Measured activity in the physical cue region, with an L19 reference; all main histories retain native concentration units.",ha="center",fontsize=9.5)
    for col,model in enumerate(("wave_pinning","goryachev")):
        color=COLORS[model];c=cases[model,39];ref=cases[model,19];w=c["withdrawal"]
        ax=axes[0,col]
        ax.axvspan(0,1000,color=PURPLE,alpha=.07)
        ax.axvline(1000,color=PURPLE,lw=.9,ls=":")
        ax.plot(ref["time"],ref["cuepeak"],color=".52",lw=1.2,ls="--",label="L19 reference")
        ax.plot(c["time"],c["cuepeak"],color=color,lw=2.,label="L39")
        ax.scatter([1000,14000],c["cuepeak"][[w,-1]],s=24,c=[PURPLE,color],zorder=5)
        ax.set_xscale("symlog",linthresh=1000)
        ax.set_xlim(0,14000);ax.set_ylim(0,1.68 if model=="wave_pinning" else 67)
        ax.set_xticks((0,500,1000,3000,14000),labels=("0","500","1000","3000","14000"))
        ax.set_ylabel("Peak activity inside physical cue support\n(native concentration)",fontsize=9)
        ax.set_xlabel("Time since cue onset (native; linear ≤1000, then log)",fontsize=8.5)
        ax.set_title(f"{'A' if col==0 else 'B'}   {LABELS[model]}",loc="left",fontsize=11,weight="bold")
        ax.legend(frameon=False,loc="lower left" if col==0 else "upper left",fontsize=8.5)
        text="Robust challenger retained" if col==0 else "Positive challenger, then decay"
        ax.annotate(f"{text}\ncue-end peak {c['cuepeak'][w]:.3f}",xy=(1000,c["cuepeak"][w]),
                    xytext=(.54,.56 if col==0 else .38),textcoords="axes fraction",fontsize=8.5,color=color,
                    arrowprops=dict(arrowstyle="-",color=PURPLE,lw=.9),
                    bbox=dict(facecolor="white",edgecolor="none",alpha=.9,pad=2))
        ax.grid(alpha=.13)
        ax=axes[1,col]
        ax.axvspan(36,39,color=PURPLE,alpha=.08)
        for state,idx,cc,ls in (("Initial",0,".3",":"),("Cue withdrawal",w,PURPLE,"--"),("After 13,000 off",-1,color,"-")):
            ax.plot(c["x"],c["u"][idx],color=cc,ls=ls,lw=1.5,label=state)
        ax.axhline(c["threshold"],color=".7",ls=":",lw=.7)
        ax.set_xlim(0,39);ax.set_ylim(0,1.68 if col==0 else 127)
        ax.set_xticks((0,10,20,30,39));ax.set_xlabel("Physical position along the L39 interval",fontsize=9)
        ax.set_ylabel("Native active concentration",fontsize=9)
        ax.set_title(f"{'C' if col==0 else 'D'}   Actual L39 spatial profiles",loc="left",fontsize=10.5,weight="bold")
        if col==0:
            box="Initial-side half-contrast width:\n11.92 → 6.65 → 6.74\nChallenged side: 0 → 4.61 → 5.24"
        else:
            box="Initial-side peak:\n117.39 → 100.05 → 117.39\nCue-region peak: 0 → 11.76 → 0"
        ax.text(.38,.97,box,transform=ax.transAxes,va="top",fontsize=8.5,
                bbox=dict(facecolor="white",edgecolor="none",alpha=.9,pad=2))
        ax.text(.987,.65,"cue region",transform=ax.transAxes,rotation=90,ha="right",fontsize=7.8,color=PURPLE)
        ax.legend(frameon=False,loc="upper left" if col==0 else "lower center",fontsize=7.8)
        ax.grid(alpha=.12)
    fig.text(.27,.477,"Broader incumbent; two domains retained after withdrawal",ha="center",fontsize=9.4,color=COLORS["wave_pinning"],weight="bold")
    fig.text(.75,.477,"Taller pole, lower substrate; weaker challenger, then decay",ha="center",fontsize=9.4,color=COLORS["goryachev"],weight="bold")
    fig.text(.5,.035,"Both native inputs transfer stimulus × inactive into active; finite responses include coupled reaction and transport.\n"
             "Native peaks: WP 0.1; Goryachev 6.4. These are not matched optical doses or biological clocks.\n"
             "Widths use each initial half-contrast threshold; zero threshold-defined extent does not mean absent activity. Finite patterns are not tracked site identities or a movement law.",
             ha="center",va="bottom",fontsize=8.3,linespacing=1.5)
    return fig


def organization_figure(values,excluded):
    fig,axes=plt.subplots(1,3,figsize=(11.2,5.0))
    fig.subplots_adjust(left=.075,right=.98,bottom=.27,top=.73,wspace=.30)
    fig.text(.5,.95,"Incumbent material increases through different spatial changes",ha="center",fontsize=16,weight="bold")
    fig.text(.5,.877,"All 97 eligible basal checkpoints · explicit normalization to each model's own reference · unchanged native measurements",ha="center",fontsize=9.3)
    specs=(("excess_amount_over_reference","E   Excess active amount"),("width_over_reference","F   Incumbent width"),("peak_over_reference","G   Peak concentration"))
    for ax,(key,title) in zip(axes,specs):
        ax.axvline(1,color=".8",lw=.7,ls=":");ax.axhline(1,color=".8",lw=.7,ls=":")
        for model in LABELS:
            rr=sorted((r for r in values if r["model"]==model),key=lambda r:r["size_native"])
            ax.plot([r["size_over_reference"] for r in rr],[r[key] for r in rr],color=COLORS[model],lw=1.4,marker="o",ms=2.6,label=LABELS[model])
        ax.set_title(title,loc="left",fontsize=10.5,weight="bold")
        ax.set_ylim(bottom=0);ax.set_xlim(0,2.7)
        ax.set_xlabel("Size / own reference size",fontsize=9)
        ax.set_ylabel("Measurement / own reference value",fontsize=9)
        ax.grid(alpha=.13)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="lower center",bbox_to_anchor=(.5,.12),ncol=5,frameon=False,fontsize=9)
    fig.text(.5,.018,"References: L19 for WP/Goryachev/deBelly; L20 for Holmes; resting L₀=1 for Spring.\n"
             "Amount is excess activity in the connected basal half-contrast core. Spring widths/amounts use its physical domain; deBelly retains nominal solver coordinates.\n"
             "Thirteen preparation-ineligible sizes are recorded separately, not assigned zero. Original full native-unit organization artwork remains linked.",
             ha="center",va="bottom",fontsize=8.1,linespacing=1.5)
    return fig


def save(fig,name):
    result=[]
    for ext in ("png","pdf","svg"):
        p=(OUT if ext=="png" else OUT/"supporting/vectors")/(name+"."+ext)
        fig.savefig(p,dpi=240,facecolor="white");result.append(p)
    return result


def write_docs(summary):
    original=source(REPO/"current manuscript figs/initial_polarity_length.png")
    assert sha(original)=="b94d05c1caabeb136a2cb744782ad56ac4e366bf2a00ca266caebd6af534f9ba"
    original_link="<"+os.path.relpath(original,OUT)+">"
    caption="""# Provisional Figure 2 candidate: finite challenge responses and incumbent organization

Two continued sheets for author review, separate from the rejected prototype. The organizing claim is that similar size-dependent reversal behavior can arise through distinct mechanisms predicting different challenge responses. The reversal trends remain in Figure1; they are not replotted here.

A–B: Actual peak activity in the physical right-hand cue footprint [L−3,L] at L19 (gray reference) and L39 (model colour). Each model uses the same requested physical support3, peak and1000-unit duration across those sizes: WP0.1, Goryachev6.4, followed by13000 cue-free native units. Both inputs add stimulus×inactive to active, with equal loss from inactive. Time is linear through1000 then logarithmic, with withdrawal marked at1000; the final coordinate14000 is13000 after withdrawal. Concentrations are native and are never independently peak-normalized. Goryachev39 develops a positive peak11.7598; it is10.581 at10 off and5.2796 at100 off, then approximately zero by1000 off. Its L19 reference instead grows to61.676 after an initial small withdrawal dip. WP's two sizes have nearly equal challenger height, including cue-end peaks1.5460/1.5434, despite different remaining incumbents.

C–D: Exact L39 spatial fields initially, at withdrawal, and after13000 cue-free units. Purple shading is the physical cue footprint, not activity extent. Dotted horizontal lines give each case's initial background-plus-half-range threshold. WP39 retains two connected above-threshold domains in every saved off-frame (2027 samples); their final initial-side/challenged-side widths are6.7382/5.2372. Goryachev39 never forms a challenger above its initial half-contrast58.695, but its positive cue-driven activity is explicitly visible in B and in the withdrawal profile. Its initial-side peak recovers from100.052 at withdrawal to117.390. Profile values are measured fields, not evidence inferred from a scalar-strength gate. The endpoint two-field comparisons in the source table remain distinct from any scalar criterion. Connected regions refer to spatial patterns, not tracked identities or protrusions.

E–G: All97 retained basal checkpoints, with excess amount, connected half-contrast width and peak concentration divided by the value at each model's own reference size. Sizes use L/Lref (WP/Goryachev/deBelly Lref19; Holmes20; Spring resting L₀ref1). Amount is the one-dimensional excess-concentration integral within the basal half-contrast core, not calibrated molecule number. Spring uses its actual physical domain for widths/amounts; deBelly retains nominal chemical-solver coordinates. Thirteen ineligible preparations remain in the linked table, without artificial zero values. At L19→39, WP's excess amount/width/peak folds are2.124/2.042/1.000, versus Goryachev2.152/1.119/1.903. All five curves and the original full native-unit organization artwork remain accessible.

The labels linking organization to response name conditional model predictions, not unique molecular identification or a causal allocation of reversal time. The native channels, doses, spatial discretization and clocks are not optical-actuator or reporter calibrations. The chemical grid has201 cells at both sizes; the sampled input peak is0.999388/0.997422 of the requested peak at L19/L39 (about0.2% mismatch), while its physical integral is nearly unchanged. This small cue-sampling mismatch does not bound basal or other discretization errors; Goryachev39's initial core has four sample centres and spatial convergence is not established by the retained checks. Thirteen thousand units of off-follow-up do not prove an asymptotic attractor. No velocity law or cell displacement is inferred. Small-loss recovery, analytic gain and evolved interventions remain supporting evidence, not the experimental headline.
"""
    fixes={"Figure1":"Figure 1","support3":"support 3","and1000":"and 1000","WP0.1":"WP 0.1","Goryachev6.4":"Goryachev 6.4","by13000":"by 13000","through1000":"through 1000","at1000":"at 1000","coordinate14000":"coordinate 14000","is13000":"is 13000","Goryachev39":"Goryachev L39","peak11.7598":"peak 11.7598","is10.581":"is 10.581","at10 off":"at 10 off","and5.2796":"and 5.2796","at100 off":"at 100 off","by1000 off":"by 1000 off","to61.676":"to 61.676","peaks1.5460":"peaks 1.5460","after13000":"after 13000","are6.7382":"are 6.7382","half-contrast58.695":"half-contrast 58.695","from100.052":"from 100.052","to117.390":"to 117.390","All97":"All 97","Lref19":"Lref 19","Holmes20":"Holmes 20","L₀ref1":"L₀ref 1","are2.124":"are 2.124","Goryachev2.152":"Goryachev 2.152","has201":"has 201","is0.999388":"is 0.999388","about0.2%":"about 0.2%"}
    for a,b in fixes.items():caption=caption.replace(a,b)
    (OUT/"CAPTION.md").write_text(caption)
    readme=f"""# Provisional Figure 2: finite responses first

This separate candidate follows the stated provisional finite-response-first emphasis. It is for author review, not author acceptance; the preceding prototype and all retained sources remain unchanged. No new integration was run.

![Finite native responses](figure2_finite_response.png)

![All-five incumbent organization](figure2_organization.png)

Read **A–B** for actual challenger activity during and after the fixed local cue: WP reaches almost the same challenger height at both sizes; Goryachev L39 develops a positive but smaller response that decays after release. **C–D** show what remains spatially: WP L39 retains two active domains, whereas Goryachev L39 restores its initial-side pattern. **E–G** restore the all-five starting-state comparison: active amount grows, chiefly through width in four models and peak height in Goryachev. The labels connect those organizations to conditional native-channel predictions; they do not uniquely identify a molecular mechanism.

[Two-sheet PDF](supporting/vectors/figure2_finite_candidate.pdf) · [Finite-response SVG](supporting/vectors/figure2_finite_response.svg) · [Organization SVG](supporting/vectors/figure2_organization.svg) · [Caption](CAPTION.md)

The local native curves expose a positive Goryachev39 challenger that was visually faint on the prior full-field colour scale. The selected endpoint fields distinguish WP39's co-retained domains from Goryachev39's incumbent return. L19 is a same-protocol within-model reference, not a matched cross-model optical stimulus. The restored all-five organization comparison makes peak-size differences visible alongside amount and width.

[Original unchanged all-five native-unit organization artwork]({original_link}) · [Exact97 basal source values](supporting/data/organization_all97.csv) · [Unavailable preparations](supporting/data/preparation_unavailable.csv)

[All four native histories and regional/edge reads](supporting/data/finite_native_readouts.csv) · [L39 native profile snapshots](supporting/data/native_profile_snapshots.csv) · [Protocols, native archives and separate full-field checks](supporting/data/selected_cases.csv) · [Validation/source hashes](supporting/metadata/validation.json)

Supporting evidence stays linked: [finite-window and amplitude audit of small-loss recovery](../optogenetic_revision/discussion_sheet/README.md), [analytic native gain and common recovery data](../supporting/data/all97_gain_recovery.csv), [endpoint-qualified evolved interventions](../supporting/data/selected_evolved_interventions.csv), and the [retained WP L15 short/long withdrawal illustration](../../../../revision_20261010/editorial/assets/wave_pinning_withdrawal.png). In that latter pair both cue-end domains remain inside support3; only the longer-pulse domain expands beyond the cue after withdrawal. It is a one-size duration illustration, not a size comparison.

The selected WP mobility intervention improves retained orientation at a fixed dose but both native and4× conditions fail the stronger full-polarity-strength gate; diffusion changes during cue plus recovery. None of these supporting controls provides an additive allocation of the shared reversal-size trend.
"""
    for a,b in {"Goryachev39":"Goryachev L39","WP39":"WP L39","Exact97":"Exact 97","support3":"support 3","and4×":"and 4×"}.items():readme=readme.replace(a,b)
    (OUT/"README.md").write_text(readme)
    rp=lambda p:str(p.relative_to(REPO))
    content=dict(title="Figure 2. Distinct finite responses to local polarity challenges",
                 role="Similar size-dependent reversal behavior can arise through distinct mechanisms, which predict different observable responses to polarity challenges.",
                 status="provisional author-review candidate; main emphasis remains revisable",
                 assets=[dict(path=rp(OUT/"figure2_finite_response.png"),caption="Finite native activation and withdrawal responses",page_size=[11.2,8.6],page_title="Figure 2. Finite challenge responses",minimum_printed_width_inches=10),
                         dict(path=rp(OUT/"figure2_organization.png"),caption="All-five incumbent organization",page_size=[11.2,5.0],page_title="Figure 2 continued. Incumbent organization",minimum_printed_width_inches=10)],
                 legend=caption.split("\n\n")[2:],
                 narrative=["Incumbent growth takes different spatial forms: most models widen; Goryachev grows taller and lowers available inactive substrate.",
                            "Under the native stimulus-times-inactive channel, Goryachev has weaker size-dependent challenger activity; its L39 response is positive and later decays. WP has nearly size-invariant challenger height but a larger remaining incumbent domain at L39.",
                            "After withdrawal WP39 retains two observed active domains while Goryachev39 restores the initial-side pattern over the finite follow-up. These are native response routes, not unique mechanistic classifications."],
                 measurements=["Activity peak/integrated increment inside the actual physical cue support, with native scale and a fixed within-model cue schedule.",
                               "Incumbent concentration, half-contrast extent/edge retreat and activity beyond the cue, followed after withdrawal.",
                               "Basal excess amount, connected width and peak concentration, preserving all97 native checkpoints and explicit reference normalization."],
                 limitations=["Native stimulus channel, dose and timescale are not experimentally calibrated; no optical actuator or reporter model is selected.",
                              "Finite 13000-unit off observations do not establish asymptotic stability, site ancestry, whole-field replacement from a scalar gate, or a movement law.",
                              "Goryachev39 has four native sample centres in its incumbent core; the retained evidence is not a spatial-convergence demonstration.",
                              "Observed response signatures and organization do not uniquely identify molecular mechanisms or causally partition reversal resistance."],
                 sources=[dict(label="Methods, exact source arrays and checks",path=rp(OUT/"supporting/metadata/validation.json")),
                          dict(label="Native protocols and endpoint comparisons",path=rp(OUT/"supporting/data/selected_cases.csv")),
                          dict(label="Original unchanged organization artwork",path=rp(original))],
                 supplementary=[dict(title="Gain and small-loss recovery: supporting evidence",assets=[],legend=["Calculated native input gain is separate from observed post-cue activity. The small-loss recovery contrast persists through stored finite trajectories but weakens with time; actual damage is only about0.22–0.285% of core activity, with no reporter calibration."],narrative=[],measurements=[],limitations=["Small-loss assays are not the finite optical challenge or a universal stability measure."],sources=[dict(label="Finite-time and amplitude audit",path=rp(F2/"optogenetic_revision/discussion_sheet/README.md")),dict(label="All97 gain/recovery values",path=rp(F2/"supporting/data/all97_gain_recovery.csv"))]),
                                dict(title="Selected native pulse and causal controls",assets=[],legend=["WP L15 T20/T60 starts from the same basal state and illustrates withdrawal-dependent domain survival. Evolved interventions have endpoint- and cue-specific limits."],narrative=[],measurements=[],limitations=["WP4× mobility passes orientation only at one dose; both conditions fail the stronger full-polarity-strength gate. Mobility is changed during cue and recovery."],sources=[dict(label="Retained WP short/long pulse artwork",path=rp(ROOT/"figures/length_mechanisms_20261008/revision_20261010/editorial/assets/wave_pinning_withdrawal.png")),dict(label="Endpoint-qualified controls",path=rp(F2/"supporting/data/selected_evolved_interventions.csv"))])])
    (OUT/"supporting/metadata/content.json").write_text(json.dumps(content,indent=2)+"\n")


def main():
    for p in (OUT,OUT/"supporting/data",OUT/"supporting/metadata",OUT/"supporting/vectors"):p.mkdir(parents=True,exist_ok=True)
    prior={str(p):sha(p) for p in F2.rglob("*") if p.is_file() and OUT not in p.parents}
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
    cases,traces,snapshots,selected=native_cases();org,excluded=organization()
    artifacts=[csvwrite(n,v) for n,v in (("finite_native_readouts.csv",traces),("native_profile_snapshots.csv",snapshots),
                ("selected_cases.csv",selected),("organization_all97.csv",org),("preparation_unavailable.csv",excluded))]
    first=finite_figure(cases);second=organization_figure(org,excluded)
    artifacts+=save(first,"figure2_finite_response")+save(second,"figure2_organization")
    combined=OUT/"supporting/vectors/figure2_finite_candidate.pdf"
    with PdfPages(combined) as pdf:pdf.savefig(first);pdf.savefig(second)
    artifacts.append(combined);plt.close(first);plt.close(second)
    write_docs(selected)
    artifacts += [OUT/"README.md",OUT/"CAPTION.md",OUT/"supporting/metadata/content.json"]
    source(Path(__file__))
    assert all(sha(Path(p))==h for p,h in prior.items())
    verification=dict(generated_at_utc=datetime.now(timezone.utc).isoformat(),passed=True,no_new_integrations=True,
                      status="provisional author-review candidate; independent review required before outline integration",
                      preserved_prior_figure2_files=len(prior),sources=SOURCES,
                      artifacts={str(p.relative_to(OUT)):sha(p) for p in artifacts},
                      selected_native_cases=selected,organization_checkpoints=97,unavailable_preparations=13,
                      finite_curve_rows=len(traces),native_profile_snapshot_rows=len(snapshots),
                      physical_cue_roi="[L-3,L]; peak over native centres; amount uses exact finite-volume overlap of physical cells with ROI",
                      halfcontrast="Frozen initial field minimum + half its range; boundary-connected contour crossings linearly interpolate native centres",
                      wp39_both_domains_in_all_saved_off_frames=2027,
                      inference="Field histories are inspected directly; scalar-strength gates do not prove complete field replacement.")
    (OUT/"supporting/metadata/validation.json").write_text(json.dumps(verification,indent=2)+"\n")
    print(json.dumps(dict(output=str(OUT),native_histories=4,organization_checkpoints=97,curve_rows=len(traces),prior_files_unchanged=len(prior)),indent=2))


if __name__=="__main__":main()
