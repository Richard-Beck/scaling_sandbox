"""Compact external-cue comparison: exact RHS on twelve existing full histories."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import sys
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import revision_arc_narrow_collect_20261019 as api

ROOT = Path(__file__).resolve().parents[2]
TOPIC = ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/narrow_arc_experiment"
OUT = TOPIC/"cue_on_core"
RUN = ROOT/"data/runs/revision_arc_narrow_grid_20261019"
NAMES = dict(wave_pinning="Wave-pinning", goryachev="Goryachev", holmes_model3="Holmes", spring="Spring")


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    for leaf in ("data", "metadata", "vectors"):
        (OUT/"supporting"/leaf).mkdir(parents=True, exist_ok=True)
    if "--render-only" in sys.argv:
        data=pd.read_csv(OUT/"supporting/data/cue_on_core_native_rates.csv",float_precision="round_trip")
        selected=pd.read_csv(OUT/"supporting/data/cue_on_core_conditions.csv",float_precision="round_trip")
        metadata_path=OUT/"supporting/metadata/verification.json"
        report=json.loads(metadata_path.read_text())
        report.setdefault("native_derivation_renderer_sha256",report["renderer_sha256"])
        draw(data,selected)
        report.update(presentation_updated_at=datetime.now(timezone.utc).isoformat(),renderer_sha256=sha(__file__),
                      presentation_update="Render only from unchanged exact source CSVs; no RHS reevaluation.")
        report["artifact_sha256"]={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob("*") if p.is_file() and p.name!="verification.json"}
        metadata_path.write_text(json.dumps(report,indent=2)+"\n")
        print(OUT/"cue_on_external_size_comparison.png")
        return
    summary_path = TOPIC/"supporting/data/narrow_arc_summary.csv"
    table = pd.read_csv(summary_path)
    assert len(table) == 645 and (table.status == "complete").all()
    selected = table[(table.cue_duration == 1000)&(((table.model!="spring")&(table.physical_distance_at_onset==11))|((table.model=="spring")&(table.physical_distance_at_onset==1.05)))].copy()
    assert len(selected) == 12 and (selected.execution == "new_native_branch").all()
    assert (selected.footprint_relative_to_core_10 == "outside").all()
    assert (selected.full_field_tail_status == "passes_finite_tail_gate").all()
    rows = []; sources = {str(summary_path):sha(summary_path)}; max_identity_error = 0.
    manifest = json.loads((RUN/"manifest.json").read_text())
    lookup = {r["condition_id"]:r for r in manifest["cases"]}
    for r in selected.to_dict("records"):
        cid = int(r["condition_id"]); case = lookup[cid]; spring = r["model"] == "spring"; n = int(r["cells"])
        p = Path(r["actual_cue_path"]); assert sha(p) == r["actual_cue_sha256"]; sources[str(p)] = sha(p)
        cue_input = RUN/r["cue_input"]; assert sha(cue_input) == r["cue_input_sha256"]; sources[str(cue_input)] = sha(cue_input)
        with np.load(cue_input) as z: signal = z["initial_cue"].copy()
        with np.load(p) as z:
            states, time = z["states"].copy(), z["time"].copy()
            q = z["material_positions"].copy() if spring else z["x"].copy()/r["circumference"]
        model = api.spring_api.Closure(r["resting_circumference"], n) if spring else api.base.ring.native_model(r["model"])
        grid = None if spring else api.base.ring.RingGrid(r["circumference"], n)
        w = np.exp(2j*np.pi*q)
        initial = complex(r["M1_initial_real"],r["M1_initial_imag"])
        for j, (t, state) in enumerate(zip(time, states)):
            if spring:
                u, C = state[:n], state[-1]
                current = api.spring_api.physical_raised_cosine(model,state,r["cue_centre_fraction"],r["requested_peak"],r["full_physical_support"])[0]
                derivative = model.rhs(state,current)
                dM1 = C*np.mean(derivative[:n]*w)+derivative[-1]*np.mean(u*w)
            else:
                u, C = state[0],r["circumference"]
                derivative = api.base.ring.native_rhs(model,grid,state,signal)
                dM1 = C*np.mean(derivative[0]*w)
            M1 = C*np.mean(u*w); P = M1/initial; R = abs(P); heading = np.degrees(np.angle(P))
            target = r["target_angle_in_P_coordinates_deg"]
            closure = 1-abs((heading-target+180)%360-180)/abs((-target+180)%360-180) if R>=.1 else np.nan
            omega = float(np.degrees(np.imag(np.conj(M1)*dM1)/abs(M1)**2)) if R>=.1 else np.nan
            rows.append(dict(condition_id=cid,model=r["model"],size_parameter=r["size_parameter"],initial_arc=r["physical_distance_at_onset"],
                             actual_target_deg=target,frame=j,cue_time=float(t),P_real=P.real,P_imag=P.imag,R=R,heading_deg=heading if R>=.1 else np.nan,
                             target_error_removed=closure,native_RHS_angular_rate_deg_per_time=omega,native_RHS_target_relative_rate=omega/abs(target),
                             M0_ratio=C*np.mean(u)/r["M0_initial"],C=float(C),prescribed_physical_arc=r["physical_distance_at_onset"]*C/r["circumference"],
                             M1_rhs_real=float(dM1.real),M1_rhs_imag=float(dM1.imag)))
        last = rows[-1]
        error = max(abs(last["P_real"]-r["P_withdrawal_real"]),abs(last["P_imag"]-r["P_withdrawal_imag"]))
        max_identity_error=max(max_identity_error,error)
        assert error < 1e-12
    data = pd.DataFrame(rows)
    data.to_csv(OUT/"supporting/data/cue_on_core_native_rates.csv",index=False)
    selected.to_csv(OUT/"supporting/data/cue_on_core_conditions.csv",index=False)
    draw(data, selected)
    stamp=datetime.now(timezone.utc).isoformat()
    (OUT/"README.md").write_text(f"""# Common external initial arc: compact cue-on comparison

Generated {stamp}. These twelve selected histories are part of the complete645-condition experiment; no integration is performed here.

[PNG](cue_on_external_size_comparison.png) · [PDF](supporting/vectors/cue_on_external_size_comparison.pdf) · [SVG](supporting/vectors/cue_on_external_size_comparison.svg) · [Every actual saved-state rate](supporting/data/cue_on_core_native_rates.csv) · [Selected conditions and all source identities](supporting/data/cue_on_core_conditions.csv).

Selection is declared by geometry and coverage: chemical d0=11 and Spring d0=1.05 are the nearest common tested distances whose entire physical supports clear the initial10%-contrast region in every model and all three tested sizes. All twelve are fresh full-stateT1000 branches, all have R≥0.400 during the cue, and all pass the separately measured finite field-tail gate. They are not chosen by a minimum-pulse or retained-success threshold. Sizes are chemicals30/35/40 and Spring rest2.5/2.75/3. Cue support, coefficient peak, native kinetics/closure, solver and basals remain unchanged. Initial target angles vary across sizes and are printed explicitly. Holding d0 fixed changes both angular geometry and native basal architecture.

The first column plots branch-safe fraction of initial target-angle error removed, not a binary success or speed criterion. The second shows R=|M1(t)/M1(0)|. The third shows actual native-RHS angular rate at every saved full state, in degrees per native time; it is not a sampled secant or derivative with respect to displayed logarithmic time. WeakR<0.1 would mask phase and rate, but no selected cue-on curve has that gap. Physical M1 includes instantaneousC(t), including Spring's Cdot term. P axes are rotated once by the initialM1 phase; rates concern the global activity resultant rather than local-peak travel or cell movement. Onset rates depend on finite-width input projection and initial moment magnitude, so rate comparisons alone do not isolate an incumbent-remodeling mechanism. Native clocks are not calibrated across models. The full supplement retains every other placement, overlap regime, duration and finite endpoint.
""")
    artifacts=[p for p in OUT.rglob("*") if p.is_file() and p.name!="verification.json"]
    result=dict(generated_at=stamp,passed=True,selected_conditions=12,saved_full_states=len(data),minimum_R=float(data.R.min()),
                maximum_withdrawal_identity_error=max_identity_error,source_sha256=sources,renderer_sha256=sha(__file__),
                selection="Nearest common tested fully10%-external initial offset in every model/size; all full native states, no weak cue-on vector gap.",
                new_integrations=0,broader_RHS_scan=False,artifact_sha256={str(p.relative_to(OUT)):sha(p) for p in artifacts})
    (OUT/"supporting/metadata/verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:result[k] for k in ("generated_at","selected_conditions","saved_full_states","minimum_R","maximum_withdrawal_identity_error")},indent=2))


def draw(data, selected):
    plt.rcParams.update({"font.size":10,"pdf.fonttype":42,"svg.fonttype":"none"})
    colors=["#c98a20","#3c8894","#564d86"]
    fig,axes=plt.subplots(4,3,figsize=(14,13))
    fig.subplots_adjust(left=.085,right=.98,top=.86,bottom=.13,hspace=.54,wspace=.27)
    fig.text(.04,.965,"Cue-on reorientation at a common external physical offset",fontsize=19,weight="bold")
    fig.text(.04,.934,"T=1000 · chemical d₀=11; Spring d₀=1.05 · full supports outside every initial 10%-contrast region",fontsize=11)
    for col,label in enumerate(["Angular target error removed","Vector magnitude R","Native instantaneous angular rate"]):
        axes[0,col].text(0,1.30,label,transform=axes[0,col].transAxes,fontsize=12,weight="bold")
    for row,(model,name) in enumerate(NAMES.items()):
        part=data[data.model==model];cases=selected[selected.model==model].sort_values("size_parameter")
        axes[row,0].text(-.18,1.18,name,transform=axes[row,0].transAxes,fontsize=12,weight="bold")
        handles=[]
        for color,r in zip(colors,cases.to_dict("records")):
            values=part[part.condition_id==r["condition_id"]]
            label=("C₀=" if model=="spring" else "C=")+f"{r['size_parameter']:g}; θ*={r['target_angle_in_P_coordinates_deg']:.2f}°"
            handles.append(Line2D([0],[0],color=color,lw=1.8,label=label))
            for col,key in enumerate(["target_error_removed","R","native_RHS_angular_rate_deg_per_time"]):
                axes[row,col].plot(values.cue_time,values[key],color=color,lw=1.5)
        axes[row,0].legend(handles=handles,fontsize=7.8,frameon=False,loc="upper left",handlelength=1.5)
        axes[row,0].set_ylim(-.03,1.04);axes[row,0].axhline(1,color="#aab7bd",ls=":",lw=.7)
        axes[row,1].set_ylim(0,max(1.15,float(part.R.max())*1.04));axes[row,1].axhline(.1,color="#999",ls="--",lw=.7)
        axes[row,2].set_yscale("symlog",linthresh=.001)
        axes[row,2].set_ylabel("degrees / native time",fontsize=9)
        axes[row,2].axhline(0,color="#aaa",lw=.6)
        for ax in axes[row,:]:
            ax.set_xscale("symlog",linthresh=.1);ax.set_xlim(0,1000);ax.set_xticks([0,.1,1,10,100,1000],["0",".1","1","10","100","1000"])
            ax.spines[["top","right"]].set_visible(False);ax.tick_params(labelsize=8.5);ax.minorticks_off()
            if row==3:ax.set_xlabel("Time since cue onset (native units)",fontsize=10)
    fig.text(.04,.083,"Twelve actual full-state histories. Error removal uses each case’s actual target; raw headings and physical geometry remain in the supporting table.",fontsize=9.4)
    fig.text(.04,.062,"Time is linear to 0.1 then logarithmic. Rate axes are linear within ±0.001 then logarithmic in magnitude; this is not a calibrated cross-model clock.",fontsize=9.3)
    fig.text(.04,.041,"Global angular rate includes finite-width input projection and initial moment magnitude. Native profiles and full outcome/coverage maps are linked in the supplement.",fontsize=9.1)
    for ext in ("png","pdf","svg"):
        path=OUT/"cue_on_external_size_comparison.png" if ext=="png" else OUT/"supporting/vectors"/("cue_on_external_size_comparison."+ext)
        fig.savefig(path,dpi=180,facecolor="white")
    plt.close(fig)


if __name__=="__main__":main()
