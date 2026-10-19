"""Complete-scientific-pilot review; no integration or protocol change.

The accepted 229-output snapshot and its renderer remain unchanged. This
version imports that exact renderer and adds explicit scientific coverage,
an adjacent raw/qualified HTML table and a separate native-profile context page.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import html
import importlib.util
import json
import os

for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = "1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

OLD = Path(__file__).with_name("revision_figure4_broad_strength_review_20261019.py")
assert hashlib.sha256(OLD.read_bytes()).hexdigest() == "61dae8ab6b51833cc3b4ac93c45b1d9aeab6e983076171a60536c7347ab1ca33"
spec = importlib.util.spec_from_file_location("accepted_partial_renderer", OLD)
base = importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
OUT = base.OUT.with_name("review_complete_scientific")


def transition_profiles(run, out):
    cache = json.loads((run / "results/partial/case_cache.json").read_text())
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 10.8))
    fig.subplots_adjust(left=.075, right=.97, bottom=.15, top=.85, hspace=.52, wspace=.29)
    fig.suptitle("Native profile context across the sampled change in organization", x=.075, ha="left", y=.974, fontsize=18)
    fig.text(.075, .939, "Largest original size per model; common f=0.01,0.03,0.1. First saved seed per condition, at fixed native midwindow time.", fontsize=10.8)
    fig.text(.075, .914, "Profiles remain unshifted and keep actual native amplitudes. Independent paths across strengths; these snapshots are not occupancy estimates.", fontsize=10.4)
    colors = ["#2676a6", "#bd7b23", "#b04d63"]
    selected, cells = [], []
    for ax, (model, name) in zip(axes.flat, base.MODELS.items()):
        size = 3. if model == "spring_ring" else 60.
        curves = []
        for f, color in zip((.01, .03, .1), colors):
            cases = [c for c in cache.values() if c["descriptors"][0]["model"] == model and c["descriptors"][0]["size"] == size and c["descriptors"][0]["fraction"] == f and c["descriptors"][0]["kind"] == "noise"]
            c = min(cases, key=lambda c:c["descriptors"][0]["seed"])
            meta = c["descriptors"][0]
            with np.load(c["profile"]) as z:
                t, u = z["time"], z["primary"]
                C = np.broadcast_to(z["circumference"], len(t))
                x = z["material_positions"] if model == "spring_ring" else z["x"]
                extent = 1. if model == "spring_ring" else size
            k, = np.flatnonzero(t == t[-1]/2)
            with np.load(c["trace"]) as z:
                j, = np.flatnonzero(z["time"] == t[k])
                R = float(abs(z["P_real"][j]+1j*z["P_imag"][j]))
                A = float(z["M0"][j]/z["M0"][0]); Q = float(z["q1"][j]/z["q1"][0])
                regions = int(z["current_half_contrast_regions"][j])
                contrast = float(z["contrast"][j]/z["contrast"][0])
            xx = np.r_[x[-1]-extent, x, x[0]+extent]
            if not curves:
                basal = u[0].copy()
                ax.plot(xx,np.r_[basal[-1],basal,basal[0]],color="#8c969d",lw=1.5,label="Exact basal",zorder=1)
            else:
                np.testing.assert_array_equal(basal,u[0])
            ax.plot(xx,np.r_[u[k,-1],u[k],u[k,0]],color=color,lw=1.5,
                    label=f"f={f:g}: {regions} regions, R={R:.2f}",zorder=2)
            curves.append((x,u[k].copy(),color))
            selected.append(dict(model=model,size=size,fraction=f,seed=meta["seed"],case_id=meta["case_id"],time=float(t[k]),
                selection="Smallest saved seed number in condition; exact native midpoint, no field conditioning or alignment",
                native_circumference=float(C[k]),coordinate="material_fraction" if model=="spring_ring" else "physical_arc",
                regions=regions,R=R,A=A,Q=Q,contrast_ratio=contrast,profile=c["profile"],profile_sha256=c["profile_sha256"],trace=c["trace"],trace_sha256=c["trace_sha256"]))
            cells.extend(dict(model=model,fraction=f,time=float(t[k]),coordinate=float(a),primary=float(b),initial=float(d),circumference=float(C[k])) for a,b,d in zip(x,u[k],u[0]))
        ax.set(xlim=(0,extent),ylim=(0,None),xlabel="Material position q" if model=="spring_ring" else "Native physical arc coordinate",ylabel="Native active concentration")
        ax.set_title(f"{name} · {'C₀=3, t=1600' if model=='spring_ring' else 'C=60, t=16000'}",loc="left",fontsize=11.5)
        ax.spines[["top","right"]].set_visible(False)
        ax.legend(frameon=False,fontsize=8.7,loc="upper right")
        if model == "goryachev":
            inset = ax.inset_axes([.37,.31,.59,.3])
            for x,u,color in curves[1:]:inset.plot(x,u,color=color,lw=1.)
            inset.set_xlim(0,extent);inset.tick_params(labelsize=7)
            inset.set_title("f=0.03/0.1 native-scale enlargement",fontsize=7.5)
    fig.text(.075,.092,"The chosen strengths bracket observed occupancy changes on the finite grid; they are not phase boundaries. Actual sampled region counts are printed in each legend.",fontsize=9.6)
    fig.text(.075,.067,"Goryachev’s large-cell high-input profiles have little contrast relative to basal; half-current-contrast regions may still be counted. Native-scale enlargement preserves that distinction.",fontsize=9.4)
    fig.text(.075,.042,"Spring uses fixed material coordinates and records the actual C(t) for physical M₀/M₁. The original 0.01 level has12 seeds; 0.03/0.1 currently have2. No peak normalization.",fontsize=9.4)
    pd.DataFrame(selected).to_csv(out/"supporting/data/transition_profile_selection.csv",index=False)
    pd.DataFrame(cells).to_csv(out/"supporting/data/transition_profile_cells.csv",index=False)
    return fig


def report(df,status,out):
    esc=html.escape
    lines=["<!doctype html><html><head><meta charset='utf-8'><title>Complete scientific broad-strength pilot</title><style>body{font:16px system-ui;max-width:1200px;margin:2rem auto;padding:0 1rem;color:#263640}table{border-collapse:collapse;width:100%;font-size:14px;margin:1rem 0 2rem}th,td{border-bottom:1px solid #d8dfe4;text-align:right;padding:.4rem}td:first-child,th:first-child{text-align:left}.note{background:#fff1d2;padding:1rem;line-height:1.5}</style></head><body>",
        "<h1>Complete scientific pilot: broader input strengths</h1>",
        f"<p>Report generated {base.now()} · native status observed {esc(status['data_observed_at'])}: {esc(status['producer'])}. Scientific cases {status['scientific_cases_available']}/576; numerical companions {status['numerical_companions_complete']}/4 complete.</p>",
        "<p class='note'>New levels have two independent seeds; the original0.001/0.01 retain12. This is exploratory native evidence. Numerical qualification/replication remains subject to root and independent review. This timestamped export has no timed refresh; its saved status does not claim the producer remains in that state.</p>",
        "<p><a href='broad_strength_review.pdf'>Two-page review</a> · <a href='broad_strength_transition_profiles.png'>Native change-context profiles</a> · <a href='README.md'>Interpretation, methods and limits</a> · <a href='supporting/data/condition_window_snapshot.csv'>All collected windows and seed ranges</a></p>",
        "<p><b>Raw</b>: more than one connected arc above current minimum plus half current contrast. <b>Qualified</b>: raw multiple regions and contrast at least10% of that case’s initial basal contrast. <b>Low contrast</b>: below that10% level. All fractions retain all saved scalar times in their denominator; these diagnostics may overlap.</p>"]
    def val(r,k,pct=False):
        v=r[k];return "—" if not np.isfinite(v) else f"{100*v:.2f}%" if pct else f"{v:.3f}"
    for model,name in base.MODELS.items():
        for size in ((1.,2.,3.) if model=="spring_ring" else (20.,40.,60.)):
            lines.append(f"<h2>{name}, {'C₀' if model=='spring_ring' else 'C'}={size:g}</h2><table><tr><th>f</th><th>Noisy n</th><th>RAW multiple regions</th><th>Contrast-qualified multiple regions</th><th>Low contrast</th><th>Mean R</th><th>Mean A</th><th>Mean Q</th><th>Mean-input R</th></tr>")
            for f in base.GRID:
                mask=(df.model==model)&(df['size']==size)&(df.fraction==f)&(df.window=='full')
                r=df[mask&(df.kind=='noise')].iloc[0];c=df[mask&(df.kind=='mean')].iloc[0]
                lines.append(f"<tr><td>{f:g}</td><td>{int(r.available_cases)}</td><td>{val(r,'multiregion_fraction',True)}</td><td>{val(r,'qualified_multiregion_fraction',True)}</td><td>{val(r,'low_contrast_fraction',True)}</td><td>{val(r,'mean_R')}</td><td>{val(r,'mean_A')}</td><td>{val(r,'mean_Q')}</td><td>{val(c,'mean_R')}</td></tr>")
            lines.append("</table>")
    lines.append("<p>A=M₀/M₀(0), Q=q₁/q₁(0), R=AQ at each frame; mean R is mean of paired products. Full chemical windows0–32000/cadence5; Spring0–3200/cadence0.5. Native units/laws are not calibrated across models. Counts do not establish competition, event classes or phase boundaries.</p></body></html>")
    (out/"index.html").write_text("\n".join(lines))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,default=base.RUN);p.add_argument('--out',type=Path,default=OUT);a=p.parse_args()
    for sub in ('data','metadata','vectors'):(a.out/'supporting'/sub).mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42})
    df,status,csvsha=base.snapshot(a.run,a.out)
    assert status['scientific_cases_available']==576
    fig=base.main_sheet(df,status)
    for text in fig.texts:
        if text.get_text().startswith('Native status observed'):
            text.set_text(f"Scientific coverage576/576 · numerical companions{status['numerical_companions_complete']}/4 complete · native observed {status['data_observed_at'].replace('T',' ')[:19]} UTC")
        if text.get_text().startswith('Filled circles:'):
            text.set_text('All scientific conditions complete: circles mark n=2 new levels; squares mark n=12 original levels. Deterministic mean controls have one path.')
    base.export(fig,a.out,'broad_strength_exploratory')
    profiles,selection=base.profile_sheet(a.run,a.out);base.export(profiles,a.out,'broad_strength_native_profiles')
    with PdfPages(a.out/'broad_strength_review.pdf') as pdf:pdf.savefig(fig);pdf.savefig(profiles)
    context=transition_profiles(a.run,a.out);base.export(context,a.out,'broad_strength_transition_profiles')
    report(df,status,a.out);plt.close('all')
    df[(df.window=='full')&(df.kind=='noise')].to_csv(a.out/'supporting/data/all_strengths_full_window.csv',index=False)
    result=dict(generated_at=base.now(),source=str(Path(__file__).resolve()),source_sha256=base.sha(__file__),accepted_renderer=str(OLD),accepted_renderer_sha256=base.sha(OLD),
        condition_snapshot_sha256=csvsha,scientific_cases=576,numerical_companions_complete=status['numerical_companions_complete'],qualification='pending_root_and_independent_review',
        native_observed_at=status['data_observed_at'],native_producer_last_observed=status['producer'],main_figure_changed=False)
    result['outputs']={str(q.relative_to(a.out)):base.sha(q) for q in sorted(a.out.rglob('*')) if q.is_file() and q.name!='RENDER_RECEIPT.json'}
    base.savejson(a.out/'supporting/metadata/RENDER_RECEIPT.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
