"""Author-requested existing weak/strong medoid histories and native snapshots.

Separate discussion pages, without native integration or main-Figure4 changes.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

ROOT = Path(__file__).resolve().parents[4]
FIG = ROOT / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4"
INVENTORY = FIG / "regime_review_20261019"
OUT = FIG / "noise_strength_review"
DATA, META, VECTORS = (OUT / "supporting" / p for p in ("data", "metadata", "vectors"))
CASES = (("wave_pinning", 20., "Wave-pinning", "C=20"),
         ("holmes_model3", 60., "Holmes 3", "C=60"),
         ("spring_ring", 3., "Spring ring", "C₀=3"),
         ("goryachev", 20., "Goryachev", "C=20"))
COLORS = dict(initial="#777777", weak="#2679ad", strong="#d66b22", mean="#292929")


def rows(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8*1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def csvwrite(path, records):
    keys = list(dict.fromkeys(k for r in records for k in r))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(records)


def one(records, **query):
    found = [r for r in records if all((float(r[k]) == v if isinstance(v, float) else r[k] == v) for k, v in query.items())]
    assert len(found) == 1, (query, len(found))
    return found[0]


def load(row, selected_time):
    p = Path(row["profile_archive"])
    assert sha(p) == row["profile_sha256"]
    assert sha(row["archive"]) == row["scalar_sha256"]
    with np.load(p) as z:
        t, u = z["time"].copy(), z["primary"].copy()
        C = np.broadcast_to(z["circumference"], len(t)).copy()
        q = z["material_positions"].copy() if "material_positions" in z else z["x"].copy()/C[0]
    with np.load(row["archive"]) as z:
        st = z["time"].copy()
        ix = np.searchsorted(st, t)
        np.testing.assert_array_equal(st[ix], t)
        P = (z["P_real"]+1j*z["P_imag"])[ix]
        regions = z["current_half_contrast_regions"][ix]
        peaks = z["peaks_prominence_ge_005"][ix]
    j, = np.flatnonzero(t == selected_time)
    M0 = C*np.mean(u, axis=1)
    M1 = C*(u@np.exp(2j*np.pi*q))/len(q)
    error = float(np.max(abs(M1/M1[0]-P)))
    assert error < 1e-10
    return dict(row=row, time=t, primary=u, C=C, q=q, P=P, R=abs(P), M0=M0,
                M1=M1, regions=regions, peaks=peaks, selected_index=int(j), replay_error=error)


def number(v):
    return f"{v:.3f}" if v >= .001 else f"{v:.1e}"


def main():
    for p in (OUT, DATA, META, VECTORS):
        p.mkdir(parents=True, exist_ok=True)
    cases = rows(INVENTORY / "supporting/data/case_inventory.csv")
    medoids = rows(INVENTORY / "supporting/data/existing_medoid_inventory.csv")
    examples = rows(INVENTORY / "supporting/data/native_medoid_examples.csv")
    records, profile_rows, archived, artifacts = [], [], {}, []
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10,
                         "axes.spines.top":False, "axes.spines.right":False})
    combined = OUT / "noise_strength_native_comparison.pdf"
    with PdfPages(combined) as pdf:
        for model, size, name, size_label in CASES:
            if model == "wave_pinning":
                snapshot = 16000.
                selected = dict(first_positive_time=np.nan, last_positive_time=np.nan, observed_span=np.nan,
                                selection="Mid-window saved field; the unchanged strong C20 medoid has no multiregion sample")
            else:
                selected = one(examples, model=model, size=size, fraction=.01, kind="resolved_multiregion")
                snapshot = float(selected["profile_time"])
            source_rows = {
                "weak": one(medoids, model=model, size=size, fraction=.001),
                "strong": one(medoids, model=model, size=size, fraction=.01),
                "mean": one(cases, model=model, size=size, fraction=.01, condition="mean")}
            native = {role: load(r, snapshot) for role, r in source_rows.items()}
            np.testing.assert_array_equal(native["weak"]["primary"][0], native["strong"]["primary"][0])
            np.testing.assert_array_equal(native["weak"]["primary"][0], native["mean"]["primary"][0])
            np.testing.assert_array_equal(native["weak"]["q"], native["strong"]["q"])
            for role, z in native.items():
                j, r = z["selected_index"], z["row"]
                summary = dict(model=model, size=size, role=role, fraction=float(r["fraction"]),
                    slot=int(r["slot"]), seed=int(r["seed"]), selected_native_time=float(z["time"][j]),
                    scalar_cadence=float(r["scalar_cadence"]), profile_cadence=float(r["profile_cadence"]),
                    regions_at_snapshot=int(z["regions"][j]), prominent_peaks_at_snapshot=int(z["peaks"][j]),
                    R_at_snapshot=float(z["R"][j]), physical_C_at_snapshot=float(z["C"][j]),
                    physical_active_M0_at_snapshot=float(z["M0"][j]), physical_M1_real_at_snapshot=float(z["M1"][j].real),
                    physical_M1_imag_at_snapshot=float(z["M1"][j].imag), initial_physical_C=float(z["C"][0]),
                    existing_condition_medoid=role!="mean", scalar_archive=r["archive"], scalar_sha256=r["scalar_sha256"],
                    profile_archive=r["profile_archive"], profile_sha256=r["profile_sha256"],
                    selection=selected["selection"], strong_run_first_time=float(selected["first_positive_time"]),
                    strong_run_last_time=float(selected["last_positive_time"]), strong_run_observed_span=float(selected["observed_span"]),
                    physical_P_replay_maximum_error=z["replay_error"])
                records.append(summary)
                prefix=model+"_"+role
                for key in ("time", "primary", "q", "C", "P", "R", "regions", "peaks"):
                    archived[prefix+"_"+key] = z[key]
                for k, q in enumerate(z["q"]):
                    profile_rows.append(dict(model=model,size=size,role=role,time=snapshot,cell=k,
                        material_position=float(q),display_position=float((q+.5)%1-.5),physical_C=float(z["C"][j]),
                        physical_ds=float(z["C"][j]/len(z["q"])),native_active=float(z["primary"][j,k])))
            ini = native["weak"]
            for k, q in enumerate(ini["q"]):
                profile_rows.append(dict(model=model,size=size,role="initial",time=0.,cell=k,
                    material_position=float(q),display_position=float((q+.5)%1-.5),physical_C=float(ini["C"][0]),
                    physical_ds=float(ini["C"][0]/len(ini["q"])),native_active=float(ini["primary"][0,k])))
            vmax=max(float(z["primary"].max()) for z in native.values())
            fig=plt.figure(figsize=(11.4,9.6))
            gs=fig.add_gridspec(2,2,height_ratios=(1.2,.9),left=.085,right=.90,bottom=.20,top=.87,hspace=.39,wspace=.25)
            axs=[fig.add_subplot(gs[0,k]) for k in range(2)]
            for ax,role,label in zip(axs,("weak","strong"),("A","B")):
                z=native[role];pos=(z["q"]+.5)%1-.5;order=np.argsort(pos)
                mesh=ax.pcolormesh(pos[order],z["time"]/1000,z["primary"][:,order],vmin=0,vmax=vmax,
                                   cmap="viridis",shading="nearest",rasterized=True)
                ax.set_xlim(-.5,.5);ax.set_ylim(z["time"][-1]/1000,0)
                ax.set_title(f"{label}  {'Weaker' if role=='weak' else 'Stronger'} input: f={float(z['row']['fraction']):g}",loc="left",fontweight="bold",fontsize=11)
                ax.set_xlabel("Material boundary position / cycle")
                ax.set_ylabel("Time (10³ native units)")
                line=ax.axhline(snapshot/1000,color="white",ls="--",lw=1.0)
                line.set_path_effects([pe.Stroke(linewidth=2,foreground="black"),pe.Normal()])
                ax.text(.025,.975,f"seed {z['row']['seed']}\nsnapshot t={snapshot:g}",transform=ax.transAxes,va="top",fontsize=8.4,
                        color="white",bbox=dict(facecolor="black",alpha=.50,edgecolor="none",pad=3))
            cax=fig.add_axes([.92,.535,.016,.328])
            cb=fig.colorbar(mesh,cax=cax);cb.set_label("Native active concentration",fontsize=10)
            ax=fig.add_subplot(gs[1,:]);q=ini["q"];pos=(q+.5)%1-.5;order=np.argsort(pos)
            curves=[("initial",ini["primary"][0],"Initial, t=0", "-",1.4)]
            curves += [(role,z["primary"][z["selected_index"]],
                        f"{'Weak' if role=='weak' else 'Strong' if role=='strong' else 'Uniform mean f=.01'}: R={number(z['R'][z['selected_index']])}",
                        "--" if role=="mean" else "-",1.5 if role!="mean" else 1.15)
                       for role,z in native.items()]
            for role, u, label, style, lw in curves:
                ax.plot(pos[order],u[order],style,color=COLORS[role],label=label,lw=lw,zorder=3 if role=="strong" else 2)
            ax.set_xlim(-.5,.5);ax.set_ylim(0,vmax*1.06)
            ax.set_xlabel("Material boundary position / cycle");ax.set_ylabel("Native active concentration")
            ax.set_title(f"C  Actual profiles at t={snapshot:g}; initial reference at t=0",loc="left",fontweight="bold",fontsize=11)
            ax.legend(loc="upper center",bbox_to_anchor=(.5,-.24),ncol=2,frameon=False,fontsize=9)
            ns=[int(native[k]["regions"][native[k]["selected_index"]]) for k in ("weak","strong","mean")]
            ax.text(.015,.96,f"Half-contrast regions at snapshot: weak {ns[0]}, strong {ns[1]}, mean {ns[2]}",
                    transform=ax.transAxes,va="top",fontsize=8.5,
                    bbox=dict(facecolor="white",edgecolor="none",alpha=.85,pad=2))
            if model=="goryachev":
                zoom=inset_axes(ax,width="37%",height="58%",loc="upper right",borderpad=1.25)
                for role,u,label,style,lw in curves:
                    zoom.plot(pos[order],u[order],style,color=COLORS[role],lw=1.0)
                detailmax=max(float(native[k]["primary"][native[k]["selected_index"]].max()) for k in ("strong","mean"))
                zoom.set_ylim(0,np.ceil(detailmax*1.1*2)/2);zoom.set_xlim(-.5,.5)
                zoom.set_title("Native-scale zoom (same curves)",fontsize=8)
                zoom.tick_params(labelsize=7)
            fig.suptitle(f"Existing noise regimes · {name}, {size_label}",fontsize=16,y=.975)
            fig.text(.5,.932,"Independent weak/strong condition medoids · identical starting basal · shared zero-based native scale",ha="center",fontsize=9.7)
            fig.text(.5,.907,"Positive input: increasing f raises mean and standard deviation together; this is not a zero-mean variance sweep.",ha="center",fontsize=8.6)
            if model=="wave_pinning":
                selection="Snapshot: middle saved time; this strong-input medoid has no multiple-region sample."
            else:
                selection=f"Snapshot: saved field in the strong medoid's qualified multiregion sample run t={float(selected['first_positive_time']):g}–{float(selected['last_positive_time']):g}."
            fig.text(.085,.070,selection,fontsize=8.8)
            fig.text(.085,.045,"The weak snapshot uses the identical native time; these are separate realizations, not paired forcing paths.",fontsize=8.8)
            if model=="goryachev":
                fig.text(.085,.020,"Mean-input control loses polarity: the stronger noisy case is a forced-pattern regime, not intact-pole wandering.",fontsize=8.8,color="#9a3824")
            else:
                fig.text(.085,.020,"Multiple threshold regions describe local organization; counts do not establish competition or site lineage. Native units are uncalibrated.",fontsize=8.6)
            stem="noise_strength_"+model
            for suffix in ("png","pdf","svg"):
                path=(OUT if suffix=="png" else VECTORS)/(stem+"."+suffix)
                fig.savefig(path,dpi=240)
                artifacts.append(path)
            pdf.savefig(fig)
            plt.close(fig)
            print(str(OUT/(stem+".png")),flush=True)
    artifacts.append(combined)
    csvwrite(DATA/"snapshot_sources_and_selection.csv",records)
    csvwrite(DATA/"native_snapshot_profiles.csv",profile_rows)
    np.savez_compressed(DATA/"native_comparison_arrays.npz",**archived)
    inputs=[INVENTORY/"supporting/data/case_inventory.csv",INVENTORY/"supporting/data/existing_medoid_inventory.csv",
            INVENTORY/"supporting/data/native_medoid_examples.csv"]
    main_art=json.loads((FIG/"supporting/metadata/figure4_with_spring_render_validation.json").read_text())["output_hashes"]
    main_art=[r for r in main_art if Path(r["path"]).stem in ("figure4_with_spring_dynamics","figure4_with_spring_geometry_patterns")]
    for r in main_art:
        assert sha(r["path"])==r["sha256"]
    audit=dict(generated_at=datetime.now(timezone.utc).isoformat(),source=str(Path(__file__).relative_to(ROOT)),
        source_sha256=sha(__file__),native_integrations=0,main_artifacts_unchanged=main_art,
        comparison_models=4,existing_noisy_condition_medoids=8,mean_controls=4,
        physical_P_replay_maximum_error=max(r["physical_P_replay_maximum_error"] for r in records),
        snapshot_profiles=len(profile_rows),same_native_snapshot_times_in_weak_and_strong=True,
        representative_selection="All eight noisy histories are unchanged condition medoids",
        snapshot_selection="Existing medoid's longest qualified positive sample run midpoint-nearest saved field; WP mid-window because no multiregion sample",
        input_hashes=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for p in inputs],
        artifacts=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for p in artifacts],
        compact_data=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for p in sorted(DATA.iterdir())])
    (META/"validation.json").write_text(json.dumps(audit,indent=2)+"\n")
    print(str(combined),flush=True)


if __name__=="__main__":
    main()
