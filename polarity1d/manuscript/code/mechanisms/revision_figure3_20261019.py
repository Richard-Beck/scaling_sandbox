"""Reanalyse retained native rings for the 19 October geometry revision.

No integration or model changes. Main figure: coordinate geometry, controlled
placement contrast, and the distinct information in local activity and its first
angular moment. Descriptor transfer is an explicitly descriptive supplement.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os

for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, Patch
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parents[1]
RUN = ROOT / "data/runs/chemical_ring_functional_v2_20261017"
MOM = ROOT / "data/runs/vision_ring_moments_20261017"
OUT = ROOT / "figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3"
MODELS = ("wave_pinning", "goryachev", "holmes_model3")
NAMES = dict(zip(MODELS, ("Wave-pinning", "Goryachev", "Holmes 3")))
CLASSES = {
    "one old-core-associated peak; incumbent identity unresolved": 0,
    "one peak associated with cue patch": 1,
    "distinct peaks associated with old core and cue patch": 2,
    "one connected region spanning old and cue (merged or bridged)": 3,
    "one shifted or merged peak; original identity unresolved": 4,
}
LABELS = ("Old-associated single", "Cue-associated single", "Distinct old + cue",
          "One spanning region", "Shifted / unresolved")
COLORS = ("#4477AA", "#CC7339", "#8764A6", "#369B90", "#999345")
SOURCES = {}
ARTIFACTS = []


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def source(path):
    path = Path(path).resolve()
    if str(path) not in SOURCES:
        SOURCES[str(path)] = sha(path)
    return path


def dest(name, group="data"):
    parent = OUT if name.endswith(".png") else OUT / "supporting" / group
    parent.mkdir(parents=True, exist_ok=True)
    path = parent / name
    ARTIFACTS.append(path)
    return path


def export(stem, fig):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(dest(stem + "." + ext, "vectors"), dpi=220,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def table(name, frame):
    frame.to_csv(dest(name), index=False)


def jsonwrite(name, value):
    dest(name, "metadata").write_text(json.dumps(value, indent=2) + "\n")


def load_phase(path, offset):
    with np.load(source(path)) as z:
        primary = z["states"][:, 0].copy() if "states" in z else z["primary"].copy()
        return z["time"].copy() + offset, z["x"].copy(), primary


def history(row):
    tc, x, pc = load_phase(row.cue_path, 0)
    tr, xr, pr = load_phase(row.recovery_path, 1000)
    np.testing.assert_array_equal(x, xr)
    np.testing.assert_array_equal(pc[-1], pr[0])
    t = np.r_[tc, tr[1:]]
    p = np.concatenate((pc, pr[1:]))
    assert t[0] == 0 and t[-1] == 14000 and np.all(np.diff(t) > 0)
    C = row.circumference
    theta = 2 * np.pi * (x - x[0]) / C
    dx = C / len(x)
    m1 = dx * (p @ np.exp(1j * theta))
    m2 = dx * (p @ np.exp(2j * theta))
    m0 = dx * p.sum(axis=1)
    vector = m1 / m1[0]
    return dict(t=t, x=x-x[0], p=p, P=vector, R=abs(vector), q1=abs(m1)/m0,
                q2=abs(m2)/m0, M0=m0)


def descriptor_transfer(frame):
    """Fixed k NN, leave one entire circumference out, ties included.

    No hyperparameter tuning. Training-only feature SD sets scales. Majority
    ties break toward the more frequent training class, then ascending class.
    The 3/15-neighbour repeats assess rule sensitivity, not uncertainty.
    """
    specs = [("Model majority", []), ("Physical separation", ["short_route"]),
             ("Angular separation", ["angular_separation_deg"]),
             ("Physical edge gap", ["short_edge_to_cue_gap"]),
             ("Gap / incumbent width", ["gap_over_width"]),
             ("Physical separation + width", ["short_route", "basal_width"]),
             ("Physical separation + C", ["short_route", "circumference"]),
             ("Edge gap + width", ["short_edge_to_cue_gap", "basal_width"])]
    predictions, folds, summaries = [], [], []
    for model in MODELS:
        g = frame[(frame.model == model) & (frame.status == "complete")].copy()
        y = g.outcome_code.to_numpy(int)
        C = g.circumference.to_numpy()
        for k in (3, 7, 15):
            for name, keys in specs:
                pred = np.zeros(len(g), dtype=int)
                for c in np.unique(C):
                    training, test = C != c, C == c
                    prior = np.bincount(y[training], minlength=5)
                    if not keys:
                        pred[test] = prior.argmax()
                    else:
                        X = g[keys].to_numpy(float)
                        sd = X[training].std(axis=0)
                        sd[sd < 1e-12] = 1
                        distances = (((X[test, None, :] - X[None, training, :]) / sd)**2).sum(axis=2)
                        for idx, distances_i in zip(np.flatnonzero(test), distances):
                            cut = np.partition(distances_i, min(k-1, len(distances_i)-1))[min(k-1, len(distances_i)-1)]
                            votes = np.bincount(y[training][distances_i <= cut + 1e-12], minlength=5)
                            tied = np.flatnonzero(votes == votes.max())
                            pred[idx] = tied[np.argmax(prior[tied])]
                    folds.append(dict(model=model, descriptor=name, k=k, held_out_C=int(c),
                                      test_conditions=int(test.sum()), training_conditions=int(training.sum()),
                                      accuracy=float(np.mean(pred[test] == y[test]))))
                for (_, row), actual, predicted in zip(g.iterrows(), y, pred):
                    predictions.append(dict(model=model, descriptor=name, k=k,
                        circumference=row.circumference, fraction=row.fraction,
                        actual=int(actual), predicted=int(predicted), correct=bool(actual == predicted)))
                summaries.append(dict(model=model, descriptor=name, k=k, conditions=len(g),
                    held_out_sizes=len(np.unique(C)),
                    equal_size_mean_accuracy=float(np.mean([np.mean(pred[C == c] == y[C == c]) for c in np.unique(C)])),
                    condition_weighted_accuracy=float(np.mean(pred == y)),
                    balanced_class_accuracy=float(np.mean([np.mean(pred[y == j] == j) for j in np.unique(y)]))))
    return pd.DataFrame(predictions), pd.DataFrame(folds), pd.DataFrame(summaries)


def setup():
    plt.rcParams.update({"font.size":9, "axes.titlesize":10, "axes.labelsize":9,
                         "xtick.labelsize":8, "ytick.labelsize":8, "svg.fonttype":"none",
                         "pdf.fonttype":42, "axes.spines.top":False, "axes.spines.right":False})


def ring_sketch(ax):
    for cx, C in ((-4.5, 30), (8, 40)):
        radius = C / (2*np.pi)
        ax.add_patch(Circle((cx, 0), radius, fill=False, color="#a5a5a5", lw=1.1))
        for d, color, ls in ((0, "#333333", "-"), (12, COLORS[1], "-")):
            angle = 2*np.pi*d/C
            point = np.array([cx, 0]) + radius * np.array([np.cos(angle), np.sin(angle)])
            ax.plot(*point, "o", color=color, ms=4)
            ax.annotate("", xy=point+1.6*np.array([np.cos(angle), np.sin(angle)]), xytext=point,
                        arrowprops=dict(arrowstyle="->", color=color, lw=1.1))
        ax.add_patch(Arc((cx, 0), 2*radius, 2*radius, theta1=0, theta2=360*12/C,
                         color=COLORS[1], lw=2))
        ax.text(cx, -.2, f"C = {C}\n{360*12/C:.0f}°", ha="center", va="center", fontsize=9)
        ax.text(cx, -radius-1.3, "d = 12", ha="center", color=COLORS[1], fontsize=9)
    ax.set_aspect("equal")
    ax.set(xlim=(-11.5, 17.2), ylim=(-8.5, 8.6))
    ax.axis("off")
    ax.set_title("A  Equal prescribed distance, different normal angles", loc="left", pad=12)
    ax.text(.5, -.04, r"$\Delta\theta=2\pi d/C$", transform=ax.transAxes, ha="center", fontsize=12)


def make_main(frame, histories):
    fig = plt.figure(figsize=(11.8, 9.3))
    fig.suptitle("Figure 3 · Geometry connects local polarity organization to an activity vector",
                 fontsize=15, y=.985)
    fig.text(.05, .946, "A larger perimeter converts the same physical separation into a smaller angular separation.", fontsize=10.8)
    a = fig.add_axes([.04, .59, .32, .31]); ring_sketch(a)
    b = fig.add_axes([.405, .615, .245, .265])
    for d, color, name in ((12, COLORS[3], "Spanning region"), (16, COLORS[2], "Distinct sites"),
                           (20, "#424242", "Opposing sites")):
        h = histories[(40, d)]
        physical_x = (h["x"]+10) % 40 - 10
        order = np.argsort(physical_x)
        b.plot(physical_x[order], h["p"][-1,order], color=color, lw=1.5, label=f"d = {d}: {name}")
        b.axvline(d, color=color, ls=":", lw=.7, alpha=.8)
    b.set(xlim=(-10, 30), ylim=(0, None), xlabel="Physical position from original peak", ylabel="Active Rac concentration")
    b.set_title("B  Placement changes local organization", loc="left", pad=13)
    b.legend(frameon=False, fontsize=7.5, loc="lower left", bbox_to_anchor=(-.03, .01))
    b.text(.02, -.30, "Holmes C40 · d: initial peak to cue centre\nReadout: 13,000 cue-free units", transform=b.transAxes, fontsize=8.5)
    c = fig.add_axes([.725, .615, .235, .265])
    g = frame[(frame.model == "holmes_model3") & (frame.circumference == 40) & (frame.status == "complete")]
    c.plot(g.angular_separation_deg, g.R_t14000, color="#444444", marker="o", ms=3.3, label=r"Relative vector magnitude $R$")
    c.plot(g.angular_separation_deg, g.q2_t14000, color="#A25F35", marker="s", ms=3, label=r"Second moment $q_2=|M_2|/M_0$")
    c.axhline(.1, color="#a5a5a5", ls=":", lw=.8)
    for d in (12,16,20):
        row = g[np.isclose(g.short_route, d)].iloc[0]
        c.annotate(f"d{d}", (row.angular_separation_deg,row.R_t14000), xytext=(4,5), textcoords="offset points", fontsize=8)
    c.set(xlim=(85, 185), ylim=(0, 1.16), xticks=[90, 120, 150, 180], xlabel="Initial-peak–cue normal angle (°)", ylabel="Relative moment magnitude")
    c.set_title("C  Opposing activity reduces the vector", loc="left", pad=13)
    c.legend(frameon=False, fontsize=7.5, loc="upper right")
    c.text(.02, -.30, "At 180°: two local peaks, R = 0.040\nDirection masked below R = 0.1", transform=c.transAxes, fontsize=8.5)
    fig.text(.05, .514, "The same physical distance or the same angle can accompany different native pattern fates.", fontsize=10.8)
    specs = [(30,12,COLORS[1],"Cue-associated single peak"),
             (40,12,COLORS[3],"One spanning region"),
             (40,16,COLORS[2],"Distinct old + cue peaks")]
    vmax = max(histories[(C,d)]["p"].max() for C,d,*_ in specs)
    for i, (C, d, color, label) in enumerate(specs):
        left = .06 + .31*i
        ax = fig.add_axes([left, .165, .25, .278])
        h = histories[(C,d)]
        keep = h["t"] <= 4000
        physical_x = (h["x"]+10) % C - 10
        order = np.argsort(physical_x)
        mesh = ax.pcolormesh(physical_x[order], h["t"][keep]/1000, h["p"][keep][:,order],
                             shading="nearest", cmap="magma", vmin=0, vmax=vmax, rasterized=True)
        ax.axvline(0, color="white", ls=":", lw=.8)
        ax.axvline(d, color="#75deeb", ls=":", lw=.8)
        ax.axhline(1, color="#75deeb", ls="--", lw=.9)
        ax.set(xlim=(-10,30), ylim=(4,0), yticks=[0,1,2,4], xlabel="Position from original peak", ylabel="Time / cue duration" if i==0 else "")
        if C < 40:
            ax.axvspan(C-10,30,color="#e1e1e1")
            ax.text(25,2.3,"Outside C30",ha="center",va="center",rotation=90,color="#666666",fontsize=8)
        ax.set_title(f"D{i+1}  C{C}, d{d}, {360*d/C:.0f}°\n{label}", color=color, loc="left", pad=12)
        row = frame[(frame.model=="holmes_model3") & (frame.circumference==C) & np.isclose(frame.short_route,d)].iloc[0]
        fig.text(left,.096, f"Final vector angle {row['angle_deg_t14000_gate_0p1']:.1f}°; R = {row.R_t14000:.3f}\nIncumbent width {row.basal_width:.2f}; edge gap {row.short_edge_to_cue_gap:.2f}", fontsize=8.5)
    cb = fig.colorbar(mesh, cax=fig.add_axes([.962,.165,.012,.278]))
    cb.set_label("Active Rac concentration", fontsize=8)
    fig.text(.05, .033, "Native Holmes histories · cue 0–1000 (cyan withdrawal), fixed full support 3 and peak 1. Field display ends at 4000; final readout is 14000.\nActivity vectors describe spatial activity, without assuming a cell-movement law. Complete maps and descriptor tests are linked in the evidence.", fontsize=8.5)
    export("revision_figure3_geometry", fig)


def make_support(frame, transfer):
    fig, axes = plt.subplots(3, 2, figsize=(11.5, 11), constrained_layout=True)
    for i, model in enumerate(MODELS):
        g = frame[frame.model == model]
        for j, descriptor in enumerate(("short_route", "angular_separation_deg")):
            ax = axes[i,j]
            for code, color in enumerate(COLORS):
                z = g[(g.status == "complete") & (g.outcome_code == code)]
                ax.scatter(z[descriptor],z.circumference,c=color,s=26,marker="s",edgecolor="none")
            excluded = g[g.status != "complete"]
            ax.scatter(excluded[descriptor], excluded.circumference, c="#cccccc", marker="x", s=15, linewidths=.7)
            ax.set(title=NAMES[model], xlabel="Prescribed initial-peak–cue distance d" if j==0 else "Initial-peak–cue normal angle (°)", ylabel="Circumference C")
    handles = [Patch(facecolor=color,label=label) for color,label in zip(COLORS,LABELS)]
    handles += [Line2D([],[],marker="x",ls="",color="#bbbbbb",label="Excluded / unqualified")]
    fig.legend(handles=handles,loc="outside lower center",ncol=3,frameon=False,fontsize=9)
    fig.suptitle("Chemical-only retained-ring support · placement outcomes\n652 complete / 867 requested; Spring counterpart is linked separately",fontsize=13)
    export("revision_figure3_all_outcomes",fig)
    fig, axes = plt.subplots(1,3,figsize=(11.5,4.4),sharey=True,constrained_layout=True)
    selected = ["Model majority","Angular separation","Physical separation","Physical separation + width"]
    for ax,model in zip(axes,MODELS):
        z=transfer[(transfer.model==model)&(transfer.k==7)].set_index("descriptor").loc[selected]
        values=z.equal_size_mean_accuracy.to_numpy()
        ax.barh(range(4),values,color=["#bdbdbd","#9F7C62","#568AAD","#5A9D88"])
        for i,v in enumerate(values):ax.text(v+.015,i,f"{v:.2f}",va="center",fontsize=9)
        ax.set(xlim=(0,1.08),yticks=range(4),yticklabels=selected,xlabel="Mean accuracy across held-out sizes",title=NAMES[model])
        ax.invert_yaxis()
    fig.suptitle("Chemical-only descriptor support · adding basal width improves held-out-size classification\nFixed 7-neighbour rule; Spring long-cue class is invariant and not ranked",fontsize=12)
    export("revision_figure3_descriptor_transfer",fig)


def main():
    setup()
    matrix = pd.read_csv(source(RUN / "results/matrix.csv"))
    moment = pd.read_csv(source(MOM / "results/matrix_summary.csv"))
    native_audit = json.loads(source(RUN / "results/audit.json").read_text())
    moment_audit = json.loads(source(MOM / "audit.json").read_text())
    assert native_audit["passed"] and moment_audit["passed"]
    source(RUN / "manifest.json")
    source(ROOT / "code/mechanisms/global_polarity_vector_20261017.py")
    source(__file__)
    assert len(matrix)==867 and (matrix.status=="complete").sum()==652
    assert not matrix.duplicated(["model","circumference","fraction"]).any()
    extra = [c for c in moment.columns if c not in matrix.columns or c in ("model","circumference","fraction")]
    frame = matrix.merge(moment[extra],on=["model","circumference","fraction"],validate="one_to_one").copy()
    expected = json.loads(source(RUN / "results/source_sha256.json").read_text())
    assert sha(MOM / "results/matrix_summary.csv") == moment_audit["outputs_sha256"]["results/matrix_summary.csv"]
    # Planned overlap cases still have a realized distance. Unqualified starts
    # have only a requested fraction; preserve that distinction explicitly.
    missing = frame.short_route.isna()
    frame["distance_for_geometry"] = frame.short_route.fillna(frame.circumference * frame.fraction)
    frame["geometry_coordinate_status"] = np.where(missing,"requested, no qualified native placement","realized native placement")
    frame["angular_separation_deg"] = 360*frame.distance_for_geometry/frame.circumference
    frame["outcome_code"] = frame.final_location_association.map(CLASSES)
    frame["gap_over_width"] = frame.short_edge_to_cue_gap/frame.basal_width
    complete = frame[frame.status=="complete"]
    assert complete.outcome_code.notna().all()
    for _, group in complete.groupby(["model", "circumference"]):
        assert group.checkpoint_sha256.nunique() == 1
        assert group.basal_width.nunique() == 1
    np.testing.assert_allclose(complete.angular_separation_deg,360*complete.realized_fraction,atol=1e-12)
    table("geometry_placement_outcomes.csv",frame)
    within = complete.groupby(["model","circumference"]).agg(conditions=("fraction","size"),outcome_classes=("outcome_code","nunique"),basal_width=("basal_width","first")).reset_index()
    within["physical_angular_rank_identical"] = True
    within["identical_basal_checkpoint_across_placements"] = True
    table("geometry_within_size_controls.csv",within)
    predictions,folds,transfer=descriptor_transfer(frame)
    table("geometry_descriptor_predictions.csv",predictions)
    table("geometry_descriptor_held_out_sizes.csv",folds)
    table("geometry_descriptor_transfer_summary.csv",transfer)
    histories={}
    selected=[]
    checks=[]
    for C,d in ((30,12),(40,12),(40,16),(40,20)):
        row=complete[(complete.model=="holmes_model3")&(complete.circumference==C)&np.isclose(complete.short_route,d)].iloc[0]
        h=history(row);histories[(C,d)]=h
        assert SOURCES[str(Path(row.cue_path).resolve())] == expected[str(Path(row.cue_path))]
        assert SOURCES[str(Path(row.recovery_path).resolve())] == expected[str(Path(row.recovery_path))]
        check={"C":C,"d":d}
        for t in (1000,4000,9000,14000):
            idx,=np.flatnonzero(h["t"]==t)
            for key in ("R","q1","q2","M0"):
                error=float(abs(h[key][idx]-row[f"{key}_t{t}"]))
                assert error<1e-10
                check[f"{key}_t{t}_absolute_error"]=error
        selected.append(row.to_dict());checks.append(check)
        np.savez_compressed(dest(f"Holmes_C{C}_d{d}_display.npz"),**h)
    for d in (16,20):
        np.testing.assert_array_equal(histories[(40,12)]["p"][0],histories[(40,d)]["p"][0])
    table("geometry_selected_cases.csv",pd.DataFrame(selected))
    make_main(frame,histories)
    make_support(frame,transfer)
    write_docs(frame,within,transfer,checks)


def write_docs(frame,within,transfer,checks):
    generated=datetime.now(timezone.utc).isoformat()
    result_lines=[]
    for model in MODELS:
        w=within[within.model==model]
        z=transfer[(transfer.model==model)&(transfer.k==7)].set_index("descriptor")
        result_lines.append(f"- {NAMES[model]}: placement changes the finite class in {(w.outcome_classes>1).sum()}/{len(w)} qualified sizes; held-out-size accuracy is {z.loc['Angular separation','equal_size_mean_accuracy']:.3f} for angle, {z.loc['Physical separation','equal_size_mean_accuracy']:.3f} for physical separation and {z.loc['Physical separation + width','equal_size_mean_accuracy']:.3f} for physical separation plus basal width.")
    caption="""Figure 3. Geometry connects local polarity organization to an activity vector.

A, Circular boundaries drawn at a common physical scale. A prescribed shorter-arc separation d of 12 native length units from the initial incumbent peak to the challenger cue centre gives an outward-normal difference of 144° at circumference 30 and 108° at circumference 40. B, Native Holmes Rac profiles 13,000 units after withdrawal at C40, from the identical qualified basal state, with a raised-cosine challenger of full support 3, native peak 1 and duration 1000. Increasing prescribed cue-centre separation changes a connected spanning region (d12) to distinct old- and cue-associated peaks (d16 and d20). Dotted vertical lines mark challenged positions. C, Normalized first-angular-moment magnitude R and second-moment ratio q₂ from every eligible C40 placement at the same readout, plotted against the normal angle between prescribed initial-peak and cue-centre positions. Opposing local peaks can nearly cancel the first moment: at d20 the two resolved peaks coexist with R=0.0399 and q₂=0.4006. The direction of the activity vector is masked below R=0.1; a small vector is not evidence of absent local activity. D1–D3, Exact native histories selected in the preceding checkpoint: C30,d12 gives a cue-associated single peak; C40,d12 gives one connected region spanning old and cue; C40,d16 retains two distinct associated peaks. Thus matching prescribed physical separation (D1/D2) or prescribed angular separation (D1/D3) does not match finite outcomes when circumference and incumbent architecture also change. D2/D3 provide the controlled within-size placement contrast. Cyan marks cue withdrawal; histories end at 4000 for display while printed readouts use time 14000.

The activity first moment integrates local activity weighted by the outward normal. P is the complex first moment divided by its initial complex value, R=|P|, and q₂ is the absolute second angular moment divided by total active amount. The first moment is an activity summary, without an assumed cell-velocity law or a unique active-site identity. Region association is a spatial endpoint descriptor, not a lineage assignment or proof of asymptotic coexistence. Geometry coordinates and chemical competition are related measurements rather than a causal partition.
"""
    dest("CAPTION.md","metadata").write_text(caption)
    evidence=r"""# Figure 3: claims and retained evidence

Generated: {generated}. Reanalysis finished; no simulations were commissioned.

The geometry claim is exact: circular shorter-arc separation d maps to outward-normal angle Δθ=2πd/C. In the competition grid d is the prescribed distance from the basal incumbent peak to the challenger cue centre, and the angle is between those prescribed positions. Neither is an observed separation between evolved sites. This adds a continuously valued two-dimensional directional consequence to a spatial polarity pattern; the established interval experiments do not supply these circular normal directions. This is a coordinate consequence, not an extra biochemical mechanism or a cell-movement prediction.

The controlled competition claim is narrower: in Holmes C40, changing only challenger cue-centre position from d12 to d16 changes the finite endpoint from one spanning region to two associated peaks. Original basal, local parameters, input peak/support, input duration and readout time agree. The C30,d12 comparison additionally changes basal organization and cannot isolate geometry from chemistry. Same-distance and same-angle contrasts show that neither one scalar coordinate completely describes the grid.

The circular-vector consequence is supported by the full native C40,d20 profile and independently replayed moments: two antipodal activity peaks have normalized first-moment magnitude 0.03993 while q₂ is 0.40065. This is cancellation of a spatial moment. The direction is unavailable by the prespecified R=0.1 gate, and the local field remains active.

## Descriptor analysis, including useful negative results

All 652 complete conditions are retained from 867 requested placements: 198 overlap exclusions and 17 unqualified preparations. Exclusions are not missing outcomes or zero activity. Within any fixed circumference d and Δθ are exact positive multiples, hence their ranks and one-dimensional ordering are identical. No within-size regression can distinguish their independent causal influence. The grid contains 17 circumferences in Wave-pinning/Goryachev and 16 with qualified, nonoverlapping Holmes outcomes.

For a descriptive transfer check, fit the same fixed 7-neighbour classifier separately within each model while withholding every placement at one entire circumference. Distances use only the training set's coordinate standard deviations. All tied neighbours at the kth distance enter; vote ties use the training class frequency and then stable class order. No tuning or random split is used. Accuracy is averaged equally over held-out sizes; condition-weighted and balanced-class accuracies, fold predictions, and k3/k15 sensitivity are exported. Deterministic placements are not independent biological samples and these accuracies carry no sampling confidence interval or significance test.

{result_lines}

Physical separation transfers better than angle in this fixed analysis, but adding incumbent width improves transfer further. This is not causal superiority of physical distance: circumference covaries with available pool, incumbent organization and the unmeasured interactions of route lengths. The equally accessible physical-separation-plus-C control performs similarly; width is not isolated by prediction. Finite class rarity and excluded placements affect transfer, and the classifier does not establish asymptotic stability. No ring-vs-interval matched intervention identifies which competition mechanisms are uniquely circular. Therefore the classifier remains supplemental and the main figure focuses on exact geometry, controlled position contrasts and the activity-moment consequence.

## Definitions

For a uniform native ring with activity a(s), define its angular moments by

$$ M_n(t)=\int_0^C a(s,t)\exp\!\left(\mathrm{{i}}n\frac{{2\pi s}}{{C}}\right)\,\mathrm{{d}}s. \tag{{1}} $$

Use P(t)=M₁(t)/M₁(0), R(t)=|P(t)| and qₙ(t)=|Mₙ(t)|/M₀(t). The angle is reported only where R≥0.1, following the reviewed existing analysis. The native finite-volume quadrature is retained, with the original peak cell as angular origin. Uniform background contributes zero first moment. Neither a global angle nor a local maximum identifies cell movement.

## Reproduction and retained sources

Run `OPENBLAS_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/revision_figure3_20261019.py` from repository root. This reads audited existing simulation and moment tables, replays four selected full native histories, validates exact phase continuity and four time-point moment summaries, exports complete placement and held-out-size tables and renders main/support sheets. Work is serial, below 1 GB resident memory and requires no new integrations.

The prior gain/recovery and standardized finite-cue ring sheets remain applicable supporting evidence, not replaced measurements. Sources include `dev/vision_aligned_checkpoint_20261017/vision_aligned_figure2_response.png`, `vision_aligned_figure2_competition.png`, and the common recovery tables under `data/runs/vision_ring_recovery_20261017/results`. These support local response measurements while the revised main figure concentrates on geometry.

See `validation.json` for numerical checks, measured allocation and artifact hashes, and `source_sha256.json` for immutable input hashes.
""".format(generated=generated,result_lines="\n".join(result_lines))
    duration_evidence=OUT/"distance_response/supporting/metadata/DURATION_EVIDENCE.md"
    if duration_evidence.exists():
        evidence += "\n## Subsequent four-model distance-duration revision\n\nThis geometry sheet is the first of two main sheets. The [four-model continuation](../../distance_response/revision_figure3_distance_duration.png) adds controlled finite-duration responses and prominent native Spring profile evidence. Its [evidence and limits](../../distance_response/supporting/metadata/DURATION_EVIDENCE.md), [150 requested conditions with actual archive provenance](../../distance_response/supporting/data/distance_duration_all_models.csv) and [endpoint sensitivity table](../../distance_response/supporting/data/distance_duration_shortest_tested_sensitivity.csv) distinguish retained angle, vector magnitude and fixed-region activity loss. The [complete Spring size/placement sheet](../../spring_distance/spring_distance_duration.png) and [Spring field/geometry audit](../../spring_distance/supporting/metadata/validation.json) retain the full three-size extension. The first sheet's native competition illustration is explicitly Holmes-specific; no all-model comparison is inferred from those selected histories.\n"
        evidence += "\n## Scope of historical chemical support and the current Spring counterpart\n\nThe retained 652-condition placement grid, held-out-size descriptor fits and exact interval/ring reflection maps cover Wave-pinning, Goryachev and Holmes only. They are explicitly chemical-only historical support. The [Spring standard-cue coverage record](../../spring_distance/supporting/metadata/descriptor_coverage.json) and [14-row table](../../spring_distance/supporting/data/spring_standard_long_cue_descriptor_coverage.csv) evaluate the same finite-duration question at all three qualified Spring sizes. All 14 eligible T=1000 geometries yield one cue-associated peak; continuous angle and old-core amount variation remain exported, but categorical descriptor ranking is uninformative for a constant class. No generic classifier accuracy is reported.\n\nExact nonuniform interval reflection is not applicable to the declared Spring mean-force closure: the native interval samples endpoint Rac forces, whereas the ring averages the nonlinear force around the perimeter. The validated uniform-state derivative and native stretching/dilution/diffusion carryover do not equate these nonuniform mechanical rules. The current Spring counterparts are its three-size basal organization/function measurements and its actual distance-duration trajectories, in explicitly native units.\n"
    dest("EVIDENCE.md","metadata").write_text(evidence)
    repo_path=lambda p:str(p.relative_to(REPO))
    content=dict(title="Figure 3. Geometry connects local polarity organization to an activity vector",
        role="Exact physical-to-angular geometry, controlled site placement and cancellation of the activity first moment.",
        assets=[dict(path=repo_path(OUT/"revision_figure3_geometry.png"),caption="Revised compact geometry sheet; all plotted fields are retained native simulations.",page_size=[11.8,10.4],page_title="Figure 3. Geometry and the activity vector",minimum_printed_width_inches=10.0)],
        legend=[p for p in caption.strip().split("\n\n")][1:],
        narrative=["Equal physical separation gives a smaller outward-normal angular difference on a larger ring.","At fixed Holmes circumference, changing challenger position alters the finite spatial endpoint from one spanning region to distinct old- and cue-associated peaks.","Strong, opposing local activity can yield a weak global activity vector. The activity angle requires local-pattern and magnitude context."],
        measurements=["Measure physical perimeter and the arc separation of local activity/protrusion regions, and compare their outward-normal angle separately from their local activity strength.","Track local activity profiles alongside first- and second-angular-moment summaries; a weak first moment can accompany multiple strong regions. A cell-velocity law is not assumed."],
        limitations=["Physical and angular separation are collinear at fixed circumference. Cross-size comparisons also change incumbent architecture and pool organization, and cannot isolate causal geometry.","Region association is a finite spatial descriptor; no unique site lineage or asymptotic coexistence is established."],
        sources=[dict(label="Geometry evidence and negative results",path=repo_path(OUT/"supporting/metadata/EVIDENCE.md")),dict(label="Complete placement outcomes and moments",path=repo_path(OUT/"supporting/data/geometry_placement_outcomes.csv")),dict(label="Held-out-size descriptor transfer",path=repo_path(OUT/"supporting/data/geometry_descriptor_transfer_summary.csv")),dict(label="Retained gain/recovery ring sheet",path=repo_path(ROOT/"figures/length_mechanisms_20261008/dev/vision_aligned_checkpoint_20261017/vision_aligned_figure2_response.png")),dict(label="Numerical and provenance validation",path=repo_path(OUT/"supporting/metadata/validation.json"))],
        supplementary=[dict(title="Historical chemical-only placement grid; Spring counterpart linked separately",path=repo_path(OUT/"revision_figure3_all_outcomes.png")),dict(title="Historical chemical-only descriptor transfer; no invariant-class Spring score",path=repo_path(OUT/"revision_figure3_descriptor_transfer.png"))])
    if duration_evidence.exists():
        content["sources"].extend([
            dict(label="Four-model finite-duration continuation and endpoint limits",path=repo_path(duration_evidence)),
            dict(label="Complete duration grid with actual native archive paths",path=repo_path(OUT/"distance_response/supporting/data/distance_duration_all_models.csv")),
            dict(label="Independent chemical/Spring retained-profile checks",path=repo_path(OUT.parent/"supporting/metadata/REQUIREMENT_EVIDENCE_AUDIT.json")),
            dict(label="Spring descriptor coverage and nonuniform interval-correspondence scope",path=repo_path(OUT/"spring_distance/supporting/metadata/descriptor_coverage.json"))])
        content["supplementary"].append(dict(title="Current Spring long-cue invariant class and continuous descriptor coverage",path=repo_path(OUT/"spring_distance/supporting/data/spring_standard_long_cue_descriptor_coverage.csv")))
        for item in content["sources"]:
            if item["label"]=="Retained gain/recovery ring sheet":
                item["label"]="Historical chemical-only gain/recovery and interval/ring sheet"
    jsonwrite("content.json",content)
    jsonwrite("source_sha256.json",SOURCES)
    allocated=sum(p.stat().st_blocks*512 for p in OUT.rglob("*") if p.is_file())
    validation=dict(generated_at=generated,passed=True,new_integrations=0,planned_conditions=867,complete_conditions=652,
        excluded_conditions=frame[frame.status!="complete"].groupby("status").size().to_dict(),
        exact_phase_continuity=True,selected_source_hashes_match_archived=True,moment_summary_hash_matches_archived=True,selected_native_history_checks=checks,
        angular_coordinate_identity_max_error=float(abs(frame[frame.status=="complete"].angular_separation_deg-360*frame[frame.status=="complete"].realized_fraction).max()),
        descriptor_configurations=int(len(transfer)),held_out_size_folds=1200,classifier_primary_k=7,classifier_sensitivity_k=[3,15],
        scientific_limits="Deterministic descriptive transfer, not independent biological replication, causal partition or movement prediction.",
        output_allocated_bytes_before_validation=int(allocated),large_output_manifest_required=bool(allocated>=2*1024**3),
        artifacts_sha256={repo_path(p):sha(p) for p in ARTIFACTS if p.exists()})
    jsonwrite("validation.json",validation)
    print(json.dumps(dict(passed=True,output=str(OUT),allocated_bytes=allocated,artifacts=len(ARTIFACTS)),indent=2))


if __name__=="__main__":
    main()
