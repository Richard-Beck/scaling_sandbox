"""Two-state teaching sketch from exact retained Holmes ring profiles; no solver."""
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
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

REPO = Path(__file__).resolve().parents[4]
BASE = REPO / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/standard_angles"
CONCEPT = BASE / "vector_concept"
OUT = CONCEPT / "teaching_example"
N_BINS = 12
BLUE = "#327aa5"
RED = "#c64250"
GREEN = "#338369"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def csvwrite(path, values):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


def main():
    for folder in (OUT, OUT / "supporting/data", OUT / "supporting/metadata", OUT / "supporting/vectors"):
        folder.mkdir(parents=True, exist_ok=True)
    old_files = {str(p): sha(p) for p in CONCEPT.rglob("*") if p.is_file() and OUT not in p.parents}
    summary = BASE / "supporting/data/standard_angles_summary.csv"
    selected_path = CONCEPT / "supporting/data/selected_case_measurements.csv"
    profiles_path = CONCEPT / "supporting/data/selected_native_profiles.csv"
    row, = [r for r in rows(summary) if r["model"] == "holmes_model3" and r["angle_label"] == "pi" and float(r["cue_duration"]) == 1000.]
    selected, = [r for r in rows(selected_path) if r["model"] == "holmes_model3"]
    assert row["condition_id"] == selected["condition_id"] == "44"
    C = float(row["circumference"])
    assert C == 40.
    cue = Path(row["actual_cue_path"])
    recovery = Path(row["actual_recovery_path"])
    assert sha(cue) == row["actual_cue_sha256"]
    assert sha(recovery) == row["actual_recovery_sha256"]
    with np.load(cue) as z:
        x = z["x"].copy()
        initial = z["states"][0, 0].copy()
        withdrawal = z["states"][-1, 0].copy()
        assert z["time"][0] == 0. and z["time"][-1] == 1000.
    with np.load(recovery) as z:
        np.testing.assert_array_equal(z["x"], x)
        np.testing.assert_array_equal(z["states"][0, 0], withdrawal)
        final = z["states"][-1, 0].copy()
        assert z["time"][0] == 0. and z["time"][-1] == 13000.
    n = len(x)
    ds = C / n
    theta = 2 * np.pi * x / C
    normal = np.exp(1j * theta)
    fields = np.stack((initial, final))
    assert np.isrealobj(fields) and np.all(fields >= 0)
    cell_contribution = fields * ds * normal
    M1 = cell_contribution.sum(axis=1)
    M0 = fields.sum(axis=1) * ds
    P = M1 / M1[0]
    R = abs(P)
    normalized_cell = cell_contribution / M1[0]
    angle0 = np.angle(M1[0])
    alpha = np.angle(np.exp(1j * (theta-angle0)))
    bin_index = (np.floor((alpha + np.pi/N_BINS)/(2*np.pi/N_BINS)).astype(int)) % N_BINS
    binned = np.array([[normalized_cell[k, bin_index == j].sum() for j in range(N_BINS)] for k in range(2)])
    bin_error = float(abs(binned.sum(axis=1)-P).max())
    assert bin_error < 1e-14
    assert abs(M0[0]-float(selected["initial_M0"])) < 1e-12
    assert abs(M1[0]-(float(selected["initial_M1_real"])+1j*float(selected["initial_M1_imag"]))) < 1e-12
    assert abs(R[1]-float(selected["final_R"])) < 1e-12
    assert abs(M0[1]/M0[0]-float(selected["final_M0_ratio"])) < 1e-12
    profiles = rows(profiles_path)
    profile_error = 0.
    for k, snapshot in enumerate(("initial", "final")):
        rr = [r for r in profiles if r["model"] == "holmes_model3" and r["snapshot"] == snapshot]
        assert len(rr) == n
        np.testing.assert_array_equal(np.array([float(r["native_active"]) for r in rr]), fields[k])
        profile_error = max(profile_error, float(abs(np.array([float(r["physical_ds"]) for r in rr])-ds).max()))
    assert profile_error < 1e-14

    cells, bins, stats = [], [], []
    for k, state in enumerate(("initial", "final")):
        for j in range(n):
            c, p = cell_contribution[k, j], normalized_cell[k, j]
            cells.append(dict(state=state, cell=j, physical_position=float(x[j]), physical_segment_length=ds,
                              outward_normal_angle_raw_rad=float(theta[j]), outward_normal_angle_initial_aligned_rad=float(alpha[j]),
                              native_active_real=float(fields[k, j]), contribution_real=float(c.real), contribution_imag=float(c.imag),
                              contribution_initial_normalized_real=float(p.real), contribution_initial_normalized_imag=float(p.imag), bin=int(bin_index[j])))
        for j in range(N_BINS):
            c = binned[k, j]
            bins.append(dict(state=state, bin=j, centre_initial_aligned_degrees=360*j/N_BINS,
                             native_cells=int(np.sum(bin_index == j)), physical_bin_length=float(np.sum(bin_index == j)*ds),
                             activity_amount=float((fields[k, bin_index == j]*ds).sum()),
                             contribution_initial_normalized_real=float(c.real), contribution_initial_normalized_imag=float(c.imag)))
        stats.append(dict(state=state, time_after_cue_on_native=0. if k == 0 else 14000., time_after_withdrawal_native=None if k == 0 else 13000.,
                          M0_physical=float(M0[k]), M0_over_initial=float(M0[k]/M0[0]), M1_physical_real=float(M1[k].real), M1_physical_imag=float(M1[k].imag),
                          P_real=float(P[k].real), P_imag=float(P[k].imag), R=float(R[k]), heading_interpreted=bool(R[k] >= .1),
                          native_activity_min=float(fields[k].min()), native_activity_max=float(fields[k].max())))
    csvwrite(OUT / "supporting/data/native_cell_contributions.csv", cells)
    csvwrite(OUT / "supporting/data/exact_bin_sums.csv", bins)
    csvwrite(OUT / "supporting/data/two_state_measurements.csv", stats)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(9.5, 6.4), facecolor="white")
    fig.text(.5, .966, "Two activity peaks can leave a weak resultant", ha="center", va="top", fontsize=17, weight="bold")
    fig.text(.5, .913, "Activity u is real; x + iy stores its weighted normal’s two coordinates.", ha="center", fontsize=10.4)
    fig.text(.5, .864,
             r"$M_1=\sum_j u_j\,\Delta s_j\,(\cos\theta_j+i\sin\theta_j),\qquad M_0=\sum_j u_j\,\Delta s_j$",
             ha="center", fontsize=16)
    fig.text(.5, .792, r"$P=M_1/M_1(0)$ rotates and scales the initial resultant to 1.", ha="center", fontsize=11)
    norm = Normalize(0., float(fields.max()))
    cmap = plt.get_cmap("viridis")
    vector_scale = .83
    axes = [fig.add_axes((.025,.278,.455,.465)), fig.add_axes((.52,.278,.455,.465))]
    for k, ax in enumerate(axes):
        # Every coloured arc is the measured value of one physical native cell.
        half = np.pi/n
        points0 = np.column_stack((np.cos(alpha-half), np.sin(alpha-half)))
        points1 = np.column_stack((np.cos(alpha+half), np.sin(alpha+half)))
        collection = LineCollection(np.stack((points0, points1), axis=1), array=fields[k], cmap=cmap, norm=norm, linewidth=13, capstyle="butt", zorder=3)
        ax.add_collection(collection)
        ax.plot([-1.38,1.38], [0,0], color=".9", lw=.6, zorder=0)
        ax.plot([0,0], [-1.38,1.38], color=".9", lw=.6, zorder=0)
        ax.text(1.4,.02,"x",fontsize=9,color=".55")
        ax.text(.02,1.35,"y",fontsize=9,color=".55")
        for j in range(N_BINS):
            c = binned[k,j]
            tail = 1.105*np.exp(2j*np.pi*j/N_BINS)
            arrow = vector_scale*c
            length = abs(arrow)
            ax.arrow(tail.real,tail.imag,arrow.real,arrow.imag,
                     width=min(.009,.18*length),head_width=min(.055,.55*length),
                     head_length=min(.075,.45*length),length_includes_head=True,
                     color=BLUE,zorder=4)
        p = vector_scale*P[k]
        headscale = min(1., abs(P[k])*5)
        ax.arrow(0,0,p.real,p.imag,width=.014,head_width=.075*headscale,head_length=.11*headscale,
                 length_includes_head=True,color=RED,zorder=6)
        ax.scatter([0],[0],s=10,color=RED,zorder=7)
        ax.text(0,.28,"resultant P",ha="center",color=RED,fontsize=10.5)
        ax.set_xlim(-1.52,1.52);ax.set_ylim(-1.52,1.52);ax.set_aspect("equal");ax.axis("off")
        label = "A   Initial basal activity" if k == 0 else "B   Retained endpoint: opposing peaks"
        ax.set_title(label,fontsize=11,weight="bold",pad=7)
    fig.text(.25,.255,"Before the antipodal cue",ha="center",fontsize=9.2,color=".35")
    fig.text(.75,.255,"1000 cue units + 13,000 cue-free units",ha="center",fontsize=9.2,color=".35")
    for k, x0 in enumerate((.25,.75)):
        fig.text(x0,.207,rf"Total activity  $M_0/M_0(0)={M0[k]/M0[0]:.4f}$",ha="center",fontsize=11.5,color=GREEN)
        fig.text(x0,.163,rf"Resultant magnitude  $R=|P|={R[k]:.4f}$",ha="center",fontsize=11.5,color=RED)
    fig.text(.75,.124,"Heading is not interpreted when R < 0.1.",ha="center",fontsize=8.6,color=".4")
    fig.text(.25,.124,"Initial normalization sets P(0) = 1.",ha="center",fontsize=8.6,color=".4")
    cax = fig.add_axes((.375,.067,.25,.016))
    cb = fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation="horizontal")
    cb.set_label("Native active concentration u (shared scale)",fontsize=8.4,labelpad=2)
    cb.ax.tick_params(labelsize=8,pad=1)
    fig.text(.018,.066,"Blue arrows: exact sums of native\ncell contributions in 30° sectors.",fontsize=8.3,va="top",color=BLUE)
    fig.text(.69,.066,"All arrows use one vector scale.\nHolmes 3 · fixed circumference C = 40.",fontsize=8.3,va="top",color=".35")
    artifact_paths = []
    for ext in ("png","pdf","svg"):
        p = (OUT if ext == "png" else OUT / "supporting/vectors") / ("holmes_two_state_vector."+ext)
        fig.savefig(p,dpi=240,facecolor="white")
        artifact_paths.append(p)
    plt.close(fig)

    caption = """# Two-state activity-vector teaching example

Holmes 3, fixed circumference 40, antipodal stationary native cue of duration 1000, followed by 13000 cue-free native units. Both rings show exact saved active-concentration profiles on one common zero-based native intensity scale. The circle represents the model ring, not a measured cell outline. Each coloured arc is one native cell; its physical segment length is 40/800=0.05. Activity is a real scalar. Multiplying it by segment length and the outward unit normal gives real x and y contributions, written together as a complex number x+iy. Their sum is the first physical activity moment M₁; the scalar total activity is M₀. Complex division by the initial M₁ rotates and scales the initial resultant to P(0)=1. It is not division by current total activity.

Outer blue arrows sum the exact native complex contributions in twelve 30° sectors after that initial normalization. Bin populations need not be equal; native cell weights and directions are retained. Arrows are translated to their sector positions for illustration. The red arrow is their vector sum drawn at the origin. Blue and red arrows have one shared vector scale, identical in both states. The final short red arrow is the measured weak resultant; its heading is not interpreted because R<0.1. Initial and final physical total active amounts are 33.1293111 and 33.0450818; their ratio is 0.9974576, whereas final R is 0.0399348. Thus nearly retained total activity can accompany nearly cancelling activity directions. M₀ is total active amount, not the conserved active-plus-inactive pool. No velocity, cell displacement, site ancestry or asymptotic-stability claim is assigned.
"""
    # Substitute recomputed amounts rather than retaining rounded hand entries.
    caption = caption.replace("33.0450818",f"{M0[1]:.7f}")
    (OUT / "supporting/metadata/CAPTION.md").write_text(caption)
    readme = r"""# A two-state worked example of the activity vector

[PNG](holmes_two_state_vector.png) · [PDF](supporting/vectors/holmes_two_state_vector.pdf) · [SVG](supporting/vectors/holmes_two_state_vector.svg)

This small teaching sketch uses one already selected native Holmes history. It adds no integration and leaves the preceding concept sheets and Figure 3 artwork unchanged. The measured profile is shown as coloured activity on a schematic circle. Activity itself is real; complex numbers package its two Cartesian directional sums:

$$
M_1=\sum_j u_j\Delta s_j\cos\theta_j+i\sum_j u_j\Delta s_j\sin\theta_j,
\qquad M_0=\sum_j u_j\Delta s_j.
\tag{1}
$$

$$
P(t)=\frac{M_1(t)}{M_1(0)},\qquad R(t)=|P(t)|,\qquad\Delta s_j=\frac{40}{800}=0.05.
\tag{2}
$$

Initial complex normalization rotates and scales the vector so P(0)=1. Opposing activity can cancel in M₁ while adding positively in M₀. The final R≈0.040 and M₀/M₀(0)≈0.9975 demonstrate that distinction in the actual saved profile. Raw u is neither complex nor independently normalized between states. This is an activity descriptor, without a movement law.

[Caption and limits](supporting/metadata/CAPTION.md) · [Two-state values](supporting/data/two_state_measurements.csv) · [Every physical cell contribution](supporting/data/native_cell_contributions.csv) · [Exact sector sums](supporting/data/exact_bin_sums.csv) · [Sources and validation](supporting/metadata/validation.json)
"""
    (OUT / "README.md").write_text(readme)
    for path in (OUT / "README.md",OUT / "supporting/metadata/CAPTION.md",OUT / "supporting/data/native_cell_contributions.csv",OUT / "supporting/data/exact_bin_sums.csv",OUT / "supporting/data/two_state_measurements.csv"):
        artifact_paths.append(path)
    sources = [summary, selected_path, profiles_path, cue, recovery, Path(__file__)]
    assert all(sha(Path(p)) == h for p,h in old_files.items())
    metadata = dict(generated_at_utc=datetime.now(timezone.utc).isoformat(), passed=True,
                    scope="separate pedagogical sketch from two exact retained native states; no integrations or main-figure rewrite",
                    selection=dict(model="holmes_model3",condition_id="44",circumference=40.,angle_degrees=180.,cue_duration_native=1000.,withdrawal_native=13000.),
                    native_cells=n,physical_ds=ds,initial_rotation_radians=float(angle0),
                    initial_M1_physical_complex=[float(M1[0].real),float(M1[0].imag)],shared_native_activity_limits=[0.,float(fields.max())],
                    bin_count=N_BINS,bin_width_degrees=360/N_BINS,bin_sum_max_complex_error=bin_error,
                    arrows="exact sums of cell contribution/M1(initial) in initial-aligned 30deg sectors; common plotted vector scale 0.83, translated tails",
                    source_profiles_match_selected_csv=True,physical_ds_match_selected_csv=True,
                    measurements=stats,existing_concept_files_unchanged=len(old_files),
                    sources={str(p.relative_to(REPO)):sha(p) for p in sources},artifacts={str(p.relative_to(OUT)):sha(p) for p in artifact_paths})
    (OUT / "supporting/metadata/validation.json").write_text(json.dumps(metadata,indent=2)+"\n")
    print(json.dumps(dict(output=str(OUT),R_final=float(R[-1]),M0_final_over_initial=float(M0[-1]/M0[0]),
                          bin_sum_error=bin_error,old_files_unchanged=len(old_files)),indent=2))


if __name__ == "__main__":
    main()
