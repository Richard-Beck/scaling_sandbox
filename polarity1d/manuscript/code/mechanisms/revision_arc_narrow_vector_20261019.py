"""Measured activity-vector construction within the new stationary-cue experiment."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from revision_standard_angles_report_20261019 import periodic_counts

ROOT=Path(__file__).resolve().parents[2]
TOPIC=ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/narrow_arc_experiment"
OUT=TOPIC/"vector_construction"


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    for leaf in ("data","metadata","vectors"):(OUT/"supporting"/leaf).mkdir(parents=True,exist_ok=True)
    path=TOPIC/"supporting/data/narrow_arc_summary.csv"
    row=pd.read_csv(path).query("condition_id==383").iloc[0]
    assert row.model=="holmes_model3" and row.circumference==30 and row.cue_duration==300 and row.physical_distance_at_onset==14
    assert row.footprint_relative_to_core_10=="outside"
    cue=Path(row.actual_cue_path);assert sha(cue)==row.actual_cue_sha256
    with np.load(cue) as z:
        x=z["x"].copy();fields=z["states"][[0,-1],0].copy();time=z["time"][[0,-1]].copy()
    assert time.tolist()==[0.,300.]
    C=30.;n=len(x);ds=C/n;theta=2*np.pi*x/C
    physical=fields*ds*np.exp(1j*theta);M1=physical.sum(axis=1);M0=fields.sum(axis=1)*ds
    P=M1/M1[0];R=abs(P);alpha=np.angle(np.exp(1j*(theta-np.angle(M1[0]))));normalized=physical/M1[0]
    bins=(np.floor((alpha+np.pi/12)/(2*np.pi/12)).astype(int))%12
    sectors=np.array([[normalized[k,bins==j].sum() for j in range(12)] for k in range(2)])
    error=float(abs(sectors.sum(axis=1)-P).max());assert error<1e-14
    assert abs(R[1]-row.R_withdrawal)<1e-13 and abs(M0[1]/M0[0]-row.M0_ratio_withdrawal)<1e-13
    cells=[];stats=[];grouped=[]
    for k,label in enumerate(["initial","cue_withdrawal"]):
        count,regions=periodic_counts(fields[k])
        stats.append(dict(state=label,time=float(time[k]),M0=float(M0[k]),M0_ratio=float(M0[k]/M0[0]),M1_real=float(M1[k].real),M1_imag=float(M1[k].imag),
                          P_real=float(P[k].real),P_imag=float(P[k].imag),R=float(R[k]),angle_deg=float(np.degrees(np.angle(P[k]))),peaks=count,halfcontrast_regions=regions))
        cells.extend(dict(state=label,cell=j,physical_position=float(x[j]),active=float(fields[k,j]),physical_ds=ds,normal_angle_deg=float(np.degrees(alpha[j])),
                          M1_contribution_real=float(physical[k,j].real),M1_contribution_imag=float(physical[k,j].imag),P_contribution_real=float(normalized[k,j].real),P_contribution_imag=float(normalized[k,j].imag)) for j in range(n))
        grouped.extend(dict(state=label,sector=j,normalized_contribution_real=float(sectors[k,j].real),normalized_contribution_imag=float(sectors[k,j].imag)) for j in range(12))
    pd.DataFrame(cells).to_csv(OUT/"supporting/data/native_cell_contributions.csv",index=False)
    pd.DataFrame(stats).to_csv(OUT/"supporting/data/two_state_measurements.csv",index=False)
    pd.DataFrame(grouped).to_csv(OUT/"supporting/data/exact_sector_sums.csv",index=False)
    pd.DataFrame([row]).to_csv(OUT/"supporting/data/selected_condition.csv",index=False)
    plt.rcParams.update({"font.size":11,"pdf.fonttype":42,"svg.fonttype":"none"})
    fig=plt.figure(figsize=(12,8))
    fig.text(.5,.966,"Local scalar activity adds to a global vector",ha="center",fontsize=20,weight="bold")
    fig.text(.5,.922,"One actual narrow-grid experiment: Holmes C=30 · initial physical d₀=14 · stationary cue T=300",ha="center",fontsize=11)
    fig.text(.5,.868,r"$M_0=\sum_j u_j\Delta s,\qquad M_1=\sum_j u_j\Delta s\,\mathbf{n}_j,\qquad P=M_1/M_1(0)$",ha="center",fontsize=16)
    fig.text(.5,.827,"Complex numbers package two real Cartesian components. Initial normalization sets P(0)=1.",ha="center",fontsize=10.5)
    norm=Normalize(0,float(fields.max()));cmap=plt.get_cmap("viridis");scale=.83
    axes=[fig.add_axes((.02,.26,.46,.50)),fig.add_axes((.52,.26,.46,.50))]
    for k,ax in enumerate(axes):
        half=np.pi/n
        a=np.column_stack((np.cos(alpha-half),np.sin(alpha-half)));b=np.column_stack((np.cos(alpha+half),np.sin(alpha+half)))
        ax.add_collection(LineCollection(np.stack((a,b),axis=1),array=fields[k],cmap=cmap,norm=norm,linewidth=14,capstyle="butt",zorder=3))
        ax.plot([-1.4,1.4],[0,0],color=".9",lw=.6);ax.plot([0,0],[-1.4,1.4],color=".9",lw=.6)
        target=np.deg2rad(row.target_angle_in_P_coordinates_deg)
        angle=np.linspace(target-np.pi*row.full_physical_support/C,target+np.pi*row.full_physical_support/C,100)
        ax.plot(1.09*np.cos(angle),1.09*np.sin(angle),color="#cf852b",lw=3,alpha=.9)
        if k==1:
            ax.annotate("cue footprint",xy=(1.09*np.cos(target),1.09*np.sin(target)),xytext=(-.55,.55),color="#b47622",fontsize=8,ha="center",arrowprops=dict(arrowstyle="-",color="#b47622",lw=.8))
        for j,value in enumerate(sectors[k]):
            tail=1.12*np.exp(2j*np.pi*j/12);v=scale*value;length=abs(v)
            ax.arrow(tail.real,tail.imag,v.real,v.imag,width=min(.009,.18*length),head_width=min(.055,.55*length),head_length=min(.075,.45*length),length_includes_head=True,color="#327aa5",zorder=4)
        v=scale*P[k]
        ax.arrow(0,0,v.real,v.imag,width=.014,head_width=.065,head_length=.095,length_includes_head=True,color="#b93e50",zorder=6)
        ax.scatter([0],[0],s=10,color="#b93e50",zorder=7)
        ax.text(.03,-.25,"resultant P",ha="left",fontsize=10,color="#b93e50")
        ax.set_xlim(-1.52,1.52);ax.set_ylim(-1.52,1.52);ax.set_aspect("equal");ax.axis("off")
        ax.set_title("Initial native basal" if k==0 else "At cue withdrawal: two active regions",fontsize=12,weight="bold",pad=8)
        fig.text(.25+.5*k,.235,f"Active amount M₀/M₀(0) = {M0[k]/M0[0]:.4f}",ha="center",fontsize=12,color="#2d7864")
        fig.text(.25+.5*k,.192,f"Vector magnitude R = {R[k]:.4f}",ha="center",fontsize=12,color="#b93e50")
    cax=fig.add_axes((.38,.116,.24,.018));cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation="horizontal")
    cb.set_label("Native active concentration (shared scale)",fontsize=9,labelpad=2);cb.ax.tick_params(labelsize=8)
    fig.text(.025,.12,"Blue arrows: exact normalized sums\nin twelve 30° sectors. All arrows share one scale.",fontsize=8.5,va="top",color="#327aa5")
    fig.text(.67,.12,"A smaller resultant can accompany strong local activity.\nThis is an activity descriptor, without a movement law.",fontsize=8.5,va="top")
    fig.text(.025,.025,"Cue-on states only; the later finite field tail fails and is not used as stability evidence. No site ancestry or intrinsic remodeling-speed claim.",fontsize=9)
    for ext in ("png","pdf","svg"):
        p=OUT/"measured_activity_vector.png" if ext=="png" else OUT/"supporting/vectors"/("measured_activity_vector."+ext)
        fig.savefig(p,dpi=200,facecolor="white")
    plt.close(fig)
    stamp=datetime.now(timezone.utc).isoformat()
    (OUT/"README.md").write_text(f"""# Activity-vector construction from the current experiment

Generated {stamp}. [PNG](measured_activity_vector.png) · [PDF](supporting/vectors/measured_activity_vector.pdf) · [SVG](supporting/vectors/measured_activity_vector.svg) · [every physical native-cell contribution](supporting/data/native_cell_contributions.csv) · [exact sector sums](supporting/data/exact_sector_sums.csv) · [two measured states](supporting/data/two_state_measurements.csv).

Selected condition383 is Holmes C30, prescribed initial physical centre distance14 (target168°), cue duration300, with full support3 and native coefficient peak1. Its whole support is outside the initial10%-contrast region. These are its actual initial and withdrawal active profiles, without further integration. The second state contains {stats[1]['peaks']} resolved peaks and {stats[1]['halfcontrast_regions']} half-contrast regions. All colour scales are shared, zero-based native concentration scales; the circle is a schematic representation of the periodic domain, not a measured cell outline.

Each physical cell contributes its real scalar activity times physical arc length times outward unit normal, represented as a complex Cartesian pair. Exact twelve-sector sums of these native contributions are divided by the same initial complexM1. Translating sector-arrow tails for display does not change their vectors. All sector and resultant arrows use the common0.83vector scale. Normalized axes are rotated once by the initialM1 phase; raw physicalM1 is separately exported. Sector sums reproduce P with maximum error {error:.3g}. Total active amount changes to {M0[1]/M0[0]:.7f} of its initial value, whereas R is {R[1]:.7f}. Thus strong local activity and almost retained total active amount can yield a much smaller resultant because directional contributions partly cancel. This measures global activity organization, not peak travel, cell velocity or site ancestry.

Selection illustrates the vector definition using the SAME645-condition experiment and a clearly external near-antipodal placement, rather than adding another protocol. No general mechanism is identified from this example. These are cue-on states only; the later13,000-unit cue-free field tail fails the declared finite gate and is retained as such in the full outcome supplement, without an asymptotic-stability claim.
""")
    artifacts={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob("*") if p.is_file() and p.name!="verification.json"}
    result=dict(generated_at=stamp,passed=True,condition_id=383,states=2,native_cells=n,sector_count=12,physical_ds=ds,common_vector_scale=scale,
                exact_sector_sum_error=error,measurements=stats,new_integrations=0,source_sha256={str(path):sha(path),str(cue):sha(cue)},renderer_sha256=sha(__file__),artifact_sha256=artifacts)
    (OUT/"supporting/metadata/verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"output":str(OUT/"measured_activity_vector.png"),"measurements":stats,"sector_error":error},indent=2))


if __name__=="__main__":main()
