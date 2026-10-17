"""Render cohesive supported-ring evidence from audited retained/new measurements.

No native integration. Two main component sheets keep function/correspondence
and finite competition prominent. Full-grid moments and the Goryachev return
pair remain separate readable support rather than a gallery of unrelated cases.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,csv,hashlib,json,os,re,sys
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
sys.dont_write_bytecode=True
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,ListedColormap
from matplotlib.patches import Rectangle,Patch
from matplotlib.lines import Line2D
import length_dev_output_layout as layout

ROOT=Path(__file__).resolve().parents[2]
REPO=ROOT.parents[1]
TOPIC=ROOT/"figures/length_mechanisms_20261008/dev/vision_ring_20261017"
SOURCE=ROOT/"data/runs/chemical_ring_functional_v2_20261017"
MOMENTS=ROOT/"data/runs/vision_ring_moments_20261017"
RECOVERY=ROOT/"data/runs/vision_ring_recovery_20261017"
DURATION=ROOT/"data/runs/vision_ring_duration_20261017"
LINEAR=ROOT/"data/runs/vision_linear_20261017/recovery/results/recovery_summary.csv"
COMPONENT=ROOT/"figures/length_mechanisms_20261008/dev/expanded_ring_components_20261017"
MODELS=("wave_pinning","goryachev","holmes_model3")
NAMES=dict(wave_pinning="Wave-pinning",goryachev="Goryachev",holmes_model3="Holmes 3")
COLORS=("#4678a3","#d77c35","#8b65a9","#69a7a0","#bcb069")
LABELS=("Single peak near old; identity unresolved","Cue-associated peak","Distinct old + cue peaks","One region spans old + cue","Shifted / merged; unresolved")
EXAMPLES=(("1","holmes_model3",30,12),("2","holmes_model3",40,12),("3","holmes_model3",40,16),("4","goryachev",40,20),("5","goryachev",80,40))
ARTIFACTS=[];SOURCES={}


def polish(text):
    """Spacing for review prose, without touching code or link destinations."""
    pieces=re.split(r'(\]\([^)]+\)|```[\s\S]*?```|`[^`]*`|C=2L)',text)
    for i,s in enumerate(pieces):
        if s.startswith(('](', '```', '`')) or s=='C=2L':continue
        s=re.sub(r'(Holmes|Goryachev|WP)(?=C[0-9])',r'\1 ',s)
        s=re.sub(r'([A-Za-z]{2,}(?:-[A-Za-z]{2,})*)(?=[0-9])',r'\1 ',s)
        s=re.sub(r'(?<=[0-9])(?=[A-Za-z])',' ',s)
        s=re.sub(r'(?<=%)(?=[A-Za-z])',' ',s)
        s=re.sub(r'\bFigure ([0-9]+) ([A-Z])',r'Figure \1\2',s)
        s=re.sub(r'(?<=[0-9]) e(?=[+-][0-9])','e',s)
        s=re.sub(r'([0-9]{4}-[0-9]{2}-[0-9]{2}) T',r'\1T',s)
        s=re.sub(r'\bC ([0-9]+)',r'C\1',s)
        s=re.sub(r',(?=[A-Za-z0-9])',', ',s);s=re.sub(r';(?=\S)','; ',s)
        for a,b in (('belowR=.1','below R = 0.1'),('primaryR=.1','primary R = 0.1'),('R=.1','R = 0.1'),
                    ('R<.1','R < 0.1'),('R≥.1','R ≥ 0.1'),('gates.05and.2','gates 0.05 and 0.2'),
                    ('gates.05 and.2','gates 0.05 and 0.2'),('and.2','and 0.2'),('Blackdot','Black dot'),
                    ('BlackR','Black R'),('orangeangle','orange angle'),('darkgray','dark gray'),
                    ('nativeRac','native Rac'),('rawM','raw M'),('RawM','Raw M'),('h.05','h = 0.05'),
                    (' on/',' on / '),('twoC','two C'),('threeC','three C'),
                    ('unperturbedt ','unperturbed t = '),('R = .1','R = 0.1'),('newresponses','new responses')):s=s.replace(a,b)
        pieces[i]=s
    return ''.join(pieces)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda:stream.read(2**20),b""):h.update(chunk)
    return h.hexdigest()


def source(path,expected=None):
    path=Path(path).resolve();digest=sha(path)
    if expected is not None:assert digest==expected,path
    SOURCES[str(path)]=digest;return path


def read(path):
    with source(path).open(newline="") as stream:return list(csv.DictReader(stream))


def number(row,key,default=np.nan):
    try:return float(row[key]) if row[key] not in ("",None) else default
    except (KeyError,TypeError,ValueError):return default


def truth(row,key):return str(row.get(key,"")).lower() in ("true","1")
def key(row):return row["model"],number(row,"circumference"),number(row,"fraction")
def primary(rows):return [r for r in rows if number(r,"epsilon")==.005]


def csvwrite(name,rows):
    path=layout.dev_output_path(TOPIC,name)
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));writer.writeheader();writer.writerows(rows)
    ARTIFACTS.append(path)


def export(fig,stem,caption):
    for suffix in (".png",".pdf",".svg"):
        path=layout.dev_output_path(TOPIC,stem+suffix);fig.savefig(path,dpi=190 if suffix==".png" else None,bbox_inches="tight");ARTIFACTS.append(path)
    plt.close(fig)
    path=layout.dev_output_path(TOPIC,stem+"_caption.txt");path.write_text(polish(caption)+"\n");ARTIFACTS.append(path)


def category(text):
    if text.startswith("one old-core-associated"):return 0
    if text=="one peak associated with cue patch":return 1
    if text=="distinct peaks associated with old core and cue patch":return 2
    if text.startswith("one connected region spanning"):return 3
    if text.startswith("one shifted"):return 4
    raise ValueError(text)


def draw_map(ax,matrix,model,sizes,value,norm,cmap,lookups=None,regional=False,highlights=False):
    for iy,C in enumerate(sizes):
        rr=sorted([r for r in matrix if r["model"]==model and number(r,"circumference")==C],key=lambda r:number(r,"fraction"))
        xx=np.array([number(r,"realized_fraction",number(r,"fraction")) for r in rr]);edges=np.r_[xx[0]-.0125,(xx[:-1]+xx[1:])/2,xx[-1]+.0125]
        for j,row in enumerate(rr):
            valid=row["status"]=="complete";prep=row["status"]=="unqualified_preparation"
            rec=(lookups.get(key(row)) if lookups is not None else row) if valid else None
            v=value(rec) if valid and rec is not None else np.nan
            face=cmap(norm(v)) if valid and np.isfinite(v) else "white" if valid else "#777777" if prep else "#dedede"
            ax.add_patch(Rectangle((edges[j],iy-.5),edges[j+1]-edges[j],1,fc=face,ec="white",lw=.25))
            if not valid:ax.plot(xx[j],iy,"_" if prep else "x",color="white" if prep else "#888888",ms=2.8,mew=.6)
            elif not np.isfinite(v):ax.plot(xx[j],iy,"o",mfc="none",mec="#8c8c8c",ms=3,mew=.7)
            elif regional and not truth(rec,"distinct_old_peak_available"):ax.plot(xx[j],iy,".",color="black",ms=2.5)
            if highlights and valid:
                for label,m,c,d in EXAMPLES:
                    if m==model and c==C and np.isclose(number(row,"short_route"),d):
                        ax.text(xx[j],iy,label,ha="center",va="center",fontsize=7.5,fontweight="bold",bbox=dict(boxstyle="circle,pad=.16",fc="white",ec="black",lw=1))
    ax.set(xlim=(.0875,.5125),ylim=(-.5,len(sizes)-.5),xlabel="Realized shorter-route fraction",ylabel="Circumference")
    ax.set_xticks((.1,.2,.3,.4,.5));ax.set_yticks(range(len(sizes)),[f"{v:g}" for v in sizes]);ax.tick_params(labelsize=7.5)


def reference_profiles(ax,model,pairs,basal,common_pairs):
    row=next(r for r in pairs if r["model"]==model and number(r,"linear_L")==15)
    with np.load(source(row["linear_checkpoint"],row["linear_checkpoint_sha256"])) as z:x=z["x"].copy();u=z["state"][0].copy()
    with np.load(source(row["ring_checkpoint"],row["ring_checkpoint_sha256"])) as z:rx=z["x"].copy();ru=z["state"][0].copy()
    origin=rx[int(np.argmax(ru))];rx=(rx-origin+15)%30-15;order=np.argsort(rx)
    ax.plot(np.r_[-x[::-1],x],np.r_[u[::-1],u],"--",color="#242424",lw=1.5,label="Native interval reflected")
    ax.plot(rx[order],ru[order],color="#cf7b33",lw=1.15,label="Native C30 ring")
    width=number(row,"ring_full_core_width");span=min(15,width*.8)
    ax.set(xlim=(-span,span),xlabel="Distance from original peak (native units)",ylabel="Primary activity (native units)",title=NAMES[model])
    b=next(r for r in basal if r["model"]==model and number(r,"circumference")==30)
    l=next(r for r in common_pairs if r["model"]==model and number(r,"linear_length")==15)
    text=f"Full width ratio {number(row,'ring_full_core_width')/number(row,'linear_core_width'):.4g}; excess amount ratio {number(row,'ring_full_excess_core_amount')/number(row,'linear_excess_core_amount'):.4g}"
    text+=f"\nMatched support3→6 gain ratio {number(row,'matched_cue_weighted_gain_ring_over_linear'):.6g}"
    assert np.isclose(number(l,"tau"),number(b,"tau"),rtol=0,atol=1e-14)
    text+=f"\nCommon active-loss F₀ ring/line {number(l,'ring_over_linear_K'):.5g}"
    ax.text(.02,.98,text,transform=ax.transAxes,va="top",fontsize=7.3,bbox=dict(fc="white",ec="none",alpha=.85))
    ax.tick_params(labelsize=8)


def function_sheet(matrix,pairs,basal,challenged,common_pairs,sizes):
    fig=plt.figure(figsize=(13.6,13.8));gs=fig.add_gridspec(3,3,height_ratios=(.8,1.45,1.45),left=.055,right=.94,bottom=.07,top=.94,hspace=.34,wspace=.29)
    for col,model in enumerate(MODELS):reference_profiles(fig.add_subplot(gs[0,col]),model,pairs,basal,common_pairs)
    fig.axes[0].legend(frameon=False,fontsize=7,loc="lower left")
    gain_refs={m:number(next(r for r in matrix if r["model"]==m and number(r,"circumference")==30 and number(r,"fraction")==.5),"initial_cue_weighted_gain") for m in MODELS}
    lookup={key(r):r for r in challenged}
    basal_lookup={(r['model'],number(r,'circumference')):number(r,'F') for r in basal}
    for col,model in enumerate(MODELS):
        ax=fig.add_subplot(gs[1,col]);draw_map(ax,matrix,model,sizes,lambda r:number(r,"initial_cue_weighted_gain")/gain_refs[model],Normalize(0,3),plt.get_cmap("viridis"));ax.set_title(f"{NAMES[model]}: initial gain / own C30 reference",fontsize=10)
        ax=fig.add_subplot(gs[2,col]);draw_map(ax,matrix,model,sizes,lambda r:number(r,"F")/basal_lookup[r['model'],number(r,'circumference')] if truth(r,"passes") else np.nan,Normalize(0,1.1),plt.get_cmap("viridis"),lookups=lookup,regional=True)
        tau=number(next(r for r in basal if r["model"]==model),"tau")
        ax.set_title(f"t10 recovery / same-size basal; τ = {tau:.5g}",fontsize=10)
    for y,height,norm,text in ((.40,.25,Normalize(0,3),"Onset gain / own C30 antipodal"),(.08,.25,Normalize(0,1.1),"F(t10) / F₀(C), common active-loss assay")):
        cax=fig.add_axes((.957,y,.012,height));cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap="viridis"),cax=cax);cb.set_label(text,fontsize=9);cb.ax.tick_params(labelsize=8)
    fig.suptitle("Ring counterparts and local responses: native correspondence, input gain and active-loss recovery",fontsize=15,y=.987)
    fig.text(.055,.015,"Profiles use exact L15/C30 native states; full 13-pair evidence remains linked. Gain maps retain full support3, whereas matched baseline input is endpoint3→ring6.\nAll maps retain867planned cells. Gray×: overlap excluded; darkgray: HolmesC80 unavailable. Blackdot: no resolved distinct old site; response concerns fixed original core.\nRecovery removes0.5% of the frozen basal active-excess direction, transfers it to inactive, evolves all fields and compares with sham. Clocks differ across native models.",fontsize=8.5)
    export(fig,"vision_ring_function_correspondence","Three exact L15/C30 native profile overlays demonstrate the C=2L relationship without interpolating new states; display windows focus on the incumbent and retain native position/activity units. Full ring width and excess-active amount double, while matched cue-weighted native input gain uses endpoint support3 and reflected ring full support6. The actual complete-grid challenges retain full support3. The initial-gain map uses each model's own C30 antipodal reference and a shared0–3 scale so nearly null distance effects remain null. The recovery map is the new active-incumbent deactivation assay only: primary ε=.005 of a frozen cosine-tapered basal excess-active direction is transferred pointwise to inactive; all native fields evolve with the actual cue on from the exact t10 snapshot. F=1−D(τ)/D(0), where D is the fixed-original-core sham-minus-probe amount. The displayed F(t10)/F0(C) uses each size's exact basal control and a common0–1.1 scale. Clocks are fixed within each model but not cross-calibrated. All702states pass the prefix numerical/amplitude/saturation gates. Raw basal F/K size curves are immediately linked so conditional normalization cannot hide basal size effects. Blackdots mark regional responses without a resolved distinct old site. Positive deficit recovery does not imply absolute core maintenance or generic stability.")


def trace_for(summary,flat):
    index=int(number(summary,"trace_index"));a,b=flat["offsets"][index:index+2]
    return {name:flat[name][a:b] for name in flat.files if name not in ("offsets","matrix_row_index")}


def scalar_vector_plot(ax,trace):
    t=trace["time"];angle=np.angle(trace["P_x"]+1j*trace["P_y"])*180/np.pi;angle[trace["R"]<.1]=np.nan
    ax.plot(t,trace["R"],color="#202020",lw=1.35,label="Magnitude R")
    right=ax.twinx();right.plot(t,angle,color="#c57530",lw=1.1,label="Angle, R ≥ 0.1")
    ax.axvline(1000,ls="--",color="gray",lw=.7);ax.axhline(.1,color="gray",ls=":",lw=.6)
    ax.set(xlim=(0,14000),ylim=(0,1.55),xticks=(0,1000,4000,9000,14000),xlabel="Time from onset (native units)",ylabel="R = |P|")
    right.set(ylim=(-190,190),yticks=(-180,-90,0,90,180),ylabel="Activity-vector angle (°)")
    ax.tick_params(labelsize=7.5);right.tick_params(labelsize=7.5,colors="#a96228");right.yaxis.label.set_color("#a96228")
    return right


def native_kymo(ax,time,x,u,C,d,vmax):
    mesh=ax.pcolormesh(np.r_[x,x[0]+C],time,np.c_[u,u[:,0]],shading="nearest",cmap="magma",vmin=0,vmax=vmax,rasterized=True)
    ax.set_facecolor("#bdbdbd");ax.hlines(1000,-10,C-10,ls="--",lw=.75,color="#83d6e2");ax.axvline(d,ls=":",lw=.6,color="#83d6e2")
    if C==30:ax.axvline(20,color="black",lw=.7);ax.text(25,2200,"outside C30",rotation=90,ha="center",va="center",fontsize=8,color="#555")
    ax.set(xlim=(-10,30),ylim=(0,4000),xticks=(-10,0,10,20,30),yticks=(0,1000,2000,3000,4000),xlabel="Position relative to original peak",ylabel="Time from cue onset")
    ax.tick_params(labelsize=7.5);return mesh


def finite_sheet(matrix,summaries,flat,sizes):
    fig=plt.figure(figsize=(13.6,14));gs=fig.add_gridspec(4,3,height_ratios=(1.4,.22,1.15,.95),left=.055,right=.925,bottom=.07,top=.94,hspace=.27,wspace=.35)
    for col,model in enumerate(MODELS):
        ax=fig.add_subplot(gs[0,col]);draw_map(ax,matrix,model,sizes,lambda r:category(r["final_location_association"]),Normalize(-.5,4.5),ListedColormap(COLORS),highlights=True);ax.set_title(NAMES[model]+": local endpoint, t14000",fontsize=11)
    legend=fig.add_subplot(gs[1,:]);legend.axis("off");legend.legend(handles=[Patch(fc=c,label=l) for c,l in zip(COLORS,LABELS)]+[Line2D([],[],marker="x",ls="",color="#888",label="Overlap excluded"),Patch(fc="#777",label="HolmesC80 unavailable")],loc="center",ncol=3,frameon=False,fontsize=8)
    path=source(COMPONENT/"supporting/data/expanded_ring_Holmes_native_display_samples.npz")
    with np.load(path) as fields:
        for col,(label,model,C,d) in enumerate(EXAMPLES[:3]):
            ax=fig.add_subplot(gs[2,col]);mesh=native_kymo(ax,fields[f"example{label}_time"],fields[f"example{label}_x"],fields[f"example{label}_primary"],C,d,2.2962084284844333)
            ax.set_title(f"{label}  Holmes C{C}, d{d}: "+("cue single" if label=="1" else "merged / bridged" if label=="2" else "two sites"),fontsize=11)
            row=next(r for r in summaries if r["model"]==model and number(r,"circumference")==C and np.isclose(number(r,"short_route"),d))
            ax=fig.add_subplot(gs[3,col]);scalar_vector_plot(ax,trace_for(row,flat));ax.set_title("Same native history: global activity response",fontsize=10)
    cax=fig.add_axes((.947,.335,.012,.235));cb=fig.colorbar(mesh,cax=cax);cb.set_label("Native active Rac",fontsize=9)
    fig.suptitle("Substantial finite ring challenge: local patterns and the activity vector answer different questions",fontsize=14.3,y=.987)
    fig.text(.055,.012,"Maps: 1000on /13000off, final t14000. Histories: exact h.05 and linear0–4000; nativeRac shared unclipped0–2.29621. Original site is intact at0 under periodic cut−10.\n1→2 changes size/material/architecture/complementary route;2→3 relocates cue in the exact same complete C40 state.5/4 mark the linked Goryachev return/replacement pair.\nP=M₁/M₁(0) is an activity descriptor. BlackR and orangeangle retain magnitude and direction; angle is withheld forR<.1. No cell velocity or deterministic persistence is inferred.",fontsize=8.5)
    export(fig,"vision_ring_finite_fates_activity_vector","Complete three-model geographic endpoint context under the native width3 1000-on/13000-off challenge, followed by the three numbered exact Holmes histories and their same-source activity-vector responses. The maps retain all867planned cells and numerical/geometry availability. Endpoint geography at14000 is finite-window evidence, not stable-branch assignment. Native h.05 histories C30d12, C40d12 and C40d16 are displayed with a common unclipped Rac scale and common linear0–4000time; periodic reindexing keeps the original region intact at0. C30 has gray20..30 outside its shorter physical ring. The first comparison changes size and both routes/material/architecture; the two C40 cases share identical complete basal Rac/Rho fields and differ only in cue location. Global P=M1/M1(0) uses the unchanged reviewed first-moment definition; R=|P| is plotted alongside masked angle through the full14000time window. Direction is omitted belowR=.1, with .05/.2 sensitivities linked. A merged region, a unique cue-associated survivor and two distinct sites can produce different or ambiguous moment responses, so local histories remain essential. Circled4/5 identify linked Goryachev replacement/return support. No activity moment is identified with movement or migration persistence.")


def moment_maps(matrix,summaries,sizes):
    lookup={key(r):r for r in summaries if r["status"]=="complete"}
    fig,axes=plt.subplots(2,3,figsize=(13.5,10.2));fig.subplots_adjust(left=.055,right=.92,bottom=.13,top=.92,hspace=.25,wspace=.28)
    for col,model in enumerate(MODELS):
        draw_map(axes[0,col],matrix,model,sizes,lambda r:number(r,"R_t14000"),Normalize(0,1.5),plt.get_cmap("viridis"),lookups=lookup,highlights=True)
        axes[0,col].set_title(NAMES[model]+": endpoint R",fontsize=11)
        draw_map(axes[1,col],matrix,model,sizes,lambda r:number(r,"angle_deg_t14000_gate_0p1"),Normalize(-180,180),plt.get_cmap("twilight_shifted"),lookups=lookup,highlights=True)
        axes[1,col].set_title("Endpoint angle; unavailable for R < 0.1",fontsize=10)
    for y,norm,cmap,title in ((.56,Normalize(0,1.5),"viridis","R at t14000"),(.16,Normalize(-180,180),"twilight_shifted","Masked activity angle (°)")):
        cax=fig.add_axes((.943,y,.014,.30));cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax);cb.set_label(title,fontsize=9)
    fig.suptitle("Complete-grid deterministic activity moments after 13000 cue-free units",fontsize=15,y=.985)
    fig.text(.055,.035,"Circled1–3: same three Holmes histories;4–5: Goryachev C40d20/C80d40. White open circle: angle masked by low magnitude; gray×: excluded; darkgray: unavailable.\nR is normalized to the initial first moment, not total activity. RawM₀,q₁,q₂ and sensitivity gatesR=.05/.2 are retained. Finite readout does not imply universal settling.",fontsize=9)
    export(fig,"vision_ring_full_grid_moment_maps","All652complete retained conditions, with all867planned cells kept in context, use P=M1/M1(0) and R=|P| from the unchanged reviewed quadrature. Endpoint time14000 means13000cue-free units after the1000-unit challenge. Shared magnitude range0–1.5 contains all retained readout magnitudes; cyclic angle range−180..180 is masked whereverR<.1. White circles distinguish masked direction from gray geometric exclusions and unavailable preparations. Numbered cells match selected native histories. Uniform activity cancels from M1/M2 numerators; mass-normalized q1/q2 and total activity M0 remain diagnostics. R and angle do not identify a unique site, stable state, velocity or migration direction. Readout-time sensitivity remains linked because some WP/Holmes conditions change late.")


def basal_support(basal,challenged,sizes):
    fig,axes=plt.subplots(2,3,figsize=(12.8,6.4));rows=[]
    for col,model in enumerate(MODELS):
        rr=sorted([r for r in basal if r["model"]==model],key=lambda r:number(r,"circumference"))
        ax=axes[0,col]
        ax.plot([number(r,"circumference") for r in rr],[number(r,"F") for r in rr],"o-",color="#426b80",ms=3)
        ax.set(ylim=(0,.2),xlabel="Circumference",ylabel="Basal recovered deficit fraction F₀",title=NAMES[model]);ax.tick_params(labelsize=8)
        ax=axes[1,col];ax.plot([number(r,"circumference") for r in rr],[number(r,"K") for r in rr],"o-",color="#426b80",ms=3)
        ratio=number(next(r for r in rr if number(r,'circumference')==39),'K')/number(next(r for r in rr if number(r,'circumference')==19),'K')
        ax.set(xlabel="Circumference",ylabel="Basal K₀ = F₀ / τ (native time⁻¹)",title=f"C19→39: {ratio:.5g}×");ax.tick_params(labelsize=8)
        for r in rr:rows.append({k:r[k] for k in ("model","circumference","tau","F","K","D0","Dtau","E","passes")})
    fig.suptitle("Basal active-loss recovery is measured once per model/circumference",y=1.02);fig.tight_layout()
    export(fig,"vision_ring_basal_active_recovery","New common active-loss assay on all50qualified exact native basals, primaryε=.005, one fixed clock within each model. F0=1−Dτ/D0; all native fields evolve and explicit inactive pools receive removed active material pointwise. The top row uses a shared0–.2 ordinate for absolute deficit fractions. The lower row retains raw K0=F0/τ with explicitly different native axes and prints C19→39 folds so narrow axes do not imply equal size effects. Model clocks differ and do not calibrate rates across models. A basal readout has no future cue-distance dependence. These measurements import no inactive-withdrawal coefficients.")
    csvwrite("vision_ring_basal_active_recovery.csv",rows)


def readout_tables(readouts):
    lookup={(r["model"],number(r,"circumference"),number(r,"fraction"),number(r,"time_from_onset")):r for r in readouts}
    comparisons=[];table=[]
    for row in readouts:
        t=number(row,"time_from_onset");end=lookup[row["model"],number(row,"circumference"),number(row,"fraction"),14000]
        a=complex(number(row,"P_x"),number(row,"P_y"));b=complex(number(end,"P_x"),number(end,"P_y"))
        valid=truth(row,"angle_available_gate_0p1") and truth(end,"angle_available_gate_0p1")
        comparisons.append(dict(model=row["model"],circumference=number(row,"circumference"),fraction=number(row,"fraction"),
            time_from_onset=t,time_after_withdrawal=number(row,"time_after_withdrawal"),
            absolute_R_difference_from_t14000=abs(number(row,"R")-number(end,"R")),
            wrapped_angle_difference_from_t14000_deg=abs(np.angle(a/b))*180/np.pi if valid else "",
            angle_available_at_both=valid,angle_availability_changed=truth(row,"angle_available_gate_0p1")!=truth(end,"angle_available_gate_0p1")))
    for model in MODELS:
        for t in (1000,4000,9000):
            rr=[r for r in comparisons if r["model"]==model and r["time_from_onset"]==t]
            largest=max(rr,key=lambda r:r["absolute_R_difference_from_t14000"])
            aa=[number(r,"wrapped_angle_difference_from_t14000_deg") for r in rr if r["angle_available_at_both"]]
            table.append(dict(model=model,time_from_onset=t,post_withdrawal_time=t-1000,
                max_R_difference=largest["absolute_R_difference_from_t14000"],median_R_difference=float(np.median([r["absolute_R_difference_from_t14000"] for r in rr])),
                max_R_case_C=largest["circumference"],max_R_case_fraction=largest["fraction"],max_wrapped_angle_difference=max(aa),
                angle_available_at_both=sum(r["angle_available_at_both"] for r in rr),angle_availability_changes=sum(r["angle_availability_changed"] for r in rr)))
    csvwrite("vision_ring_readout_comparisons.csv",comparisons);csvwrite("vision_ring_readout_sensitivity_summary.csv",table)
    return table


def goryachev_support(matrix,summaries,flat):
    fig,axes=plt.subplots(2,2,figsize=(11.8,8.5),gridspec_kw={"height_ratios":(1.5,1.)})
    histories=[]
    for C,d in ((40,20),(80,40)):
        row=next(r for r in matrix if r['model']=='goryachev' and number(r,'circumference')==C and number(r,'short_route')==d)
        phases=[]
        for phase,offset in (('cue',0),('recovery',1000)):
            with np.load(source(row[phase+'_path'])) as z:
                phases.append((z['time'].copy()+offset,z['x'].copy(),z['states'][:,0].copy() if 'states' in z else z['primary'].copy()))
        np.testing.assert_array_equal(phases[0][2][-1],phases[1][2][0])
        t=np.r_[phases[0][0],phases[1][0][1:]];u=np.r_[phases[0][2],phases[1][2][1:]];x=phases[0][1]
        origin=x[np.argmax(u[0])];px=(x-origin+10)%C-10;order=np.argsort(px);choose=t<=4000
        histories.append((C,d,t[choose],px[order],u[choose][:,order],float(u.max())))
    vmax=max(h[-1] for h in histories)
    for col,(C,d,t,x,u,_) in enumerate(histories):
        ax=axes[0,col];mesh=ax.pcolormesh(np.r_[x,x[0]+C],t,np.c_[u,u[:,0]],cmap='magma',vmin=0,vmax=vmax,shading='nearest',rasterized=True)
        ax.set_facecolor('#bdbdbd');ax.hlines(1000,-10,C-10,color='#83d6e2',ls='--',lw=.8)
        if C==40:ax.axvline(30,color='black',lw=.7);ax.text(50,2200,'outside C40 ring',rotation=90,ha='center',va='center',fontsize=9,color='#555')
        ax.set(xlim=(-10,70),ylim=(0,4000),xlabel='Position relative to original peak',ylabel='Time from cue onset',title=f'Goryachev C{C}, antipodal d{d}')
        row=next(r for r in summaries if r['model']=='goryachev' and number(r,'circumference')==C and number(r,'fraction')==.5)
        scalar_vector_plot(axes[1,col],trace_for(row,flat))
    fig.subplots_adjust(left=.065,right=.91,bottom=.17,top=.88,hspace=.33,wspace=.35)
    cax=fig.add_axes((.938,.51,.015,.33));fig.colorbar(mesh,cax=cax,label='Native active field')
    fig.suptitle('Geographically different single-peak endpoints: replacement and retained old activity',fontsize=14,y=.98)
    fig.text(.065,.035,polish('Both cues are antipodal and use the same native input law, local width and1000on/13000off protocol. Size, material, architecture and both routes change together.\nThe C80 full history retains original-region activity while the challenger transient vanishes; C40 relocates activity toward the cue. No pure distance mechanism is isolated.\nTop: linear0–4000 display, exact native samples, periodic cut−10. Bottom: same-source activity magnitude/angle through14000; direction withheld forR<.1.'),fontsize=8.5)
    export(fig,'vision_ring_Goryachev_return_support','Selected full native Goryachev C40d20 and C80d40 histories, highlighted as4/5 on the complete source maps. The comparison keeps antipodal fraction, local width3, peak6.4 and1000on/13000off protocol, but changes circumference, conserved material, architecture and both absolute routes jointly. Periodic reindexing keeps the original region intact at0; physical position0 is not duplicated as two sites. Shared activity scale is unclipped. C40 relocates toward the cue; C80 retains activity near the original region throughout and the challenged transient vanishes after withdrawal. The same-source R and masked angle distinguish these activity outcomes without assuming movement. This direct history supports the bounded C80 retention interpretation beyond its sole-old-geographic endpoint label. It does not reclassify every sole-old endpoint as a tracked incumbent or stable state.')


def markdown_table(headers,rows):
    return '| '+' | '.join(headers)+' |\n| '+' | '.join('---' for _ in headers)+' |\n'+'\n'.join('| '+' | '.join(str(x) for x in row)+' |' for row in rows)


def fmt(value):return f'{value:.5g}' if np.isfinite(value) else 'masked'


def methods():
    text=r'''# Activity-moment and active-loss definitions

The global extraction imports the reviewed `moment()` and `analytic_checks()` unchanged from `global_polarity_vector_20261017.py`, with a preserved copy in the new extraction run. For primary active field values at native periodic cell centres,

$$
M_n(t)=h\sum_i u_i(t)\exp(\mathrm{i}n\theta_i),\qquad
\theta_i=\frac{2\pi(x_i-x_0)}{C},\qquad h=\frac{C}{N}.
\tag{1}
$$

The activity vector, its magnitude and the mass-normalized diagnostics are

$$
P(t)=\frac{M_1(t)}{M_1(0)},\qquad R(t)=|P(t)|,\qquad
q_1(t)=\frac{|M_1(t)|}{M_0(t)},\qquad q_2(t)=\frac{|M_2(t)|}{M_0(t)}.
\tag{2}
$$

Complex division puts the initial incumbent at $P(0)=1$. Primary angle is $\arg P$, withheld when $R<0.1$; sensitivities use0.05 and0.2 without changing the vector. Uniform activity cancels from the first/second nonzero angular moment numerators. It changes $M_0$ and therefore mass-normalized $q_1,q_2$, which are retained alongside the initial-moment normalization. A small first moment can arise from cancellation between strong opposed sites. Neither a nonzero vector nor its angle uniquely identifies a protrusion, tracked site, velocity or direction of migration.

The new common active-loss perturbation uses the frozen basal half-contrast core, background and cosine taper. The removed primary profile is $\varepsilon$ times taper times basal excess above background, with primary$\varepsilon=0.005$. In these three native chemical models the same material is transferred pointwise to the explicit inactive pool. All native fields evolve under the same actual cue as their sham. The perturbation direction stays the original basal direction even at the challenged t10 snapshot; there is no renewed inactive-withdrawal/return step.

$$
D(t)=A_{\mathrm{sham}}(t)-A_{\mathrm{probe}}(t),\qquad
F=1-\frac{D(\tau)}{D(0)},\qquad K=\frac{F}{\tau}.
\tag{3}
$$

All amounts in this assay integrate the fixed original basal core. It does not follow a translating/broadening site. Fixed clocks are WP/Holmes0.1778279410038923 and Goryachev0.005623413251903491 native units, identical to the matched linear assay and fixed across sizes/placements within each model. They are not a calibration across models. The main map displays $F_{10}(C,d)/F_0(C)$ with each exact same-size basal control. Raw$F,K,D(0),D(\tau)$, errors, whole-active and outside-core deficits, actual sham/probe changes and rates remain linked.

All702 ring states pass the numerical, amplitude and prefix-saturation gates;7020 native trials were independently audited. Positive paired-deficit damping is conditional recovery along this disturbance, not actual net maintenance or generic stability. A black dot flags the eight t10 states lacking a resolved distinct old-site association; those remain regional responses. The independent linear/ring13-pair audit concerns independently prepared states with different native meshes, so close agreement is not exact mirrored-state identity or agreement within numerical$E$.

The complete original matrix has867planned cells:652complete histories,198geometry exclusions and17unqualified HolmesC80 preparations. No excluded state receives a calculated vector. Phase-local recovery time is shifted by1000 before joining; the common withdrawal state is retained once. The four readouts are onset1000,4000,9000,14000, equivalent to0,3000,8000,13000after withdrawal, all exactly present in native samples. Minima and maximum unmasked angles are sampled summaries, not interpolated event times. Later source sampling differs (5 or50native units); no greater temporal precision is inferred.

The lean NPZ stores concatenated native-time traces and offsets by complete condition: $M_0$, real/imaginary$M_1,M_2$, $P_x,P_y,R,q_1,q_2$. The867-row matrix preserves unavailable statuses;2608readout rows retain the primary and sensitivity angle gates. Native sources are referenced by immutable path/hash. All extraction is read-only; no solver or deterministic persistence fit is involved.
'''
    path=layout.dev_output_path(TOPIC,'MOMENT_AND_RECOVERY_METHODS.md');path.write_text(polish(text)+'\n');ARTIFACTS.append(path)


def write_results(matrix,pairs,common_pairs,basal,challenged,summaries,readouts,readout_summary,durations):
    function_rows=[];correspondence=[];fate_rows=[];vector_rows=[];duration_rows=[];sensitivity=[]
    for model in MODELS:
        bb=[r for r in basal if r['model']==model];cc=[r for r in challenged if r['model']==model]
        b19=next(r for r in bb if number(r,'circumference')==19);b39=next(r for r in bb if number(r,'circumference')==39)
        base={number(r,'circumference'):number(r,'F') for r in bb}
        ratios=[number(r,'F')/base[number(r,'circumference')] for r in cc]
        function_rows.append((NAMES[model],len(bb),fmt(number(b19,'K')),fmt(number(b39,'K')),fmt(number(b39,'K')/number(b19,'K')),f'{min(ratios):.5g}–{max(ratios):.5g}',sum(number(r,'sham_core_change')<0 for r in cc)))
        matches=[r for r in common_pairs if r['model']==model]
        correspondence.append((NAMES[model],len(matches),f"{100*max(abs(number(r,'relative_K_difference')) for r in matches):.5g}%"))
        rr=[r for r in matrix if r['model']==model and r['status']=='complete'];counts=[sum(category(r['final_location_association'])==k for r in rr) for k in range(5)]
        fate_rows.append((NAMES[model],len(rr),*counts))
    for label,model,C,d in EXAMPLES:
        r=next(r for r in summaries if r['model']==model and number(r,'circumference')==C and number(r,'short_route')==d)
        vector_rows.append((label,NAMES[model],C,d,fmt(number(r,'sampled_minimum_R')),fmt(number(r,'R_t14000')),fmt(number(r,'angle_deg_t14000_gate_0p1')),fmt(number(r,'q1_t14000')),fmt(number(r,'q2_t14000'))))
    keys=sorted({(r['model'],number(r,'circumference'),number(r,'fraction')) for r in durations})
    for model,C,fraction in keys:
        rr=[r for r in durations if r['model']==model and number(r,'circumference')==C and number(r,'fraction')==fraction and number(r,'cue_free_time')==13000]
        rr.sort(key=lambda r:number(r,'cue_duration'));assert len(rr)==3
        classes=[r['location_association'] for r in rr];assert len(set(classes))==1
        duration_rows.append((NAMES[model],f'{C:g}',fmt(C*fraction),LABELS[category(classes[0])],*(fmt(number(r,'R')) for r in rr),*(fmt(number(r,'angle_deg')) if number(r,'R')>=.1 else 'masked' for r in rr)))
    for r in readout_summary:
        if r['time_from_onset'] in (4000,9000):
            sensitivity.append((NAMES[r['model']],int(r['post_withdrawal_time']),fmt(r['max_R_difference']),fmt(r['median_R_difference']),f"C{r['max_R_case_C']:g}, f={r['max_R_case_fraction']:g}",fmt(r['max_wrapped_angle_difference'])))
    links=lambda p:os.path.relpath(p,TOPIC).replace(' ','%20')
    text=f'''# Supported rings: local response, finite competition and global activity

Generated {datetime.now(timezone.utc).isoformat()}; all retained-data producers and the new ring recovery/duration campaigns are complete. This is one supported-ring evidence package for the three-figure checkpoint, not a new main-text variant. It follows [MANUSCRIPT_VISION.md]({links(REPO/'MANUSCRIPT_VISION.md')}) and keeps noisy persistence experiments separate.

The ring counterparts preserve native interval organization and function at C=2L, while finite competition depends on both size and cue position. The new conditional active-loss response is weakened by a challenge relative to each exact same-size basal control, yet remains positive while every early unperturbed fixed core loses active material. The substantial1000-unit cue then produces cue-associated single peaks, merged regions or distinct sites. Those local patterns and the global activity vector provide complementary readouts; neither scalar recovery nor vector angle alone determines site fate or cell movement.

![Ring correspondence, input gain and the same new active-loss recovery assay](vision_ring_function_correspondence.png)

The upper native overlays use exact L15/C30 states, with native mesh samples and local display windows. Full ring width and excess-active core amount are twice the interval values; the same Figure1B cue-weighted derivative uses endpoint support3 versus reflected ring full support6. The actual ring maps retain full support3. All13matching native pairs and source hashes remain in the [baseline data]({links(COMPONENT/'supporting/data/expanded_ring_baseline_pair_measurements.csv')}). The new common active-loss13-pair comparison, at the identical native clocks, also agrees closely:

{markdown_table(('Model','Native C=2L pairs','Maximum relative K difference'),correspondence)}

These are independently prepared states and unmatched native meshes. Goryachev's largest difference occurs at L35/C70 and exceeds its paired numerical error; the agreement is not exact identity, a grid-convergence result or a biological calibration of circumference. The [new same-assay pair audit]({links(RECOVERY/'results/ring_linear_common_recovery_pairs_audit.json')}) is distinct from the exact mirrored-operator fixture.

The gain maps retain all867planned cells and one shared0–3 scale after each model's own C30 antipodal normalization. Distance dependence is nearly null: the largest within-size range is0.0392%WP,0.0891%Goryachev and0.00208%Holmes. The recovery row uses only the new common active-deactivation assay, displayed as F₁₀(C,d)/F₀(C) on a shared0–1.1scale. [Raw basal F/K curves](vision_ring_basal_active_recovery.png) are immediate support: normalizing the challenged map must not hide basal size effects.

{markdown_table(('Model','Basals','K₀ at C19','K₀ at C39','C39/C19','F₁₀/F₀(C) range','Early sham cores declining'),function_rows)}

Each recovery probe transfers0.5% of a frozen tapered basal excess-active direction from active to inactive, evolves all native fields and compares against the paired sham. Challenged probes start exactly10native units after cue onset, with the cue continuing throughout both sham and probe. All702states pass the numerical/amplitude/prefix gates;7020native trials were independently audited. The main statistic is the fraction of the initial sham-relative fixed-core deficit reduced by the fixed clock; its definition and all raw rates are in the [methods](supporting/metadata/MOMENT_AND_RECOVERY_METHODS.md) and [independent audit]({links(RECOVERY/'results/independent_data_audit.json')}). Positive conditional recovery is a restoring contribution along this specified disturbance. It does not mean the open original core gains material, that a moving/broadening site loses all its material, or that an incumbent is generically stable. Eight early states lack a resolved distinct old-site association and are flagged as regional responses.

The [fresh audit of the new recovery archives](supporting/metadata/FRESH_RESULT_AUDIT.json) confirms the actual loss allocation: at the exact unperturbedt10state, fixed-core diffusion is negative in all652conditions. Reaction is positive in500but smaller than outward diffusion, and negative in152. This allocates the native rate in the fixed original region. The diffusion term describes redistribution across the frozen core boundary; the regional total does not measure the whole moving/broadening incumbent. No legacy recovery coefficient enters this budget. Over the probe window all652shams lose fixed-core amount, while311perturbed probes lose amount and341rebound from their own kicked starting value. Conditional deficit recovery therefore does not imply a universal sign of the actual perturbed response.

**Corresponding measurements:** image native-equivalent active-reporter concentration, half-contrast width and excess integrated signal; measure the immediate response to a standardized local input and paired recovery from a small incumbent-targeted reduction. Report local perturbation/support, removed active amount, control trajectory and observation window. These measurements test the predicted response patterns without uniquely establishing the underlying model mechanism.

![Complete finite-competition context, same-source Holmes histories and activity vectors](vision_ring_finite_fates_activity_vector.png)

This sheet gives finite local competition equal prominence with transient assays. Top maps show the complete local geographic classification after1000on/13000off. They retain198overlap exclusions and17unavailable HolmesC80 preparations. The three exact native Holmes histories supply temporal evidence beyond endpoint geography: C30d12 develops a cue-associated single survivor; C40d12 merges/bridges; C40d16 retains distinct associated sites. The C30/C40 comparison jointly changes total material, architecture and complementary route18→28. The twoC40conditions relocate the cue within the identical complete basal Rac/Rho state.

{markdown_table(('Model','Complete','Sole near old, unresolved','Cue-associated','Distinct old+cue','Spanning / bridged','Shifted / unresolved'),fate_rows)}

The linked [Goryachev C40d20/C80d40 native/vector pair](vision_ring_Goryachev_return_support.png) exposes replacement versus retained original-region activity under the same antipodal fraction. The C80history keeps the old region active while the challenger transient vanishes; endpoint geography alone would not establish that history. This pair changes size/material/architecture and both absolute routes together and does not isolate a pure distance mechanism.

The same-source activity vector is P=M₁/M₁(0); R=|P| is retained alongside direction, with angle unavailable belowR=.1. Uniform activity cancels from the moment numerator. Raw total activity M₀ and mass-normalized q₁,q₂ remain diagnostics; q₂ helps expose opposed-site structure when the first moment is small but is not a unique site-count estimator. [Complete-grid R and masked-angle maps](vision_ring_full_grid_moment_maps.png) immediately complement the local outcome maps rather than replacing them.

At the endpoint, angle is masked in5/652conditions at the primaryR=.1 gate; sensitivity gates.05and.2 mask3and18conditions, respectively. The vector itself is unchanged by the gate. Magnitude and native local patterns are therefore essential when interpreting unavailable or threshold-dependent direction.

The [previously reviewed opposed-site cancellation component](../global_polarity_vector_20261017/global_polarity_component.png) directly displays the HolmesC40d20diagnostic: strong local sites can coexist with a near-zero first moment and unavailable angle. This existing support complements the present same-source examples without adding another illustrative campaign.

{markdown_table(('Map label','Model','C','d','Sampled min R','Endpoint R','Masked endpoint angle°','Endpoint q₁','Endpoint q₂'),vector_rows)}

**Corresponding measurements:** following a declared spatial input, image local site creation, displacement/merger and loss together with the two-dimensional reporter activity moment, magnitude and direction. Pair activity imaging with local protrusion measurements to assess their relationship; no velocity law is assumed. In particular, two strong opposed sites can make direction unavailable through cancellation. These deterministic challenges do not estimate migration direction or noisy polarity persistence; the separate noise figure addresses replicated stochastic persistence.

## Why this finite protocol, and what changes with readout

The width3 raised-cosine input is a **substantial finite challenge**: peak0.1WP,6.4Goryachev or1Holmes,1000native units on, then13000off, with unchanged native kinetics and conserved-pool densities. Native times/amplitudes are not cross-calibrated. The long exposure is a controlled competition protocol, not an optimized threshold or a guarantee of commitment. Bounded500/2000checks use nine declared anchors from the already reviewed spatial patterns, retaining the exact basal/cue inputs and13000off. Every tested anchor keeps its local endpoint class across the three durations, but global R/angle can change substantially:

{markdown_table(('Model','C','d','Local class, all3durations','R500','R1000','R2000','Angle500°','Angle1000°','Angle2000°'),duration_rows)}

Thus the selected local fates have bounded duration robustness; the global response is not universally saturated. Opposed HolmesC40d20 also changes angle availability with the duration/magnitude threshold. The complete [81-row duration comparison]({links(DURATION/'results/duration_comparison.csv')}) contains all three post-withdrawal readouts, and the [duration verification]({links(DURATION/'results/verification.json')}) records the18newresponses and unchanged native input hashes.

Full-grid local and vector readouts use exactly sampled3000,8000and13000cue-free units. Relative to13000off,18/652local classes differ at3000off (12WP,6Holmes), and4/652still differ at8000off (2WP,2Holmes). The four latter cells are WP C60f.225/C80f.2 and HolmesC35f.35/C70f.275. All these late changes go from two resolved sites to one cue-associated or spanning region. They remain explicit finite-window exceptions;13000off is not a universal settling certificate.

{markdown_table(('Model','Cue-free readout','Maximum absolute ΔR vs13000','Median absolute ΔR','Max-ΔR source cell','Maximum angle change, both unmasked°'),sensitivity)}

The [full local readout changes]({links(DURATION/'results/readout_changed_cells.csv')}), [class-count audit]({links(DURATION/'results/finite_readout_audit.json')}) and [complete vector comparisons](supporting/data/vision_ring_readout_comparisons.csv) preserve those nonmonotone/late outcomes. Angle changes compare only times where both magnitudes pass the same gate. The sampled minima and maximum unmasked angles are not interpolated event times.

## Evidence and reproduction

- [867-row moment summary]({links(MOMENTS/'results/matrix_summary.csv')}), [2608 exact readouts]({links(MOMENTS/'results/readout_sensitivity.csv')}), [lean native-time moment traces]({links(MOMENTS/'results/moment_traces.npz')}) and [read-only extraction audit]({links(MOMENTS/'audit.json')}).
- [New common active-loss basal data]({links(RECOVERY/'results/basal.csv')}), [challenged data]({links(RECOVERY/'results/challenged.csv')}), [all867status rows]({links(RECOVERY/'results/challenged_full_matrix.csv')}) and [new recovery result note]({links(TOPIC/'supporting/metadata/RECOVERY_RESULTS.md')}).
- [Methods and clocks](supporting/metadata/MOMENT_AND_RECOVERY_METHODS.md), [vector exports](supporting/vectors/) and [artifact/source manifest](supporting/metadata/vision_ring_report.manifest.json). The original checkpoint and all old raw data remain preserved; no legacy inactive-withdrawal coefficient fills a new main recovery panel.

The full-grid moments cover434002native frames; background-cancellation relative error is below2.92×10⁻¹⁵ and retained legacy q₁ agrees to3.34×10⁻¹⁶. The prior eight-case summaries agree to floating-point precision. All702recovery states, exact native source snapshots/kicks and7020branches pass independent checks. These audits qualify the declared measurements, not continuum convergence, global stability or a causal decomposition of size and route effects.

Reproduce from the completed audited measurements with the existing environment:

```bash
/home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/vision_ring_report_20261017.py
```

These component sheets and measurement proposals support one updated checkpoint, alongside the all-five linear and replicated-noise results. They do not replace the reserved new WGD motility-results section in the manuscript vision. The [fresh actual-result critique](supporting/metadata/FRESH_RESULT_CRITIQUE.md) records the independent data checks, claim limits and final review disposition. No top-level navigation is changed here.
'''
    # Keep explanatory prose readable while retaining compact scientific IDs
    # and every literal linked filename.
    for prefix in ('all','All','shared','support','peak','sampled','after','before','at','through','below','above','exactly','full','case','states','branches'):
        text=re.sub(r'\b'+prefix+r'(?=[0-9])',prefix+' ',text)
    for suffix in ('on','off','planned','matching','states','branches','conditions','native','rows','frames','cue-free','and'):
        text=re.sub(r'(?<=[0-9])'+suffix+r'\b',' '+suffix,text)
    text=text.replace('twoC40conditions','two C40 conditions').replace('All13','All 13').replace('belowR=','below R = ')
    path=TOPIC/'RESULTS.md';path.write_text(polish(text));ARTIFACTS.append(path)


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--duration-comparison",type=Path,default=DURATION/'results/duration_comparison.csv');parser.add_argument("--no-report",action="store_true");args=parser.parse_args()
    TOPIC.mkdir(parents=True,exist_ok=True);source(__file__)
    matrix=read(SOURCE/"results/matrix.csv");summaries=read(MOMENTS/"results/matrix_summary.csv");readouts=read(MOMENTS/"results/readout_sensitivity.csv")
    source(MOMENTS/"audit.json");source(MOMENTS/"source_sha256.json")
    pairs=read(COMPONENT/"supporting/data/expanded_ring_baseline_pair_measurements.csv")
    basal=primary(read(RECOVERY/"results/basal.csv"));challenged=primary(read(RECOVERY/"results/challenged.csv"))
    recovery_audit=json.loads(source(RECOVERY/"results/independent_data_audit.json").read_text());source(RECOVERY/"results/recovery_audit.json")
    assert len(matrix)==867 and len(basal)==50 and len(challenged)==652 and all(truth(r,"passes") for r in basal+challenged)
    common_pairs=read(RECOVERY/"results/ring_linear_common_recovery_pairs.csv");source(RECOVERY/"results/ring_linear_common_recovery_pairs_audit.json")
    sizes=sorted({number(r,"circumference") for r in matrix})
    plt.rcParams.update({"font.size":9,"axes.titlesize":10,"axes.labelsize":8.5,"axes.spines.top":False,"axes.spines.right":False,"svg.fonttype":"none","pdf.fonttype":42})
    with np.load(source(MOMENTS/"results/moment_traces.npz")) as flat:
        function_sheet(matrix,pairs,basal,challenged,common_pairs,sizes)
        finite_sheet(matrix,summaries,flat,sizes)
        goryachev_support(matrix,summaries,flat)
    moment_maps(matrix,summaries,sizes);basal_support(basal,challenged,sizes);tables=readout_tables(readouts)
    methods()
    if not args.no_report:
        assert args.duration_comparison and args.duration_comparison.exists(),"Provide completed duration comparison"
        durations=read(args.duration_comparison);assert len(durations)==81
        source(DURATION/'results/verification.json');source(DURATION/'results/finite_readout_audit.json');source(DURATION/'results/readout_changed_cells.csv')
        write_results(matrix,pairs,common_pairs,basal,challenged,summaries,readouts,tables,durations)
    links=0
    if not args.no_report:
        for target in re.findall(r'\]\(([^)]+)\)',(TOPIC/'RESULTS.md').read_text()):
            if not target.endswith('vision_ring_report.manifest.json'):
                assert (TOPIC/target.replace('%20',' ')).exists(),target
            links+=1
    manifest=dict(generated_at=datetime.now(timezone.utc).isoformat(),artifact_root="../..",draft_components=args.no_report,
        sha256={str(p.relative_to(TOPIC)):sha(p) for p in ARTIFACTS},source_sha256=SOURCES,
        local_links_checked=links,native_sources_unchanged=True,legacy_inactive_withdrawal_metrics_imported=False,new_integrations=0,
        displayed_recovery_statistic='F_t10(C,d)/F_basal(C)',shared_recovery_color_scale=[0,1.1],
        full_moment_conditions=652,full_planned_matrix=867,main_component_sheets=2,
        primary_vector_gate=.1,vector_gate_sensitivities=[.05,.2])
    path=layout.dev_output_path(TOPIC,"vision_ring_report.manifest.json");path.write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps(dict(artifacts=len(ARTIFACTS),sources=len(SOURCES),report_written=not args.no_report,local_links=links),indent=2))


if __name__=="__main__":main()
