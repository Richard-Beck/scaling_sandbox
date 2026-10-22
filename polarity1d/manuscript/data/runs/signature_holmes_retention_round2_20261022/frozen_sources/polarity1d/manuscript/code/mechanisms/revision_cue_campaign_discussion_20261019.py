"""Two-page author discussion from prior organization and audited pulse brackets."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.text import Text

REPO=Path(__file__).resolve().parents[4]
ROOT=REPO/"polarity1d/manuscript"
F2=ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2"
OUT=F2/"cue_campaign_discussion_20261019"
OLD=F2/"finite_response_candidate_20261019"
CAMPAIGN=ROOT/"data/runs/debelly_corrected_three_test_20261008"
SOURCES={}
NAMES={"wave_pinning":"Wave-pinning","goryachev":"Goryachev"}
COLORS={"wave_pinning":"#0072b2","goryachev":"#d55e00"}
ANCHORS={"wave_pinning":(9.,13.,19.),"goryachev":(9.,13.,19.,21.,25.,39.)}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def source(path):
    path=Path(path);SOURCES[str(path.relative_to(REPO))]=sha(path);return path
def rows(path):
    with source(path).open() as stream:return list(csv.DictReader(stream))
def write_csv(name,data):
    p=OUT/"supporting/data"/name
    with p.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)


def extract():
    data=rows(CAMPAIGN/"results/endpoints.csv")
    manifest=json.loads(source(CAMPAIGN/"manifest.json").read_text())
    cases={c["case_index"]:c for c in manifest["cases"]}
    result=[]
    for model,lengths in ANCHORS.items():
        for length in lengths:
            selected=[]
            for condition in ("fixed_physical","fixed_relative"):
                r,=[r for r in data if r["model"]==model and float(r["length"])==length and r["condition"]==condition and r["definition"]=="full_polarity"]
                other,=[q for q in data if q["case_index"]==r["case_index"] and q["definition"]=="orientation"]
                assert r["status"]==other["status"]=="success"
                for key in ("endpoint","lower_failing_duration","selected_tail_field_pass_0_1_percent"):
                    assert r[key]==other[key]
                assert r["selected_tail_field_pass_0_1_percent"]=="True"
                assert r["nonmonotonic_success_pattern"]=="False" and r["duration_boundary_unresolved"]=="False"
                c=cases[int(r["case_index"])];selected.append((float(c["width"]),r,c))
            (wn,n,cn),(wb,b,cb)=sorted(selected,key=lambda item:item[0])
            assert cn["checkpoint_sha256"]==cb["checkpoint_sha256"] and cn["checkpoint"]==cb["checkpoint"]
            assert cn["cells"]==cb["cells"]==201
            assert cn["response_amplitude"]==cb["response_amplitude"]==(.1 if model=="wave_pinning" else 6.4)
            nf,np_,bf,bp=map(float,(n["lower_failing_duration"],n["endpoint"],b["lower_failing_duration"],b["endpoint"]))
            assert 0<nf<np_ and 0<bf<bp
            result.append(dict(model=model,length_native=length,narrow_support_native=wn,broad_support_native=wb,
                requested_peak_native=cn["response_amplitude"],native_cells=cn["cells"],narrow_last_tested_fail=nf,narrow_first_tested_pass=np_,
                broad_last_tested_fail=bf,broad_first_tested_pass=bp,passing_pulse_ratio=bp/np_,
                bracket_ratio_lower=bf/np_,bracket_ratio_upper=bp/nf,narrow_selected_field_tail_pass=True,broad_selected_field_tail_pass=True,
                orientation_and_strength_brackets_identical=True,narrow_case_index=int(n["case_index"]),broad_case_index=int(b["case_index"]),
                narrow_selected_trace=n["selected_trace"],narrow_trace_prefix=n["selected_trace_prefix"],
                broad_selected_trace=b["selected_trace"],broad_trace_prefix=b["selected_trace_prefix"],
                checkpoint=cn["checkpoint"],checkpoint_sha256=cn["checkpoint_sha256"],
                cue_free_followup_native=10000.,source="data/runs/debelly_corrected_three_test_20261008/results/endpoints.csv"))
    return result


def box(fig,x,y,w,h,text,edge="#8a959d",face="#f7f9fa",size=9,style="-"):
    patch=FancyBboxPatch((x,y),w,h,boxstyle="round,pad=.008",transform=fig.transFigure,
                        linewidth=.9,edgecolor=edge,facecolor=face,linestyle=style)
    fig.add_artist(patch);fig.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=size,linespacing=1.5)
def arrow(fig,a,b):
    fig.add_artist(FancyArrowPatch(a,b,transform=fig.transFigure,arrowstyle="-|>",mutation_scale=11,color="#596871",lw=1.1))
def header(fig,title,subtitle):
    fig.text(.035,.969,title,fontsize=16,weight="bold",va="top")
    fig.text(.035,.928,subtitle,fontsize=9,color="#495760",va="top")


def organization_page():
    fig=plt.figure(figsize=(11.5,8.6))
    header(fig,"Figure 2 discussion: larger incumbents and cue trade-offs",
           "Prior organization result + existing feasibility evidence + proposed campaign · author-review draft · 19 October 2026")
    ax=fig.add_axes((.018,.355,.964,.555));ax.imshow(plt.imread(source(OLD/"figure2_organization.png")));ax.axis("off")
    fig.text(.035,.337,"E–G above are prior work, reused unchanged. The cue experiment below is a proposal for discussion.",fontsize=9.2,color="#495760")
    fig.text(.035,.291,"PROPOSED CAMPAIGN: vary peak, support and pulse duration independently",fontsize=11.2,weight="bold",color="#4d4568")
    box(fig,.035,.135,.17,.105,"Exact native basal\nrestart each trial",size=10)
    arrow(fig,(.211,.188),(.247,.188))
    for left,label in ((.255,"Cue peak"),(.405,"Physical support"),(.555,"Pulse duration")):
        ax=fig.add_axes((left,.145,.12,.095));x=np.linspace(0,1,151)
        if label=="Cue peak":
            p=.5*(1+np.cos(np.pi*x));ax.plot(x,.45*p,color="#adb6bd",lw=1.8);ax.plot(x,p,color="#776590",lw=1.8)
        elif label=="Physical support":
            for width,color in ((.48,"#adb6bd"),(.95,"#776590")):
                p=np.where(x<=width,.5*(1+np.cos(np.pi*np.clip(x/width,0,1))),0.);ax.plot(x,p,color=color,lw=1.8)
        else:
            ax.plot([0,.1,.1,.42,.42,1],[0,0,1,1,0,0],color="#adb6bd",lw=1.8,ls="--")
            ax.plot([0,.1,.1,.85,.85,1],[0,0,1,1,0,0],color="#776590",lw=1.8)
        ax.set_xlim(0,1);ax.set_ylim(-.05,1.1);ax.axis("off")
        fig.text(left+.06,.121,label,ha="center",fontsize=9)
    fig.text(.467,.255,"Three independent controls; cue sketches are schematic",ha="center",fontsize=8.4,color="#5f6870")
    arrow(fig,(.689,.191),(.727,.224));arrow(fig,(.689,.184),(.727,.122))
    box(fig,.741,.194,.223,.076,"During cue: incumbent weakening\nmeasure latency after cue onset",edge="#9a86aa",face="#f9f6fc",style="--",size=8.8)
    box(fig,.741,.088,.223,.076,"After withdrawal: retained replacement\nmeasure at a declared follow-up",edge="#9a86aa",face="#f9f6fc",style="--",size=8.8)
    fig.text(.035,.060,"Regional movies across this grid have not been collected. Pulse duration T and cue-on event latency are different measurements.",fontsize=8.8)
    fig.text(.035,.029,"For author input: optical actuator/native-channel mapping, cue levels, weakening criterion and retention criterion.",fontsize=8.8,color="#4d4568")
    return fig


def bracket(a,b):return f"({a:g}, {b:g}]"


def feasibility_page(data):
    fig=plt.figure(figsize=(11.5,8.6))
    header(fig,"Existing two-width comparisons motivate a targeted cue campaign",
           "Same exact basal and requested peak within each model/length · right end-to-zero raised-cosine patches · no new simulations")
    ax=fig.add_axes((.085,.365,.465,.455))
    ax.axhspan(.53,1,color="#eff6f2");ax.axhspan(1,1.2,color="#fbf2ec");ax.axhline(1,color="#707880",lw=1,ls="--")
    for model in NAMES:
        r=[r for r in data if r["model"]==model]
        x=np.array([r["length_native"] for r in r]);y=np.array([r["passing_pulse_ratio"] for r in r])
        lo=np.array([r["bracket_ratio_lower"] for r in r]);hi=np.array([r["bracket_ratio_upper"] for r in r])
        ax.errorbar(x,y,yerr=np.stack((y-lo,hi-y)),fmt="o" if model=="wave_pinning" else "s",ms=4,markerfacecolor="none",markeredgewidth=1.1,
                    color=COLORS[model],elinewidth=1.35,capsize=3,capthick=1.1,label=NAMES[model],zorder=3)
    ax.set_xlim(7,41);ax.set_ylim(.55,1.18);ax.set_xticks((9,13,19,21,25,39));ax.set_yticks((.6,.7,.8,.9,1.,1.1))
    ax.set_xlabel("Native interval length L",fontsize=10);ax.set_ylabel("Broader / narrower passing-pulse duration",fontsize=10)
    ax.set_title("Selected qualified anchors",loc="left",fontsize=11.5,weight="bold")
    ax.legend(frameon=False,loc="lower right",fontsize=9)
    ax.text(.02,.97,"Above 1: broader requires a longer passing pulse",transform=ax.transAxes,va="top",fontsize=8.1,color="#965a36")
    ax.text(.02,.035,"Below 1: broader permits a shorter passing pulse",transform=ax.transAxes,fontsize=8.1,color="#356c51")
    ax.annotate("same tested bracket",xy=(19,1),xytext=(13.5,1.075),fontsize=8.1,arrowprops=dict(arrowstyle="-",color=".5",lw=.8))
    ax.annotate("L21, pulse 495:\nbroad passes; narrow fails",xy=(21,.99),xytext=(25.4,1.095),fontsize=8.1,
                arrowprops=dict(arrowstyle="-",color=".5",lw=.8))
    ax.grid(alpha=.12)
    fig.text(.085,.855,"Native requested peaks: WP 0.1; Goryachev 6.4. These are not matched optical doses.",fontsize=8.7)
    tableax=fig.add_axes((.594,.352,.383,.465));tableax.axis("off")
    chosen=[("wave_pinning",9), ("wave_pinning",19), ("goryachev",9), ("goryachev",19), ("goryachev",21), ("goryachev",39)]
    cells=[]
    for model,length in chosen:
        r,=[r for r in data if r["model"]==model and r["length_native"]==length]
        cells.append([f"{'WP' if model=='wave_pinning' else 'Gory'}; L{length}",f'{r["narrow_support_native"]:g} → {r["broad_support_native"]:g}',
                      bracket(r["narrow_last_tested_fail"],r["narrow_first_tested_pass"]),bracket(r["broad_last_tested_fail"],r["broad_first_tested_pass"])])
    table=tableax.table(cellText=cells,colLabels=["Model; size","Support","Narrow (F, P]","Broad (F, P]"],
                       colWidths=(.20,.17,.315,.315),cellLoc="center",loc="center",bbox=(0,.03,1,.94))
    table.auto_set_font_size(False);table.set_fontsize(7.8)
    for (i,j),cell in table.get_celld().items():
        cell.set_edgecolor("#d8dfe3");cell.set_linewidth(.6)
        cell.set_facecolor("#edf1f4" if i==0 else ("#f9fafb" if i%2 else "white"))
        if i==0:cell.get_text().set_weight("bold")
    fig.text(.594,.838,"Actual tested native pulse brackets",fontsize=11.1,weight="bold")
    fig.text(.594,.320,"Every shown endpoint: native field-tail PASS.",fontsize=9.1,weight="bold",color="#356c51")
    fig.text(.594,.269,"Both scalar gates give the same displayed brackets.\nF = last tested fail; P = first tested pass.",fontsize=8.7,linespacing=1.45)
    fig.text(.085,.262,"Point = passing-pulse ratio; whiskers use both tested brackets.\nBracket envelope, not CI. Exact construction is in the caption.",fontsize=8.4,linespacing=1.5)
    box(fig,.066,.104,.43,.111,"Existing feasibility\nBroader helps WP at the selected sizes; it delays smaller Gory,\nshares the L19 bracket, then helps by L21 and at larger sizes.",size=9)
    box(fig,.532,.104,.43,.111,"Measurement gap\nThreshold searches retain scalar histories, tail summaries and\nfinal fields, rather than full regional movies across the grid.",size=9)
    fig.text(.066,.064,"Increasing support at fixed requested peak also increases integrated cue. These comparisons sample two widths per size, rather than a complete width scan.",fontsize=8.4)
    fig.text(.066,.033,"Existing success: persistent reversed scalar polarity near its initial magnitude after 10,000 cue-free units.",fontsize=8.4)
    fig.text(.066,.012,"Field-tail qualification is separate; whole-field replacement is not implied. N = 201; no spatial-convergence claim.",fontsize=8,color="#596871")
    return fig


def check_page_bounds(fig,page):
    fig.canvas.draw();renderer=fig.canvas.get_renderer();bounds=fig.bbox
    checks=[]
    for text in fig.findobj(match=Text):
        if not text.get_visible() or not text.get_text():continue
        if text.axes is not None and not text.axes.axison:continue
        extent=text.get_window_extent(renderer)
        assert extent.x0>=bounds.x0 and extent.y0>=bounds.y0 and extent.x1<=bounds.x1 and extent.y1<=bounds.y1,(page,text.get_text(),extent.bounds,bounds.bounds)
        checks.append((extent.x0,extent.y0,bounds.x1-extent.x1,bounds.y1-extent.y1))
    return dict(page=page,native_text_boxes_checked=len(checks),all_within_page=True,
                minimum_margin_points=float(np.min(checks))*72/fig.dpi,
                prior_organization_raster_contained=page==1)


def main():
    for group in ("data","vectors","metadata","viewers"):(OUT/"supporting"/group).mkdir(parents=True,exist_ok=True)
    data=extract();write_csv("selected_width_brackets.csv",data)
    source(OLD/"supporting/data/organization_all97.csv");source(ROOT/"code/core/polarity_establishment.py")
    source(ROOT/"code/core/corrected_reversal.py");source(__file__)
    before={str(p.relative_to(REPO)):sha(p) for p in OLD.rglob("*") if p.is_file()}
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
    figures=[organization_page(),feasibility_page(data)]
    text_bounds=[check_page_bounds(fig,i) for i,fig in enumerate(figures,1)]
    pdf=OUT/"supporting/vectors/figure2_cue_campaign_discussion.pdf"
    with PdfPages(pdf) as book:
        for i,fig in enumerate(figures,1):
            fig.savefig(OUT/f"figure2_cue_campaign_page{i}.png",dpi=220)
            fig.savefig(OUT/"supporting/vectors"/f"figure2_cue_campaign_page{i}.svg")
            book.savefig(fig);plt.close(fig)
    print(str(pdf),flush=True)
    assert all(sha(REPO/p)==h for p,h in before.items())
    generated=datetime.now(timezone.utc).isoformat()
    (OUT/"README.md").write_text("# Figure 2 cue-campaign discussion\n\nA two-page author-review draft, not a finished main-figure claim. No new native integration.\n\n[Two-page PDF](supporting/vectors/figure2_cue_campaign_discussion.pdf)\n\n![Prior organization and proposed campaign](figure2_cue_campaign_page1.png)\n\n![Checked existing cue-width feasibility](figure2_cue_campaign_page2.png)\n\nPage 1 retains the prior all-five E–G organization artwork unchanged and proposes independent cue peak/support/duration controls. Dashed outcome boxes are proposed regional measurements, not completed grid coverage. Page 2 uses a few audited two-width anchors and exact tested fail/pass brackets; the ratio whiskers are bracket envelopes, not confidence intervals. The proposed optical actuator and regional weakening/retention criteria remain for author input.\n\n[Exact selected brackets and source identities](supporting/data/selected_width_brackets.csv) · [Caption](CAPTION.md) · [Source and artifact hashes](supporting/metadata/validation.json) · [Original native-unit organization art](../../../../../../../../current%20manuscript%20figs/initial_polarity_length.png)\n")
    (OUT/"CAPTION.md").write_text("Figure 2 discussion checkpoint. Page 1 reuses prior all-five basal organization E–G without changing its 97 measurements or own-model references. It motivates a proposed campaign in which cue peak, physical support and pulse duration vary independently from exact native basals. The optical actuator, levels and regional outcome criteria remain open. Cue-on weakening latency is distinct from pulse duration; full incumbent-region movies have not been retained throughout the existing threshold searches. Page 2 compares broader and narrower right raised-cosine patches from identical qualified basals, at fixed requested native peak (WP 0.1; Goryachev 6.4). Points are ratios of tested passing pulses. Each whisker is the envelope obtained from the two last-fail/first-pass brackets, not a statistical confidence interval. All shown anchors pass the separate native field-tail gate and have identical orientation/strength brackets. At Gory L21, pulse495 passes for the broader patch and fails for the narrower patch; touching bracket bounds do not erase that observed contrast. Fixed requested peak is not fixed integrated input or a cross-model optical calibration. Native sampled peaks differ slightly with support. Existing scalar success after10000 cue-free units does not establish whole-field replacement. Only two widths per size are sampled; the independent three-control campaign and regional measurements remain proposed.\n")
    with (OUT/"CAPTION.md").open("a") as stream:
        stream.write(r"""

For narrower and broader cues, subscripts n and b, the displayed point and bracket envelope are

\[
r_P=\frac{P_{\mathrm b}}{P_{\mathrm n}},\qquad
\mathcal B=\left[\frac{F_{\mathrm b}}{P_{\mathrm n}},\frac{P_{\mathrm b}}{F_{\mathrm n}}\right].
\tag{1}
\]

The original native-unit organization artwork and 97 basal measurements remain linked from the README. The selected width comparisons use the existing N = 201 native grids; the intervals are tested temporal brackets, without an independent spatial-convergence claim. The chosen WP anchors stop at L19: larger selected scalar endpoints have failed field-tail flags or fixed-width censoring, so they are not used as headline broad-cue rescues here.
""")
    artifacts={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob("*") if p.is_file() and p.name!="validation.json"}
    record=dict(generated_at=generated,discussion_draft=True,new_integrations=0,selected_pairs=len(data),all_selected_field_tail_flags_pass=True,
                both_scalar_gate_brackets_match=True,all_pairs_exact_same_checkpoint=True,old_candidate_files_preserved=len(before),
                text_bounds=text_bounds,
                source_sha256=SOURCES,artifact_sha256=artifacts)
    (OUT/"supporting/metadata/validation.json").write_text(json.dumps(record,indent=2)+"\n")
    print(json.dumps({k:v for k,v in record.items() if k not in ("source_sha256","artifact_sha256")},indent=2))


if __name__=="__main__":main()
