"""Collect the completed native atlas and publish the audited Figure 3 reader.

This is read-only analysis of frozen native evidence. It neither integrates
models nor substitutes missing workers for finite-horizon lower bounds.
"""
from pathlib import Path
from datetime import datetime, timezone
from urllib.parse import quote
import argparse
import html
import json
import os
import re
import shutil

import numpy as np
import pandas as pd

import signature_figure3_atlas_20261020 as atlas
from signature_figure3_heading_20261020 import direction_intervals

OUT, RUNS, sha, read, write = atlas.OUT, atlas.RUNS, atlas.sha, atlas.read, atlas.write
SYNTHESIS = RUNS/"signature_figure3_synthesis_20261020"
NAMES = dict(wave_pinning="Wave-pinning",goryachev="Goryachev",holmes_model3="Holmes 3",spring="Spring")


def link(path, base):
    return quote(os.path.relpath(Path(path), Path(base)), safe="/._-")


def readable_text(document):
    """Separate prose words/numbers without changing paths, code or MathML."""
    parts=re.split(r'(<[^>]+>)',document)
    protected=False
    for i,part in enumerate(parts):
        if part.startswith('<'):
            if re.match(r'<(?:code|math|script|style)(?:\s|>)',part):protected=True
            elif re.match(r'</(?:code|math|script|style)>',part):protected=False
        elif not protected:
            part=re.sub(r'(?<=[a-z])(?=\d)',' ',part)
            parts[i]=re.sub(r'(?<=\d),(?=\d)',', ',part)
    return ''.join(parts)


def interval_rates(time,vector):
    """Principal saved-interval angle changes; never unwrap across a gap."""
    t,z=np.asarray(time),np.asarray(vector)
    if len(t)!=len(z) or len(t)<2 or np.any(np.diff(t)<=0):
        raise ValueError("Aligned native moments and increasing times required")
    increment=np.angle(z[1:]*np.conj(z[:-1]))
    rates={}
    for floor in (1e-10,.05,.1,.2):
        frames,valid,chord_minimum=direction_intervals(z,floor)
        value=increment/np.diff(t)
        value[~valid]=np.nan
        rates[str(floor)]=value
    return .5*(t[:-1]+t[1:]),rates


def rate_export(case,status,source,out):
    path=source/status["readouts"]
    assert sha(path)==status["readouts_sha256"]
    with np.load(path) as z:
        time=z["time"].copy();vector=z["vector_real"]+1j*z["vector_imag"]
    middle,rates=interval_rates(time,vector)
    base=out/"supporting/viewers/vector_rates";base.mkdir(exist_ok=True)
    arrays=base/(case["case_id"]+".npz")
    np.savez_compressed(arrays,time=time,vector_real=vector.real,vector_imag=vector.imag,interval_midpoint=middle,
        **{"angular_rate_floor_"+key:value for key,value in rates.items()})
    values=np.concatenate([r[np.isfinite(r)] for r in rates.values()])
    scale=max(float(np.percentile(abs(values),75))/5 if len(values) else 0.,1e-8)
    ordinate=np.arcsinh(values/scale)
    lo,hi=(float(ordinate.min()),float(ordinate.max())) if len(values) else (-1.,1.)
    lo,hi=min(lo,0.),max(hi,0.)
    if hi-lo<1e-4:lo,hi=-1.,1.
    margin=.08*(hi-lo);lo-=margin;hi+=margin
    temporal_scale=1. if case["model"]=="spring" else 10.
    xx=90+880*np.log1p(middle/temporal_scale)/np.log1p(time[-1]/temporal_scale)
    transform=lambda value: 340-240*(np.arcsinh(value/scale)-lo)/(hi-lo)
    paths=[]
    for floor,color,width in [('1e-10','#a1a8ab',2.5),('0.05','#21689a',1.8),('0.1','#d17b2b',1.2),('0.2','#80558f',.7)]:
        value=rates[floor];commands=[];active=False
        for x,y in zip(xx,value):
            if np.isfinite(y):commands.append(f'{"L" if active else "M"}{x:.2f},{transform(y):.2f}');active=True
            else:active=False
        paths.append(f'<path d="{" ".join(commands)}" fill="none" stroke="{color}" stroke-width="{width}"/>')
    ticks=[]
    for point in np.linspace(lo,hi,7):
        y=340-240*(point-lo)/(hi-lo);value=scale*np.sinh(point)
        ticks.append(f'<path d="M85,{y:.2f}H970" stroke="#e0e4e8"/><text x="80" y="{y+5:.2f}" text-anchor="end">{value:.3g}</text>')
    for point in np.linspace(0,np.log1p(time[-1]/temporal_scale),6):
        x=90+880*point/np.log1p(time[-1]/temporal_scale);t=temporal_scale*np.expm1(point)
        ticks.append(f'<text x="{x:.2f}" y="370" text-anchor="middle">{t:.5g}</text>')
    title=html.escape(case['case_id']+' · '+NAMES[case['model']]+' · saved-interval global heading rate')
    document=f'''<svg xmlns="http://www.w3.org/2000/svg" width="1050" height="460" viewBox="0 0 1050 460"><rect width="1050" height="460" fill="white"/><g font-family="sans-serif" font-size="14" fill="#15212a"><text x="525" y="28" text-anchor="middle" font-size="20">{title}</text><text x="525" y="54" text-anchor="middle">Grey: exact direction; blue: R ≥ 0.05; orange: R ≥ 0.1; purple: R ≥ 0.2</text>{''.join(ticks)}<path d="M90,100V340H970" fill="none" stroke="#15212a"/>{''.join(paths)}<text x="525" y="405" text-anchor="middle">Native time, log(1 + t / {temporal_scale:g}) axis; all retained intervals</text><text transform="translate(20,220) rotate(-90)" text-anchor="middle">Principal Δheading / Δtime (rad / native time)</text><text x="525" y="435" text-anchor="middle" font-size="12">Signed asinh rate axis; no instantaneous RHS rate, no unobserved full turns, no unwrapping across cancellation/quality gaps.</text></g></svg>'''
    target=base/(case["case_id"]+".svg");target.write_text(document)
    return [dict(label="Global heading-rate time course",path=str(target),sha256=sha(target)),
        dict(label="Full saved-interval heading-rate arrays and quality variants",path=str(arrays),sha256=sha(arrays))]


def collect(run, out):
    plan = read(run/"manifest.json")
    assert sha(run/"frozen_renderer.py") == plan["renderer_sha256"]
    portable_audit=read(out/"supporting/metadata/lossless_pixel_audit.json")
    assert portable_audit["complete"]
    portable={row["case_id"]:row for row in portable_audit["records"]}
    cases, completed = plan["cases"], []
    for slot, indices in enumerate(plan["chunks"]):
        record = read(run/"workers"/f"{slot:04d}.json")
        assert record["complete"] and len(record["records"]) == len(indices)
        for index, result in zip(indices, record["records"]):
            case = cases[index]
            assert result["case_id"] == case["case_id"]
            assert sha(Path(result["png"])) == result["png_sha256"]
            webp=portable[case["case_id"]]
            assert webp["pixel_equal"] and webp["png_sha256"]==result["png_sha256"]
            assert sha(Path(webp["webp"]))==webp["webp_sha256"]
            assert sha(Path(case["source_status"])) == case["source_status_sha256"]
            status, source = read(Path(case["source_status"])), Path(case["source_run"])
            assert status["status"] == "complete"
            files = [dict(label="Native status, gates and event records", path=case["source_status"],sha256=case["source_status_sha256"]),
                dict(label="Frozen native campaign plan",path=str(source/"manifest.json"),sha256=case["source_manifest_sha256"])]
            for phase in status["phases"]:
                files.append(dict(label=f'Full native species: {phase["kind"]}',path=str(source/phase["path"]),sha256=phase["sha256"]))
            if "readouts" in status:
                files.append(dict(label="Saved scalar/vector diagnostics",path=str(source/status["readouts"]),sha256=status["readouts_sha256"]))
            if "region_trace" in case:
                files.append(dict(label="Corrected complete-region diagnostics",path=case["region_trace"],sha256=case["region_trace_sha256"]))
            if status.get("exact_rhs_archive"):
                files.append(dict(label="Exact cue-free native RHS at withdrawal",path=str(source/status["exact_rhs_archive"]),sha256=status["exact_rhs_sha256"]))
            files.extend(rate_export(case,status,source,out))
            for f in files:
                assert Path(f["path"]).is_file(), f["path"]
            completed.append(dict(**case, **{k:v for k,v in result.items() if k not in case},
                webp=webp["webp"],webp_sha256=webp["webp_sha256"],lossless_rgba_pixel_equal=True,files=files))
    assert len(completed) == len(cases) == len({r["case_id"] for r in completed})
    counts = {stage:sum(r["stage"] == stage for r in completed) for stage in sorted({r["stage"] for r in completed})}
    summary = dict(generated_utc=datetime.now(timezone.utc).isoformat(), producer_state="complete; all planned workers verified",
        cases=len(cases),workers=len(plan["chunks"]),stage_counts=counts,
        atlas_manifest_sha256=sha(run/"manifest.json"),renderer_sha256=plan["renderer_sha256"],
        lossless_webp_pixel_equal_images=portable_audit["pixel_equal_images"],lossless_webp_audit_sha256=sha(out/"supporting/metadata/lossless_pixel_audit.json"),
        collector_source_sha256=sha(Path(__file__)),png_and_status_hashes_verified=len(completed),
        full_archive_hash_scope="Phase hashes were checked by each atlas producer and the independent native audits; collection verifies terminal status and PNG hashes and all linked native paths.")
    write(run/"results/complete_atlas_audit.json",summary)
    write(out/"supporting/metadata/native_atlas_audit.json",summary)
    write(out/"supporting/data/native_history_index.json",dict(summary=summary,cases=completed))
    return completed, summary


def controls(out):
    main=pd.read_csv(SYNTHESIS/"B_complete_region_evidence_and_geometry.csv")
    control=pd.read_csv(SYNTHESIS/"B_control_evidence_and_geometry.csv")
    rows=[]
    for c in control.to_dict("records"):
        selected=main[(main.model==c["model"]) & (main.size_parameter==c["size_parameter"]) &
            (abs(main.initial_gap-c["initial_gap"])<1e-8) & (main.target_ratio==c["target_ratio"]) & main.control.isna()]
        assert len(selected)==1
        r=selected.iloc[0].to_dict()
        row={k:c[k] for k in ["model","size_parameter","control","initial_gap","target_ratio","slot"]}
        for key in ["seed_status","winner","primary_claim_eligible","resolution_lower","resolution_upper","resolution_status","withdrawal_time","measured_gap","finite_tail_pass"]:
            row["reference_"+key],row["control_"+key]=r[key],c[key]
        row.update(seed_status_equal=r["seed_status"]==c["seed_status"],winner_equal=r["winner"]==c["winner"],
            eligibility_equal=r["primary_claim_eligible"]==c["primary_claim_eligible"],
            relative_withdrawal_time_error=abs(c["withdrawal_time"]/r["withdrawal_time"]-1) if r["withdrawal_time"]>0 else np.nan,
            resolution_brackets_overlap=bool(max(r["resolution_lower"],c["resolution_lower"]) <= min(r["resolution_upper"],c["resolution_upper"]))
                if np.isfinite(r["resolution_upper"]) and np.isfinite(c["resolution_upper"]) else None)
        rows.append(row)
    frame=pd.DataFrame(rows)
    frame.to_csv(out/"supporting/data/B_control_convergence.csv",index=False)
    return frame


CSS="""body{margin:0;background:#f4f6f8;color:#15212a;font:17px/1.55 system-ui,sans-serif}main{max-width:1160px;margin:auto;padding:30px}h1{font-size:30px;line-height:1.2}h2{margin-top:42px;font-size:23px}p{max-width:1050px}.status,.note{background:#e5eef3;border-left:4px solid #26768e;padding:12px 18px}a{color:#12637a}img.plate{width:100%;background:white;border:1px solid #cbd5dc}nav{display:flex;gap:20px;flex-wrap:wrap}table{border-collapse:collapse;width:100%;font-size:15px}th,td{padding:8px;border:1px solid #ccd4d9;text-align:left}th{background:#e9eef1}.equation{display:flex;justify-content:center;align-items:center;margin:25px 0;font-family:serif;font-size:23px;position:relative}.eqno{position:absolute;right:0;font:16px system-ui}details{border:1px solid #ccd4d9;margin:12px 0;padding:12px}select,input{padding:8px;margin:6px;font-size:16px}code{font-size:13px;overflow-wrap:anywhere}figcaption{font-size:15px}.small{font-size:14px}button{padding:8px 14px} .viewer{max-width:1320px}ul{padding-left:24px}"""


def atlas_viewer(cases,summary,out):
    base=out/"supporting/viewers"
    entries=[]
    for case in cases:
        title=f'{case["case_id"]} · {NAMES[case["model"]]} · size {case["size_parameter"]:g} · initial gap {case["initial_gap"]:.4g}'
        if case.get("target_ratio") is not None:title+=f' · ratio {case["target_ratio"]:g}'
        if case.get("multiplier") is not None:title+=f' · pulse {case["multiplier"]:g}×'
        if case.get("control"):title+=' · '+case["control"]
        entries.append(dict(id=case["case_id"],title=title,stage=case["stage"],model=case["model"],
            classification=case["classification"],image=link(case["webp"],base),png_sha256=case["png_sha256"],
            full_horizon=case["full_horizon"],native_frames=case["native_frames"],
            rate_image=link(next(f["path"] for f in case["files"] if f["label"]=="Global heading-rate time course"),base),
            files=[dict(label=f["label"]+(" (local/HPC data)" if Path(f["path"]).suffix==".npz" or not Path(f["path"]).is_relative_to(out) else ""),url=link(f["path"],base),sha256=f["sha256"]) for f in case["files"]]))
    document='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Figure 3 · Complete native histories</title><style>'''+CSS+'''</style><main class="viewer">
<nav><a href="../../index.html">Figure 3 scientific reader</a><a href="../data/native_history_index.json">Complete source index</a></nav>
<h1>Every retained native history</h1><p class="status">All 1,150 planned histories are complete. Producer stopped; this page has no timed refresh. Generated '''+html.escape(summary["generated_utc"])+'''.</p>
<p class="note">All displayed histories use portable lossless WebP images, individually audited for exact RGBA equality against the original PNGs. The global-rate SVGs are portable. Links labeled local/HPC data point to retained native species, exact RHS, source statuses and full rate arrays; those data archives are kept separately from the compact figure checkout. <a href="../metadata/REBUILD.txt">Rebuild instructions and data requirements</a>.</p>
<p>The native field uses the active species, including active Rac for Holmes 3. Basins cover activity beyond the cue. Spatial labels describe association with initial incumbent and cue positions, not molecular ancestry. A one-region relocation is a separate topology outcome. The full saved native species archives and exact withdrawal RHS are linked below each history.</p>
<label>Arm <select id="stage"><option value="">All</option>'''+''.join(f'<option>{s}</option>' for s in summary["stage_counts"])+'''</select></label>
<label>Model <select id="model"><option value="">All</option>'''+''.join(f'<option value="{k}">{v}</option>' for k,v in NAMES.items())+'''</select></label>
<label>Find <input id="search" placeholder="gap, size, ratio or case ID"></label><br>
<label>History <select id="case" style="max-width:100%"></select></label><button id="back">Previous</button><button id="next">Next</button>
<p id="count" class="small"></p><h2 id="title"></h2><img id="image" class="plate" alt="Native field, profiles, complete region amounts and vector diagnostics"><img id="rate_image" class="plate" alt="Global heading-rate time course with all four quality gates"><p id="info" class="small"></p><div id="sources"></div>
<script>const cases='''+json.dumps(entries,separators=(',',':')).replace('</',r'<\/')+''';let selected=[];
const byId=id=>document.getElementById(id);function filter(){const s=byId('stage').value,m=byId('model').value,q=byId('search').value.toLowerCase();selected=cases.filter(c=>(!s||c.stage===s)&&(!m||c.model===m)&&c.title.toLowerCase().includes(q));byId('case').replaceChildren(...selected.map(c=>{const o=document.createElement('option');o.value=c.id;o.textContent=c.title;return o}));byId('count').textContent=selected.length+' matching histories / '+cases.length+' verified total';show()}
function show(){const c=selected.find(c=>c.id===byId('case').value);if(!c)return;byId('title').textContent=c.title;byId('image').src=c.image;byId('rate_image').src=c.rate_image;byId('info').textContent=c.native_frames+' retained full-field samples; native horizon '+c.full_horizon+'; PNG SHA256 '+c.png_sha256;const list=document.createElement('ul');for(const f of c.files){const li=document.createElement('li'),a=document.createElement('a'),small=document.createElement('code');a.href=f.url;a.textContent=f.label;small.textContent=' · SHA256 '+f.sha256;li.append(a,small);list.append(li)}byId('sources').replaceChildren(list);history.replaceState(null,'','#'+c.id)}
for(const id of ['stage','model','search'])byId(id).addEventListener(id==='search'?'input':'change',filter);byId('case').addEventListener('change',show);function step(d){const s=byId('case');s.selectedIndex=Math.max(0,Math.min(s.options.length-1,s.selectedIndex+d));show()}byId('back').onclick=()=>step(-1);byId('next').onclick=()=>step(1);const initial=location.hash.slice(1);filter();if(initial&&cases.some(c=>c.id===initial)){byId('case').value=initial;show()}
</script></main></html>'''
    (base/"native_history_viewer.html").write_text(document)
    (out/"supporting/metadata/REBUILD.txt").write_text("""Figure3 native-history and figure reconstruction

The portable reader includes all1150 lossless WebP histories, global-rate SVGs,
two main plates/PDF, seven supplementary pages and complete analysis CSVs.
Original PNGs, full native-species NPZ archives, exact withdrawal RHS arrays,
per-case native statuses and full saved-interval rate NPZs are retained locally
under the manuscript data/runs and topic directories. They are not implied to
be included in a compact Git clone. Local-data links are explicitly labeled.

Reconstruction needs the retained signature_figure3_* runroots and immutable
producer manifests, native basals and phase arrays named in native_history_index.json.
Missing raw data must be restored from that scientific archive; do not invent
or replace it with a finite-cap value. Existing live outputs are reproducible,
while native producers require the existing scaling_sandbox Conda environment
and SLURM. No package installation or retuned model parameters are required.

From the repository root, use code under polarity1d/manuscript/code/mechanisms/:
  signature_figure3_synthesis_20261020.py --include-A2 --include-B
  signature_figure3_heading_20261020.py
  signature_figure3_render_20261020.py
The complete native atlas uses signature_figure3_atlas_20261020.py freeze,
then the emitted SLURM payload. Use a fresh output run directory when the
default immutable plan exists. Its frozen renderer checks every source phase.
Then signature_figure3_lossless_atlas_20261020.py freeze/worker/collect exports
pixel-identical lossless WebPs. Run signature_figure3_reader_20261020.py after
both producers have completed and their independent audits are present.
Finally organize_length_dev.py --group revision_program_20261020/figure3_distance
updates only this topic. Never run an unscoped organizer as a reconstruction step.

The native PNG/WebP pixel audit and native phase hashes are separate evidence.
Lossless image equality does not replace the independent native-field checks.
""")


def reader(out,summary,comparison):
    base=out
    A1=pd.read_csv(SYNTHESIS/"A1_complete_region_evidence.csv")
    B=pd.read_csv(SYNTHESIS/"B_complete_region_evidence_and_geometry.csv")
    heading=pd.read_csv(SYNTHESIS/"A1_direct_heading_check.csv")
    data_links=''.join(f'<li><a href="supporting/data/{p.name}">{html.escape(p.stem.replace("_"," "))}</a></li>' for p in sorted((out/"supporting/data").glob("*.csv")))
    auditdefs=[("Original A1 full-field audit","signature_figure3_campaign_20261020/results/primary_independent_native_audit.json"),
        ("Fresh basal qualifications and long no-cue controls","signature_figure3_basals_20261020/results/basal_audit.json"),
        ("Corrected continuous handover refinements","signature_figure3_handover_refinement_v2_20261020/results/independent_audit.json"),
        ("A2 independent full-field audit","signature_figure3_regional_A2_20261020/results/all_independent_native_audit.json"),
        ("B independent full-field audit","signature_figure3_regional_B_20261020/results/all_independent_native_audit.json"),
        ("B mirrored-placement comparison","signature_figure3_regional_B_20261020/results/mirror_comparison.json"),
        ("B numerical controls independent audit","signature_figure3_regional_B_controls_20261020/results/all_independent_native_audit.json")]
    audit_links=''
    for label,relative in auditdefs:
        path=RUNS/relative
        if not path.is_file():
            # Exact refinement audit filename is kept in the producer record.
            candidates=sorted((path.parent).glob("*audit*.json"))
            assert len(candidates)==1,(str(path),[str(c) for c in candidates])
            path=candidates[0]
        audit_links+=f'<li><a href="{link(path,base)}">{label}</a> <code>{sha(path)}</code></li>'
    medians=[]
    for model,name in NAMES.items():
        p=heading[(heading.model==model)&heading.control.isna()&(heading.classification=="external")&(heading.coherence_at_handover>=.05)]
        medians.append(f'<tr><td>{name}</td><td>{p.handover_angle_fraction_of_requested.median():.3f}</td><td>{int(p.handover_angle_fraction_of_requested.notna().sum())}</td></tr>')
    main=B[B.control.isna()]
    bcount=main.groupby("model").agg(eligible=("primary_claim_eligible","sum"),total=("slot","count"))
    brows=''.join(f'<tr><td>{NAMES[m]}</td><td>{int(r.eligible)}</td><td>{int(r.total-r.eligible)}</td></tr>' for m,r in bcount.iterrows())
    unseeded=main[main.seed_lower_bound.notna()]
    seed_bounds=''.join(f'<tr><td>{NAMES[r.model]}</td><td>{r.size_parameter:g}</td><td>{r.initial_gap:.6g}</td><td>{r.target_ratio:g}</td><td>≥ {r.seed_lower_bound:g}</td></tr>' for r in unseeded.itertuples())
    text='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Figure 3 · Separation and native polarity competition</title><style>'''+CSS+'''</style><main>
<nav><a href="supporting/vectors/figure3_distance_main.pdf">Two-page main figure PDF</a><a href="supporting/vectors/figure3_distance_supplement.pdf">Seven-page supplement</a><a href="supporting/viewers/native_history_viewer.html">All 1,150 native histories</a><a href="#methods">Measurements and methods</a><a href="#sources">Source tables and audits</a></nav>
<h1>Separation affects native competition differently across models</h1>
<p class="status">Complete audited Figure 3 campaign. '''+html.escape(summary["generated_utc"])+'''. All native producers and the native-history renderer have stopped. This reader does not refresh automatically.</p>
<p>Rotating the same cue at a fixed circumference separates distance from enlargement. Under continuous forcing, Goryachev handover slows with separation, whereas Holmes 3 is nearly flat over its valid two-region placements. Wave-pinning and Spring show a near-edge increase followed by a broad plateau. Cue-free competition answers a different question: Wave-pinning, Holmes 3 and Spring can coarsen very slowly, with all completed-horizon lower bounds retained. These observations support model-dependent responses, rather than a universal absence of a distance effect in Wave-pinning.</p>
<p>The geometric translation requires measured native activity directions. Complete-region handover typically occurs around half the requested cue angle. A center-distance remap therefore does not itself measure time to turn through that angle, and enlargement does not guarantee greater directional persistence from these assays alone.</p>
<figure><a href="figure3_distance_main.png"><img class="plate" src="figure3_distance_main.png" alt="Fixed-size ring cue placement, native histories and distance responses"></a><figcaption>Continuous-cue handover and cue-free resolution use their own measured gaps. All four native models and all three sizes appear. Open symbols at negative gap are borderline placements; they do not enter the external-placement distance claim. Crosses show unavailable two-region handovers. Triangles preserve completed finite-cap lower bounds. Lines connect only adjacent, eligible uncensored placements; missing, merged and single-region outcomes are not bridged.</figcaption></figure>
<figure><a href="figure3_distance_translation.png"><img class="plate" src="figure3_distance_translation.png" alt="Measured-width geometry remaps and direct angle checks"></a><figcaption>The fixed-angle panel shows handover time, with geometry-only evaluation of the small-cell curve at larger-cell measured-width gaps. Dashed curves stop within observed support. The inset checks the actual global heading at handover. The fixed-physical-scale panel divides a conditional geometric center step by measured resolution time; it is not a measured cell angular velocity. Native endpoint headings and full-angle target times are retained in the supplement and tables.</figcaption></figure>
<h2>Standalone figure legend</h2>
<p><strong>A–C.</strong> Each model starts from a fresh, qualified single-region basal, with unchanged native kinetics. Chemical circumferences are30,37.8 and47.6; Spring resting-size labels are2.5,3.15 and3.97, with its actual circumference recorded throughout. One raised-cosine native cue is rotated around each basal. The primary initial edge gap is the shortest arc between the incumbent's initial10%-contrast edge and the cue support's near edge. Negative gaps straddle that edge. A1 begins with5,000 native chemical units and extends to20,000 where the retained continuation gate requires it; Spring uses2,000 units. Every timing bound uses that case's actual completed horizon. Complete-region handover is the first continuous crossing of challenger excess above incumbent excess while two resolved associated regions exist. A2 tests0.5,0.75,1,1.5 and2 times that placement's own handover, then13,000 cue-free chemical units or1,300 Spring units. No cap substitutes for an unavailable handover. B seeds to complete-region excess ratios0.5,1 and2, withdraws the cue, and measures the first single-resolved-peak time over a20,000-unit chemical or2,000-unit Spring tail. Its horizontal axis uses the measured current global10%-contrast edge gap at withdrawal. Two peaks with merged primary contours remain a separate ineligible outcome. Time brackets reflect retained sample intervals; unresolved runs are lower bounds at their completed caps.</p>
<p><strong>D–F.</strong> Native incumbent width and cue support enter the conversion from initial gap to requested center angle. Geometry-only remaps evaluate the reference-size handover curve at the larger cell's gap, without extrapolation or interpolation through unavailable/censored neighbors. They are distinct from measured size shifts at matched gap. Global direction comes from the native active-field first spatial harmonic. Coherence is the instantaneous first-harmonic magnitude divided by mean active field; a0.05 reporter gate is displayed, with exact-direction and0.1/0.2 alternatives tabulated. It differs from Figure4's basal-moment-normalized polarity magnitude. Full target-angle timing is unavailable across cancellation/quality gaps;50%,90% and99% targets provide sensitivity. Chemical native time is not converted to physical seconds. These polarity-activity readouts remain distinct from observed migration and conditional OU persistence estimates.</p>
<h2>Competition and angular reorientation are different measurements</h2>
<table><thead><tr><th>Model</th><th>Median actual handover angle / requested angle</th><th>External cases with usable handover heading</th></tr></thead><tbody>'''+''.join(medians)+'''</tbody></table>
<p>The near-half-angle result follows from the directly measured first moment, not from an assumed cell-velocity law. A concrete cue-free Wave-pinning example at C30, initial gap0.5, ratio1 has a potential center step of about1.54 radians but only about0.75 radians of basal-to-final vector rotation. The final peak can occupy an intermediate direction. The global-vector checks therefore qualify, rather than automatically confirm, the proposed geometry-to-persistence translation.</p>
<h2 id="methods">Measurements and native protocols</h2>
<p>The amount is native active concentration integrated over a complete spatial region, after subtracting the common instantaneous minimum. Holmes3 uses active Rac; the other models use their native active species. In a two-peak field, the periodic basins meet at the two intervening saddles. Boundary position comes from the local three-point parabolic minimum; the actual field is integrated as a piecewise-linear profile, including fractional cells. The boundary is not asserted to be the exact minimum of that piecewise-linear field. The two basins partition the complete native excess. Active flat or nonunique saddles are unscorable. A single moving region, a merged region, and cue-phase takeover receive explicit alternate outcomes.</p>
<div class="equation"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mrow><msub><mi>E</mi><mi>j</mi></msub><mo>(</mo><mi>t</mi><mo>)</mo><mo>=</mo><msub><mo>∫</mo><mrow><msub><mi>Ω</mi><mi>j</mi></msub><mo>(</mo><mi>t</mi><mo>)</mo></mrow></msub><mo>[</mo><mi>a</mi><mo>(</mo><mi>x</mi><mo>,</mo><mi>t</mi><mo>)</mo><mo>−</mo><msub><mi>a</mi><mi>min</mi></msub><mo>(</mo><mi>t</mi><mo>)</mo><mo>]</mo><mspace width=".2em"/><mi mathvariant="normal">d</mi><mi>x</mi><mo>,</mo><mspace width="1em"/><mi>D</mi><mo>=</mo><mfrac><mrow><msub><mi>E</mi><mi>C</mi></msub><mo>−</mo><msub><mi>E</mi><mi>I</mi></msub></mrow><mrow><msub><mi>E</mi><mi>C</mi></msub><mo>+</mo><msub><mi>E</mi><mi>I</mi></msub></mrow></mfrac></mrow></math><span class="eqno">(1)</span></div>
<p>Handover is the first challenger amount greater than incumbent amount, equivalently D crosses0. It is distinct from the D≥0.5 replacement threshold, incumbent disappearance, peak-resolution time, or full-angle reorientation. The B ratios0.5,1 and2 correspond to D values−1/3,0 and1/3. A2's final winner is spatial association of the final single region, with field stationarity assessed separately.</p>
<p>Incumbent50% excess loss and challenger reaching50% of basal incumbent excess remain saved native readouts. The original fixed-incumbent/cue-footprint diagnostics are retained in the immutable original statuses; they do not define the corrected handover or B seed ratio. The corrected full-field crosswalk links every original A1 case to its explicit topology and replayed continuous event. The retained A1 continuation gate stops at5,000 chemical units only when its original fixed-incumbent excess is at most10% of basal and a cue-associated single threshold region is present; otherwise it extends to20,000. This historical continuation gate does not define the corrected handover. All roots replay the original saved upper state, check full native species, and retain positivity/conservation evidence.</p>
<table><thead><tr><th>Native cue</th><th>Peak</th><th>Full support</th><th>Chemical cells / unit length</th></tr></thead><tbody><tr><td>Wave-pinning</td><td>0.1</td><td>3</td><td>20</td></tr><tr><td>Goryachev</td><td>6.4</td><td>3</td><td>20</td></tr><tr><td>Holmes3</td><td>1</td><td>3</td><td>20</td></tr><tr><td>Spring</td><td>4</td><td>0.2</td><td>Native material-cell density; actual circumference varies</td></tr></tbody></table>
<p>All12 fresh basals qualified and remained cue-free stable for53,000 chemical units or5,300 Spring units, covering the longest completed A2 protocol. Native positivity and conserved-pool checks pass. The primary A1 ensemble contains152 placements and24 mirrored controls. Of the143 external/borderline placements,109 have a valid continuous two-region handover and34 have an alternate one-region/topology outcome. A2's545 completed pulses give72 bracketed shortest tested switches and37 placements with no observed switch at the tested durations; all34 unavailable calibrations remain explicit. Nine internal references are supplementary.</p>
<p>B contains357 primary runs plus24 mirrored controls. Across all381,307 achieve eligible two-region withdrawal,39 take over during the cue as a single region,28 reach the ratio without separated primary contours, and7 do not attain the target by their completed seed cap. There are52 eligible cue-free resolution lower bounds. Primary-run inclusion by model is shown below; every alternate outcome remains available in the tables and native-history reader.</p>
<table><thead><tr><th>Model</th><th>Eligible primary B seeds</th><th>Alternate or unreached primary seeds</th></tr></thead><tbody>'''+brows+'''</tbody></table>
<details><summary>All seven unattained seed-ratio lower bounds</summary><p>These bounds concern cue-driven seeding, not cue-free resolution time. No cue-free competition time is invented when the required seed state was not attained.</p><table><thead><tr><th>Model</th><th>Size label</th><th>Initial gap</th><th>Target ratio</th><th>Seeding time bound</th></tr></thead><tbody>'''+seed_bounds+'''</tbody></table></details>
<p>At withdrawal the exchange rate is calculated from the exact cue-free right-hand side. The current watershed integration weights and background are frozen for this instantaneous budget, separating reaction and diffusion and, for Spring, dilution and changing spatial measure. This fixed-region derivative is not the material derivative of moving region boundaries. The native RHS arrays, integration weights and term totals remain linked to each B run.</p>
<p>Geometry uses both routes and measured contours. Initial placements additionally retain center distance, target angle, long route and50%-edge reclassification. B retains current10%,25%,50% and75% contour gaps. Peak disappearance uses the declared resolved-peak prominence0.05 and contrast gate, while mere contour coalescence is stored separately. Global headings are unavailable below the declared exact first-moment floor; reporter-gated timing does not unwrap through cancellation. Conservative complex-moment chord guards flag possible gaps between saved headings without claiming a continuously sampled trajectory. Every native history also links the full saved-interval global heading-rate time course: the principal angular increment divided by that interval's duration, with exact-direction and0.05/0.1/0.2 gates. These preserve cancellation gaps and do not assert instantaneous RHS rates or recover unobserved complete turns.</p>
<h2>Numerical sensitivity and limits</h2>
<p>Two models, three gaps, one size each have fresh mesh and time-step controls. Both fine-grid basals independently qualify; all12 A1 controls and36 B controls pass native audits. A1 handover and B withdrawal-time errors are at most0.0449%. All36 B controls retain the seed status, primary-gap eligibility and winner;35 of36 saved resolution brackets overlap. Goryachev C37.8, initial gap4, ratio1 shifts from1040–1060 to1080–1100 units on the finer mesh, about4% in resolution time. This sensitivity is retained explicitly. The24 B mirrored placements retain their outcomes.</p>
<p>Finite full-field stationarity and scalar/spatial fate remain separate. A finite cap supplies a lower bound, not proof that a state is permanently arrested. A2 calibrates only the five declared multipliers, so brackets are tested-duration brackets rather than an exact universal commitment threshold. No noise, dose sweep, moving cue, deBelly ring, or cell migration law enters this figure. Native model differences identify contrasting response classes; the assay does not by itself establish a unique causal molecular pathway.</p>
<h2 id="sources">Complete source tables and independent audits</h2><ul>'''+data_links+'''</ul><ul>'''+audit_links+'''</ul>
<p><a href="supporting/metadata/main_generation.json">Main-figure generation and hashes</a> · <a href="supporting/metadata/native_atlas_audit.json">Complete atlas audit</a> · <a href="supporting/data/native_history_index.json">Every native case and archive hash</a> · <a href="supporting/metadata/reader_generation.json">Reader provenance</a></p>
<p class="small">The immutable failed first discrete-boundary root refinements are retained as diagnostic evidence. The corrected continuous-subcell convention and successful v2 proofs supply the final measurements. Missing producer work is never counted as censoring.</p></main></html>'''
    text=readable_text(text)
    (out/"index.html").write_text(text)
    caption_start=text.index('<h2>Standalone figure legend</h2>')
    caption_end=text.index('<h2>Competition and angular reorientation')
    (out/"supporting/metadata/figure3_distance_caption.html").write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Figure3 legend</title><style>'+CSS+'</style><main>'+text[caption_start:caption_end]+'</main></html>')


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--run",type=Path,default=atlas.DEFAULT)
    parser.add_argument("--output",type=Path,default=OUT);args=parser.parse_args()
    out=args.output.resolve()
    cases,summary=collect(args.run.resolve(),out)
    comparison=controls(out)
    for name in ("generation.json","heading_audit.json"):
        shutil.copy2(SYNTHESIS/name,out/"supporting/metadata"/("synthesis_"+name))
    for name in ("B_control_evidence_and_geometry.csv",):
        shutil.copy2(SYNTHESIS/name,out/"supporting/data"/name)
    atlas_viewer(cases,summary,out)
    reader(out,summary,comparison)
    files=[out/"index.html",out/"supporting/viewers/native_history_viewer.html",out/"supporting/data/native_history_index.json",out/"supporting/data/B_control_convergence.csv",out/"supporting/metadata/figure3_distance_caption.html"]
    write(out/"supporting/metadata/reader_generation.json",dict(generated_utc=summary["generated_utc"],native_integrations=0,
        collector_source_sha256=sha(Path(__file__)),direction_guard_source_sha256=sha(Path(__file__).with_name("signature_figure3_heading_20261020.py")),
        heading_rate_definition="Principal saved-interval angular increment divided by native time interval; exact-direction and0.05/0.1/0.2 coherence/chord guards; NaN gaps preserved. Not an instantaneous RHS derivative or unobserved full-turn count.",
        atlas_summary=summary,outputs={str(p.relative_to(out)):sha(p) for p in files}))
    print(json.dumps(dict(reader=str(out/"index.html"),verified_histories=len(cases),B_controls=len(comparison))))


if __name__=="__main__":main()
