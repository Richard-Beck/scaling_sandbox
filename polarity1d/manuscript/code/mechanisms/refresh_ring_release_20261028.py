"""Refresh the R1/R2 revision's local release receipts and storage inventory.

Only the ring-r1r2 section of DATA_MANIFEST is replaced; no data is deleted.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess
from xml.etree import ElementTree as ET
from PIL import Image
import audit_ring_revision1_20261028 as A

ROOT=A.R1.ROOT;BASE=A.BASE;META=A.META
START='<!-- ring-r1r2-20261027 begin -->';STOP='<!-- ring-r1r2-20261027 end -->'
def now():return datetime.now(timezone.utc).isoformat()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run(args):return subprocess.run(args,cwd=ROOT,text=True,capture_output=True,check=True).stdout.strip()
def write(p,v):A.R1.write(p,v)
def refresh():
    audit=A.audit();groups=[]
    paths=['polarity1d/manuscript/data/runs/ring_r1_fast_cytosol_20261027','polarity1d/manuscript/data/runs/ring_r2_per_cell_20261027','polarity1d/manuscript/data/runs/poisson_ring_modes_20261027',str(BASE.relative_to(ROOT))]
    measured=now()
    for path in paths:
        n=int(run(['du','-sx','-B1',path]).split()[0]);groups.append(dict(path=path,exists=True,allocated_bytes=n,allocated_gib=n/2**30))
    inventory=dict(measured_utc=measured,measurement='du -sx -B1 (allocated bytes)',atomic=False,groups=groups,filesystem=run(['df','-B1',str(BASE)]),personal_quota='No new personal-quota measurement; filesystem available space is not personal quota headroom.')
    write(META/'storage_measurement.json',inventory)
    state='Both event producers and their extension checks are complete.' if audit['producers_complete'] else 'R1 is complete; existing event producers/extension checks remain active.'
    purposes=[
        'Completed 5,760 matched-noise fast-cytosol histories and 16 native solver checks, exact inputs, immutable sources and integrity receipts. Scientifically authoritative compact records; no new simulations in revision 1.',
        'Existing 144-history per-cell event-supply arm, sharing original middle sizes. Native-cadence polarity, event logs, validated low modes and restart evidence are retained. Producer completion is separate from rate censoring.',
        'Authoritative original 432-history event campaign, including full retained Fourier records, exact frozen inputs and continuation checkpoints. Its existing extension controller is preserved; revision 1 only reads these data.',
        'R1/R2 main PNG/PDF/SVG, complete compact plotted data, revised legends, methods, bootstrap sufficient statistics and integrity receipts. Same figure paths; every original panel retained.'
    ]
    block=START+'\n\n## Ring figures R1 and R2, revision 1\n\n'+f'Measured **{measured}** with `du -sx -B1` (allocated storage; non-atomic traversal). {state} Preserve all authoritative scientific records; inventory does not authorize deletion.\n\n| Path / related group | Measured allocated storage | Purpose and scientific importance/status |\n| --- | ---: | --- |\n'
    for g,purpose in zip(groups,purposes):block+=f'| `{g["path"]}/` | {g["allocated_bytes"]:,} bytes ({g["allocated_gib"]:.6f} GiB) | {purpose} |\n'
    block+='\n[Measurement receipt](polarity1d/manuscript/figures/length_mechanisms_20261008/dev/ring_r1r2_20261027/supporting/metadata/storage_measurement.json) · [Retention policy](RETENTION.md#ring-figures-r1-and-r2-27-october-2026). Existing campaign forecasts remain in their frozen manifests.\n\n'+STOP
    manifest=ROOT/'DATA_MANIFEST.md';text=manifest.read_text();a=text.index(START);b=text.index(STOP,a)+len(STOP);manifest.write_text(text[:a]+block+text[b:])
    artifacts=[];formats=[];links=0
    for fig in ('R1','R2'):
        folder=BASE/f'figure{fig}'
        for path in sorted(folder.rglob('*')):
            if path.is_file() and (path.parent==folder or path.parent==folder/'supporting/data' or path.name in ('main_legend.md','scientific_readme.md')):
                artifacts.append(dict(path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=sha(path)))
        for path in (folder/'README.md',folder/'supporting/metadata/main_legend.md',folder/'supporting/metadata/scientific_readme.md'):
            for url in re.findall(r'\]\(([^)]+)\)',path.read_text()):
                if url.startswith(('http:','https:','#')):continue
                assert (path.parent/url).exists(),(path,url);links+=1
        svg=ET.parse(folder/f'figure{fig}_main.svg').getroot();sw=float(svg.attrib['width'].removesuffix('pt'))*25.4/72
        info=run(['pdfinfo',str(folder/f'figure{fig}_main.pdf')]);pw=float(re.search(r'Page size:\s+([\d.]+)',info).group(1))*25.4/72
        with Image.open(folder/f'figure{fig}_main.png') as im:shape=list(im.size)
        assert abs(sw-180)<.01 and abs(pw-180)<.01
        formats.append(dict(figure=fig,svg_width_mm=sw,pdf_width_mm=pw,png_pixels=shape,width_180mm_pass=True))
    protected=['current manuscript figs','manuscript_outline','polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round7_20261026']
    changes=run(['git','diff','--name-only','a395934f','--',*protected]).splitlines()
    assert not changes,changes
    write(META/'release_audit.json',dict(generated_utc=now(),base_commit='88365c29',preserved_successor='a395934f',branch=run(['git','branch','--show-current']),current_commit=run(['git','rev-parse','HEAD']),protected_scope_changes=changes,artifacts=artifacts,format_checks=formats,local_links_checked=links,structural_checks_pass=True,scientific_audit='revision1_audit.json',producers_complete=audit['producers_complete'],qualification='Structural and numerical checks do not substitute for visual review or campaign completion.'))
    write(META/'completion_audit.json',dict(generated_utc=now(),goal_complete=False,scientific_work_complete=audit['producers_complete'],checkpoint_only=not audit['producers_complete'],revision=1,requirements=[dict(requirement='Every R1 A–H / R2 A–G panel retained; requested encodings and estimands',status='verified',evidence='revision1_audit.json; main_legend.md; PNG/PDF/SVG'),dict(requirement='No retuning, new simulation campaigns or frozen worker changes',status='verified',evidence='Only saved compact seeds, Fourier records and event logs reanalysed'),dict(requirement='Pending event producers and extension checks finish; rebuild again',status='verified' if audit['producers_complete'] else 'pending',evidence='Existing dependency finalizer 21391427 follows controller successors; reads the revised renderer'),dict(requirement='Commit and push',status='see integration_commit.json and final.txt')],scientific_checks=len(audit['checks']),scientific_checks_passed=audit['passed']))
    print(json.dumps(dict(storage_groups=len(groups),artifacts=len(artifacts),formats=formats,protected_scope_changes=changes,producers_complete=audit['producers_complete'])))
if __name__=='__main__':refresh()
