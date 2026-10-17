"""Report exact pre-cue basal checkpoints without integrating or replacing states."""
from pathlib import Path
import csv, hashlib, html, json, textwrap
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / 'data/runs/debelly_corrected_three_test_20261008'
OUT = ROOT / 'figures/initial_polarity_length'
DATA = ROOT / 'data/runs/initial_polarity_length_corrected_20261008'
LABELS = dict(wave_pinning='Wave-pinning', goryachev='Goryachev', debelly='deBelly', holmes_model3='Holmes 3', spring='Spring')
COLORS = ['#0072B2', '#D55E00', '#009E73']
STAMP = datetime.now(timezone.utc).isoformat()
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.titleweight': 'bold', 'pdf.fonttype': 42})

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def csvwrite(path, rows):
    with path.open('w', newline='') as out:
        w = csv.DictWriter(out, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def amount(u, q, width, length, spring=False):
    if spring:
        inside = q < width / length
        xx = np.r_[q[inside], width / length]
        yy = np.r_[u[inside], np.interp(width / length, q, u)]
        return float(np.trapezoid(yy, xx) * length)
    # Native finite-volume values are cell averages. Integrate the fractional edge cell.
    edges = np.linspace(0, length, len(u) + 1)
    return float(u @ np.clip(width - edges[:-1], 0, np.diff(edges)))

def measure(c, fields, q):
    model = c['model']; L0 = c['length']
    L = float(fields['length']) if model == 'spring' else L0
    name = 'rac' if model == 'debelly' else 'rac_active' if model == 'holmes_model3' else 'active'
    u = fields[name]; baseline = float(u.min()); threshold = baseline + .5 * np.ptp(u)
    assert u[0] >= threshold and np.ptp(u) > 0
    j = np.flatnonzero(u < threshold)
    if not len(j): width = L
    else:
        i = int(j[0]); assert i > 0
        edge = q[i-1] + (threshold-u[i-1])/(u[i]-u[i-1])*(q[i]-q[i-1])
        width = float(edge*L)
    total = amount(u, q, width, L, model == 'spring')
    excess = amount(np.maximum(u-baseline, 0), q, width, L, model == 'spring')
    inactive = None if model == 'debelly' else float(fields['rac_inactive' if model == 'holmes_model3' else 'inactive'][-1])
    weights = np.ones(len(q)) / len(q)
    if model == 'spring':
        weights[:] = 1/(len(q)-1); weights[[0,-1]] *= .5
    orientation = float((u*(2*q-1)) @ weights / (u @ weights))
    np.testing.assert_allclose(orientation, c['initial_orientation'], rtol=1e-9, atol=1e-11)
    assert width > 0 and total >= excess >= 0 and np.isfinite(list(fields[name])).all()
    row = dict(model=model, domain_length=L0, physical_length=L,
               old_pole_mean_active=total/width, old_pole_peak_active=float(u.max()),
               old_pole_width=width, old_pole_fraction=width/L, old_pole_total_active=total,
               old_pole_excess_active=excess, opposite_inactive=inactive,
               background_active=baseline, pole_threshold=threshold, active_polarity=orientation,
               length_over_L0=L/L0 if model == 'spring' else None,
               spring_elastic_tension=.01*(L-L0) if model == 'spring' else None,
               spring_extension_feedback=80*(L-L0) if model == 'spring' else None,
               mca_bound_polarity=float(fields['mca_bound'] @ (2*q-1) / fields['mca_bound'].sum()) if model == 'debelly' else None,
               checkpoint=c['checkpoint'], checkpoint_sha256=c['checkpoint_sha256'])
    np.testing.assert_allclose(total-excess, baseline*width, rtol=1e-10, atol=1e-12)
    return row

def main():
    OUT.mkdir(parents=True, exist_ok=True); DATA.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((CAMPAIGN/'manifest.json').read_text())
    status = json.loads((CAMPAIGN/'results/status.json').read_text()); assert status['complete'] and status['execution_errors'] == 0
    states = {}; rows = []; excluded = {}; selection = {}
    for model in manifest['models']:
        cases = sorted([c for c in manifest['cases'] if c['model']==model and c['condition']=='fixed_relative'], key=lambda c:c['length'])
        excluded[model] = [c['length'] for c in cases if not c['eligible']]
        selected = [c for c in cases if c['eligible']]
        selection[model] = [selected[i]['length'] for i in [0, (len(selected)-1)//2, len(selected)-1]]
        for c in selected:
            p = ROOT/c['checkpoint']; assert sha(p) == c['checkpoint_sha256']
            siblings = [s for s in manifest['cases'] if s['model']==model and s['length']==c['length']]
            assert len(siblings)==3 and all(s['checkpoint_sha256']==c['checkpoint_sha256'] for s in siblings)
            with np.load(p) as z:
                fields = {k[6:]:z[k][-1].copy() for k in z.files if k.startswith('field_')}
                q = np.linspace(0,1,len(z['x'])) if model=='spring' else z['x']/c['length']
                np.testing.assert_allclose(q, np.linspace(0,1,len(q)) if model=='spring' else (np.arange(len(q))+.5)/len(q), atol=1e-12)
                assert all(np.isfinite(v).all() for v in fields.values())
                states[model,c['length']] = (q,fields)
                rows.append(measure(c,fields,q))
    assert len(rows)==97
    csvwrite(DATA/'steady_state_metrics.csv', rows)
    snapshots=[]
    for model, lengths in selection.items():
        for L in lengths:
            q,fields=states[model,L]
            for name,value in fields.items():
                for fraction, v in zip(q, np.broadcast_to(value, q.shape)):
                    snapshots.append(dict(model=model, domain_length=L, position_fraction=float(fraction), field=name, value=float(v), spatial_field=bool(np.ndim(value)>0)))
    csvwrite(DATA/'selected_field_profiles.csv', snapshots)
    pages=[]
    pdfpath=OUT/'initial_polarity_length.pdf'
    with PdfPages(pdfpath) as pdf:
        def save(fig,title,caption,slug):
            fig.suptitle(title,fontsize=19,weight='bold',y=.975)
            caption_lines = textwrap.wrap(caption,145)
            caption_y = .045 + (len(caption_lines)*15/72 + .12)/fig.get_figheight()
            fig.text(.065,caption_y,'\n'.join(caption_lines),fontsize=10,va='top',linespacing=1.5)
            fig.text(.065,.018,'Initial polarized state · verified basal release · before the right cue',fontsize=9,color='#526674')
            fig.text(.95,.018,str(len(pages)+1),ha='right',fontsize=9)
            fig.savefig(OUT/f'{slug}.png',dpi=140)
            pdf.savefig(fig); plt.close(fig)
            pages.append(dict(title=title,caption=caption,image=f'{slug}.png'))
        fig=plt.figure(figsize=(12,10))
        text=[('The initial polarized state as cell length changes',
               'Five models; exact steady pre-cue states used in the current three-test comparison. No additional stimulus is present in these snapshots.'),
              ('Which states are shown?',
               'Each state was prepared by native rest, a finite left patch at nominal amplitude with 20% support, then verified basal release. Preparation required two stable 1,000-unit windows and a 10,000-unit holdout. These are the exact checkpoints shared by all three subsequent right-cue tests. Summaries include every eligible length; snapshots use the shortest, middle and longest eligible sampled lengths.'),
              ('What counts as the old pole?',
               'The old pole is the boundary-connected left region whose active field exceeds halfway between its minimum and maximum. Its edge is linearly interpolated between samples. Mean concentration and total material use this same region; the peak is also shown. Solid amount curves include the background; dashed curves subtract the minimum active concentration. Width is in native length units, not a fraction of the cell.'),
              ('Coordinates and fields',
               'Snapshot x axes show position as a fraction of the cell (0 = old pole, 1 = opposite end); field values are not rescaled. Domain length on summary x axes is the native prescribed domain, or resting L0 for Spring. Spring pole widths and amounts use its steady physical length. deBelly amounts use native nominal finite-volume coordinates, consistent with its chemical solver; mechanical deformation is shown separately. Amounts are one-dimensional concentration integrals, not calibrated molecule counts.'),
              ('Model-specific readouts',
               'Active means Rac for Holmes and deBelly. Opposite-end inactive means the final native sample: the right boundary node for Spring and the rightmost cell average for finite-volume models. deBelly has no explicit inactive Rac pool; Rho and unbound MCA are separate variables and are not substitutes. Its inactive panel is therefore unavailable. MCA polarity is the normalized first spatial moment of bound attachments (negative = toward the old pole). Spring shows physical L/L0 plus the native elastic term 0.01*(L-L0).'),
              ('Coverage and provenance',
               '97 polarized preparations: Wave-pinning 20, Goryachev 21, deBelly 12, Holmes 3 23, Spring 21. The 13 ineligible preparations remain excluded and are marked on summary axes. All checkpoint hashes and shared-start identities were checked. Existing preparation and response data are unchanged. Report generated '+STAMP+'. Simulation producer finished; source campaign finalized '+status['generated_at']+'.')]
        y=.94
        for i,(heading,body) in enumerate(text):
            if i==0: continue
            fig.text(.065,y,heading,fontsize=14,weight='bold',va='top'); y-=.034
            lines=textwrap.wrap(body,130)
            fig.text(.065,y,'\n'.join(lines),fontsize=11,va='top',linespacing=1.5)
            y-=len(lines)*.023+.038
        assert y>.09,y
        save(fig,text[0][0],text[0][1],'01_definitions')
        for page, metrics in enumerate([
                [('old_pole_mean_active','Active concentration under old pole','Concentration'),('old_pole_width','Old-pole width','Native length')],
                [('old_pole_total_active','Total active material under old pole','Concentration × length'),('opposite_inactive','Inactive at opposite end','Concentration')]],2):
            fig,axes=plt.subplots(5,2,figsize=(12,15))
            for i,model in enumerate(manifest['models']):
                rs=[r for r in rows if r['model']==model]; x=[r['domain_length'] for r in rs]
                for ax,(key,title,unit) in zip(axes[i],metrics):
                    if key=='opposite_inactive' and model=='debelly':
                        ax.text(.5,.5,'Unavailable: no inactive Rac field',ha='center',transform=ax.transAxes,fontsize=11)
                        ax.set_xlim(min(x),max(x))
                    else:
                        ax.plot(x,[r[key] for r in rs],'-o',color=COLORS[0],ms=3,label='Pole mean' if key=='old_pole_mean_active' else 'Total' if key=='old_pole_total_active' else None)
                        extra='old_pole_peak_active' if key=='old_pole_mean_active' else 'old_pole_excess_active' if key=='old_pole_total_active' else None
                        if extra:
                            ax.plot(x,[r[extra] for r in rs],'--',color=COLORS[1],label='Peak' if extra=='old_pole_peak_active' else 'Excess above minimum')
                            ax.legend(frameon=False,fontsize=9)
                        ax.set_ylim(bottom=0)
                    if excluded[model]:
                        ax.scatter(excluded[model],[-.04]*len(excluded[model]),marker='x',c='#888888',s=23,transform=ax.get_xaxis_transform(),clip_on=False)
                        lo,hi=ax.get_xlim(); ax.set_xlim(min(lo,min(excluded[model])-.5),hi)
                    ax.set_title(LABELS[model]+' · '+title,fontsize=12)
                    ax.set_xlabel('Resting domain length L0' if model=='spring' else 'Domain length');ax.set_ylabel(unit);ax.grid(alpha=.18)
            fig.subplots_adjust(left=.09,right=.96,top=.93,bottom=.15,hspace=.72,wspace=.3)
            save(fig,'Steady pre-cue state: '+('pole concentration and width' if page==2 else 'pole material and opposite-end inactive'),
                 'Each marker is one verified polarized checkpoint. Gray crosses below the x axis mark lengths without an eligible polarized preparation; these are not zero-valued observations. Lines connect sampled lengths within each model. Native units differ between models. Pole width is the connected half-contrast region; concentration and material use that same region.',f'{page:02d}_summaries')
        for model in manifest['models']:
            lengths=selection[model]; names=list(states[model,lengths[0]][1])
            for start in range(0,len(names),8):
                chunk=names[start:start+8]; nrows=(len(chunk)+1)//2
                fig,axes=plt.subplots(nrows,2,figsize=(12,max(8,3*nrows+2)),squeeze=False)
                for ax,name in zip(axes.flat,chunk):
                    for color,L in zip(COLORS,lengths):
                        q,fields=states[model,L]; value=fields[name]
                        ax.plot(q,np.broadcast_to(value,q.shape),color=color,lw=2,label=f'{"L0" if model=="spring" else "L"} = {L:g}')
                    ax.set(xlabel='Scaled position along cell',ylabel='Field value (native units)',xlim=(0,1))
                    ax.set_title(name.replace('_',' ')+(' · global scalar' if np.ndim(value)==0 else ''),fontsize=12)
                    ax.grid(alpha=.18);ax.legend(frameon=False,fontsize=9)
                    ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,4),useOffset=False)
                for ax in list(axes.flat)[len(chunk):]:ax.set_visible(False)
                fig.subplots_adjust(left=.09,right=.96,top=.89,bottom=.19,hspace=.65,wspace=.3)
                caption='Shortest, middle and longest eligible sampled lengths share each field facet; blue, orange and green identify length, not time. All fields stored in the exact checkpoint are included, including derived rates and scalar mechanical diagnostics. Global scalars are drawn as horizontal lines for comparison; they are not spatial fields.'
                if model=='debelly':caption+=' Displacement is the stored native displacement; scalar boundary separation is its domain-extension diagnostic. These are basal rates before any right cue.'
                if model=='spring':caption+=' Length is the steady physical domain length; snapshot position runs from the current left to right boundary.'
                save(fig,LABELS[model]+': pre-cue field snapshots'+(f' ({start//8+1})' if len(names)>8 else ''),caption,f'{len(pages)+1:02d}_{model}_fields_{start//8+1}')
        fig,axes=plt.subplots(2,2,figsize=(12,10))
        rs=[r for r in rows if r['model']=='spring'];x=[r['domain_length'] for r in rs]
        for ax,key,title,ylabel in [(axes[0,0],'length_over_L0','Spring: steady extension ratio','Physical L / resting L0'),(axes[0,1],'spring_elastic_tension','Spring: native global elastic term','0.01 × (L - L0)'),(axes[1,0],'spring_extension_feedback','Spring: extension-dependent deactivation','80 × (L - L0)')]:
            ax.plot(x,[r[key] for r in rs],'-o',color=COLORS[0],ms=4);ax.set(title=title,xlabel='Resting domain length L0',ylabel=ylabel);ax.grid(alpha=.18)
        rs=[r for r in rows if r['model']=='debelly'];ax=axes[1,1]
        ax.plot([r['domain_length'] for r in rs],[r['mca_bound_polarity'] for r in rs],'-o',color=COLORS[0],ms=4,label='Bound MCA')
        ax.axhline(0,color='#999999',lw=1);ax.set(title='deBelly: MCA polarity',xlabel='Domain length',ylabel='Normalized first moment');ax.grid(alpha=.18)
        fig.subplots_adjust(left=.1,right=.96,top=.9,bottom=.22,hspace=.5,wspace=.35)
        save(fig,'Steady pre-cue mechanical and attachment state',
             'Spring: the elastic restoring term is spatially global, with equal and opposite signs in the boundary velocities. It is a native model force/tension proxy, not a calibrated physical tension. The biochemical feedback uses a different coefficient: deactivation is 3 + 80*(L-L0). Separate panels preserve both native quantities without a second-axis rescaling. deBelly: MCA polarity integrates bound attachments on the native nominal cell coordinate; negative values indicate enrichment toward the original left pole.',f'{len(pages)+1:02d}_mechanics')
    links=''.join(f'<section><h2>{html.escape(p["title"])}</h2><img src="{p["image"]}" alt="{html.escape(p["title"])}"><p>{html.escape(p["caption"])}</p></section>' for p in pages)
    (OUT/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Initial polarity versus length</title><style>body{max-width:1150px;margin:36px auto;padding:0 24px;font:17px/1.6 system-ui;color:#202c35}img{max-width:100%}section{border-top:1px solid #ddd;margin-top:32px}a{color:#176487}</style><h1>Initial polarized state versus length</h1><p><a href="initial_polarity_length.pdf">PDF report</a> · <a href="../../data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv">All summary values (CSV)</a> · <a href="../../data/runs/initial_polarity_length_corrected_20261008/selected_field_profiles.csv">Selected field profiles (CSV)</a></p><p>Generated '+STAMP+'; source simulation producer finished.</p>'+links+'</html>')
    provenance=dict(generated_at=STAMP,producer='finished',source_campaign_generated_at=status['generated_at'],source_manifest_sha256=sha(CAMPAIGN/'manifest.json'),renderer_sha256=sha(Path(__file__)),pdf_sha256=sha(pdfpath),pages=len(pages),qualified_preparations=len(rows),selected_lengths=selection,excluded_lengths=excluded,all_checkpoint_hashes_pass=True,all_three_geometry_checkpoint_identities_pass=True,independent_active_moments_pass=True,pole_amount_background_identity_pass=True,summary_csv_sha256=sha(DATA/'steady_state_metrics.csv'),snapshots_csv_sha256=sha(DATA/'selected_field_profiles.csv'),checkpoints=[dict(model=r['model'],length=r['domain_length'],path=r['checkpoint'],sha256=r['checkpoint_sha256']) for r in rows])
    (DATA/'report_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps(dict(pdf=str(pdfpath),pages=len(pages),checkpoints=len(rows),selected_lengths=selection),indent=2))

if __name__=='__main__':main()
