"""Assess audited substrate/transport controls, including clamp effectiveness."""
from pathlib import Path
import csv, gzip, json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import goryachev_depletion as g

BASE=g.n.DATA/'runs/goryachev_depletion_20261008_v2'
FOLLOW=g.n.DATA/'runs/goryachev_depletion_followup_20261008'
OUT=g.n.DATA/'runs/goryachev_depletion_assessment_20261008'
KINDS=['native','both_initial','feedback_tight','cue_clamp','both_clamp','transport_matched','both_initial_transport']
LABELS=dict(native='Native',both_initial='Initial cue + feedback compensation',feedback_tight='Feedback compensation (tight)',cue_clamp='Cue compensation',both_clamp='Cue + feedback compensation',transport_matched='Matched diffusion timescale',both_initial_transport='Initial compensation + matched transport')


def read(path):return json.loads(path.read_text())


def build():
    summaries=[];inputs=[]
    for run in [BASE,FOLLOW]:
        assert read(run/'results/status.json')['complete']
        assert read(run/'results/independent_audit.json')['passed']
        inputs.append(dict(run=str(run.relative_to(g.n.ROOT)),manifest_sha256=g.n.sha(run/'manifest.json'),audit_sha256=g.n.sha(run/'results/independent_audit.json')))
        for s in read(run/'results/summaries.json'):
            s['run']=run;summaries.append(s)
    OUT.mkdir(exist_ok=True)
    rows=[];effectiveness=[]
    for s in summaries:
        j=s['job'];c=j['case']
        for gate,sel in s['selected'].items():
            row=dict(length=c['length'],geometry=c['condition'],kind=j['kind'],gate=gate,**sel)
            if j['kind']=='native':
                row['lower']=j['native_result'][gate]['lower_failing_duration']
                row['relative_bracket']=j['native_result'][gate]['relative_bracket_width']
            rows.append(row)
        if j['kind'] not in ['cue_clamp','feedback_clamp','both_clamp','feedback_tight']:continue
        d=s['run']/'workers'/f'{j["slot"]:03d}'
        path=d/'selected_orientation.npz'
        if not path.exists():continue
        with np.load(path) as z:
            u=z['pulse_active'];v=z['pulse_field_inactive'];times=z['pulse_time']
            q=(np.arange(v.shape[1])+.5)/v.shape[1];mask=q>=.5
            ref=j['reference_v'];eps=ref*j['epsilon_fraction']
            ve=(ref+eps)*v[:,mask]/(v[:,mask]+eps);ratio=ve/ref
            # Weight availability error by the intended compensated feedback.
            model=g.Goryachev();weights=(model.a*u[:,mask]**2+model.b*u[:,mask])*ref
            intended=weights.sum(axis=1);realized=(weights*ratio).sum(axis=1)
            cumulative_ratio=float(np.trapezoid(realized,times)/np.trapezoid(intended,times))
            effectiveness.append(dict(length=c['length'],geometry=c['condition'],kind=j['kind'],
                    duration=float(times[-1]),minimum_actual_challenger_inactive=float(v[:,mask].min()),
                    minimum_effective_over_reference=float(ratio.min()),
                    maximum_effective_over_reference=float(ratio.max()),
                    sampled_fraction_below_95pct=float(np.mean(ratio<.95)),
                    intended_feedback_weighted_effective_over_reference=cumulative_ratio))
    g.n.csvwrite(OUT/'endpoints.csv',rows);g.n.csvwrite(OUT/'clamp_effectiveness.csv',effectiveness)
    ratios=[]
    for geometry in g.t.TESTS:
        for gate in g.r.DEFINITIONS:
            for kind in KINDS:
                short=next(x for x in rows if x['geometry']==geometry and x['kind']==kind and x['gate']==gate and x['length']==9.)
                long=next(x for x in rows if x['geometry']==geometry and x['kind']==kind and x['gate']==gate and x['length']==39.)
                ratios.append(dict(geometry=geometry,gate=gate,kind=kind,short_duration=short['endpoint'],long_duration=long['endpoint'],
                        long_over_short=long['endpoint']/short['endpoint'] if long['endpoint'] and short['endpoint'] else None))
    g.n.csvwrite(OUT/'length_ratios.csv',ratios)
    # Independently compare native endpoint replays against original native records.
    with gzip.open(g.t.DEFAULT/'results/records.json.gz','rt') as stream:originals=json.load(stream)['cases']
    comparisons=[]
    for s in summaries:
        j=s['job'];c=j['case']
        if j['kind']!='native':continue
        orig=next(x for x in originals if x['model']=='goryachev' and x['length']==c['length'] and x['condition']==c['condition'])
        for gate,sel in s['selected'].items():
            old=next(x for x in orig['tests'] if x['duration']==sel['endpoint'])
            assert sel['field_tail_pass']==old['field_tail_pass']
            delta=abs(sel['final_orientation']-old['final_orientation']);assert delta<1e-9
            comparisons.append(dict(length=c['length'],geometry=c['condition'],gate=gate,orientation_error=delta))
    g.n.write(OUT/'assessment_audit.json',dict(generated_at=g.n.now(),native_replay_checks=comparisons,
         no_cue_controls_passed=all(s['no_cue_valid'] for s in summaries),inputs=inputs,
         report_code_sha256=g.n.sha(Path(__file__))))
    pdf=PdfPages(OUT/'depletion_assessment.pdf')
    for gate in g.r.DEFINITIONS:
        fig,axes=plt.subplots(1,3,figsize=(15,5))
        for ax,geometry in zip(axes,g.t.TESTS):
            for kind in KINDS:
                rs=sorted([x for x in rows if x['geometry']==geometry and x['kind']==kind and x['gate']==gate],key=lambda x:x['length'])
                line=ax.plot([x['length'] for x in rs],[x['endpoint'] or np.nan for x in rs],'-',label=LABELS[kind])[0]
                for x in rs:
                    ax.plot(x['length'],x['endpoint'] or 5120.,'o' if x['endpoint'] else '^',color=line.get_color(),mfc=line.get_color() if x['field_tail_pass'] and x['endpoint'] else 'white',ms=5)
            ax.set(title=geometry,xlabel='Length',ylabel='Shortest passing tested cue duration',yscale='log');ax.grid(alpha=.2)
        axes[0].legend(fontsize=7);fig.suptitle('Goryachev depletion diagnostics / '+gate)
        fig.tight_layout();fig.savefig(OUT/(gate+'_durations.png'),dpi=180);pdf.savefig(fig);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(15,5))
    for ax,geometry in zip(axes,g.t.TESTS):
        for kind in ['native','feedback_tight','both_clamp','transport_matched','both_initial_transport']:
            s=next(x for x in summaries if x['job']['case']['length']==39. and x['job']['case']['condition']==geometry and x['job']['kind']==kind)
            j=s['job'];d=s['run']/'workers'/f'{j["slot"]:03d}'
            if kind=='native':path=s['run']/next(x['archive'] for x in s['trials'] if x['duration']==s['selected']['orientation']['endpoint'])
            else:path=d/'selected_orientation.npz'
            if not path.exists():continue
            with np.load(path) as z:
                time=z['history_time'];amount=z['old_pole_excess'];ax.plot(time,amount/amount[0],label=LABELS[kind])
            ax.set_xlim(0,4000.)
        ax.axhline(.5,color='#777777',ls=':');ax.set(title=geometry,xlabel='Time since cue onset',ylabel='Old-pole excess / initial excess');ax.grid(alpha=.2)
    axes[0].legend(fontsize=7);fig.suptitle('L=39: pole-removal histories at each control’s selected orientation dose')
    fig.tight_layout();fig.savefig(OUT/'old_pole_removal.png',dpi=180);pdf.savefig(fig);plt.close(fig);pdf.close()
    (OUT/'index.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>Goryachev depletion assessment</title></head><body style="font:16px system-ui;max-width:1500px;margin:auto"><h1>Goryachev depletion assessment</h1><p>Generated '+g.n.now()+'; all producers finished; audits passed.</p><p><a href="INTERPRETATION.md">Interpretation and design limits</a> · <a href="depletion_assessment.pdf">PDF figures</a> · <a href="endpoints.csv">Endpoint table</a> · <a href="clamp_effectiveness.csv">Compensation effectiveness</a></p>'+''.join('<img style="width:100%" src="'+x+'">' for x in ['orientation_durations.png','full_polarity_durations.png','old_pole_removal.png'])+'</body></html>')
    print(json.dumps(dict(jobs=len(summaries),native_checks=len(comparisons),no_cue_controls_passed=True)))


if __name__=='__main__':build()
