"""Native active/inactive snapshots for relative-20% Goryachev reversal."""
from pathlib import Path
import csv, json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import goryachev_depletion as g

SOURCE=g.n.DATA/'runs/goryachev_depletion_20261008_v2'
DATA=g.n.DATA/'runs/goryachev_native_snapshots_20261008'
FIG=g.n.ROOT/'figures/goryachev_native_snapshots'


def main():
    status=g.r.read(SOURCE/'results/status.json')
    assert status['complete'] and status['audit_passed']
    summaries=g.r.read(SOURCE/'results/summaries.json')
    DATA.mkdir(exist_ok=True);FIG.mkdir(exist_ok=True)
    selected=[];rows=[]
    with PdfPages(FIG/'active_inactive_snapshots.pdf') as pdf:
        for L in [9.,39.]:
            s=next(s for s in summaries if s['job']['case']['length']==L and s['job']['case']['condition']=='fixed_relative' and s['job']['kind']=='native')
            c=s['job']['case'];D=s['selected']['orientation']['endpoint']
            record=next(r for r in s['trials'] if r['duration']==D)
            path=SOURCE/record['archive'];assert g.n.sha(path)==record['archive_sha256']
            assert g.n.sha(g.n.ROOT/c['checkpoint'])==c['checkpoint_sha256']
            with np.load(path) as z, np.load(g.n.ROOT/c['checkpoint']) as cp:
                np.testing.assert_array_equal(z['originalstate'],cp['state'])
                np.testing.assert_allclose(z['pulse_field_active'][0],cp['state'][0],rtol=0,atol=1e-12)
                np.testing.assert_allclose(z['pulse_field_inactive'][0],cp['state'][1],rtol=0,atol=1e-12)
                snaps=[]
                for label,phase,time in [('Before cue','pulse',0.),('Early cue','pulse',2.5),
                     ('25% of pulse','pulse',D/4),('50% of pulse','pulse',D/2),
                     ('75% of pulse','pulse',3*D/4),('Withdrawal','pulse',D),
                     ('50 after withdrawal','recovery',50.),('Final recovery','recovery',10000.)]:
                    i=int(np.argmin(abs(z[phase+'_time']-time)))
                    actual=float(z[phase+'_time'][i]);elapsed=actual+(D if phase=='recovery' else 0.)
                    u=z[phase+'_field_active'][i];v=z[phase+'_field_inactive'][i]
                    mass=(u+v).sum();assert abs(mass/cp['state'].sum()-1.)<1e-6
                    snaps.append(dict(label=label,phase=phase,time=elapsed,u=u.copy(),v=v.copy(),
                                      orientation=float(z[phase+'_orientation'][i])))
                    for x,active,inactive in zip(z['x'],u,v):
                        rows.append(dict(length=L,cue_duration=D,stage=label,phase=phase,
                                         time_since_cue_onset=elapsed,x=float(x),
                                         active=float(active),inactive=float(inactive)))
                fig,axes=plt.subplots(2,8,figsize=(20,6.1),sharex=True,sharey='row')
                umax=max(float(p['u'].max()) for p in snaps)*1.07
                for col,p in enumerate(snaps):
                    for row,key,color in [(0,'u','#0072B2'),(1,'v','#D55E00')]:
                        ax=axes[row,col];ax.plot(z['x'],p[key],color=color,lw=2)
                        if p['phase']=='pulse' and col>0:ax.axvspan(.8*L,L,color='#777777',alpha=.13,zorder=0)
                        ax.set_xlim(0,L);ax.set_ylim(0,umax if row==0 else .62)
                        ax.set_xticks([0,L/2,L]);ax.grid(alpha=.18);ax.tick_params(labelsize=10)
                        if row==1:ax.set_xlabel('Position x',fontsize=10)
                    axes[0,col].set_title(p['label']+'\nt = '+f'{p["time"]:g}',fontsize=10,pad=10)
                    axes[1,col].text(.5,-.38,'Polarity = '+f'{p["orientation"]:.3f}',ha='center',transform=axes[1,col].transAxes,fontsize=10)
                axes[0,0].set_ylabel('Active Cdc42\nConcentration',fontsize=12)
                axes[1,0].set_ylabel('Inactive Cdc42\nConcentration',fontsize=12)
                fig.suptitle(f'Native Goryachev reversal · L = {L:g} · right 20% patch · cue duration {D:g}',fontsize=17,y=.98)
                fig.text(.055,.045,'Gray band: cue support while cue is on; peak 6.4. Withdrawal snapshot is immediately before removal.\n'
                    'Actual native concentration and position; no field normalization. Active scale fixed across stages within each cell; inactive scale 0–0.62 in both cells.',fontsize=10,color='#43515e')
                fig.subplots_adjust(left=.055,right=.99,top=.79,bottom=.23,hspace=.28,wspace=.13)
                stem=f'active_inactive_L{int(L)}';fig.savefig(FIG/(stem+'.png'),dpi=170)
                fig.savefig(FIG/(stem+'.svg'));pdf.savefig(fig);plt.close(fig)
                selected.append(dict(length=L,duration=D,archive=str(path.relative_to(g.n.ROOT)),
                     archive_sha256=g.n.sha(path),checkpoint=c['checkpoint'],checkpoint_sha256=c['checkpoint_sha256'],
                     orientation_pass=record['orientation_success'],strength_pass=record['full_polarity_success'],
                     field_tail_pass=record['field_tail_pass'],times=[p['time'] for p in snaps]))
    g.n.csvwrite(DATA/'snapshots.csv',rows)
    stamp=g.n.now()
    g.n.write(DATA/'provenance.json',dict(generated_at=stamp,producer='finished',selected=selected,
        exact_restart_state_pass=True,initial_sample_arrays_atol_1e12_pass=True,selected_snapshot_conservation_pass=True,
        renderer_sha256=g.n.sha(Path(__file__)),csv_sha256=g.n.sha(DATA/'snapshots.csv'),
        pdf_sha256=g.n.sha(FIG/'active_inactive_snapshots.pdf')))
    (FIG/'index.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>Native Goryachev active/inactive snapshots</title></head><body style="font:16px system-ui;max-width:1800px;margin:auto"><h1>Native Goryachev reversal: active and inactive fields</h1><p>Generated '+stamp+'; producer finished. L=9 and L=39, relative-20% cue, exact native shortest passing pulses. No compensation or transport intervention.</p><p><a href="active_inactive_snapshots.pdf">Two-page PDF</a> · <a href="../../data/runs/goryachev_native_snapshots_20261008/snapshots.csv">Snapshot values</a></p>'+''.join('<img style="width:100%" src="active_inactive_L'+str(L)+'.png">' for L in [9,39])+'</body></html>')
    print(json.dumps(dict(figures=str(FIG),selected=selected),indent=2))


if __name__=='__main__':main()
