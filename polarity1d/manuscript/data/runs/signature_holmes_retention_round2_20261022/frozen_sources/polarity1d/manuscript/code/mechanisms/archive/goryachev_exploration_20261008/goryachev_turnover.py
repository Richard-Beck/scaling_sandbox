"""Gross activation turnover of native steady incumbent poles versus length."""
from pathlib import Path
import csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import goryachev_depletion as g

DATA=g.n.DATA/'runs/goryachev_incumbent_turnover_20261008'
FIG=g.n.ROOT/'figures/goryachev_incumbent_turnover'


def integral(field,width,length):
    dx=length/len(field)
    overlap=np.clip(width-np.arange(len(field))*dx,0.,dx)
    return float(field@overlap)


def main():
    manifest=g.r.read(g.t.DEFAULT/'manifest.json')
    status=g.r.read(g.t.DEFAULT/'results/status.json');assert status['complete'] and status['execution_errors']==0
    cases=sorted([c for c in manifest['cases'] if c['model']=='goryachev' and c['condition']=='fixed_relative' and c['eligible']],key=lambda c:c['length'])
    DATA.mkdir(exist_ok=True);FIG.mkdir(exist_ok=True)
    model=g.Goryachev();rows=[]
    for c in cases:
        L=c['length'];path=g.n.ROOT/c['checkpoint'];assert g.n.sha(path)==c['checkpoint_sha256']
        with np.load(path) as z:u,v=z['state'].copy()
        q=(np.arange(len(u))+.5)/len(u)
        threshold=float(u.min()+.5*np.ptp(u));j=int(np.flatnonzero(u<threshold)[0]);assert j>0
        width=float((q[j-1]+(threshold-u[j-1])/(u[j]-u[j-1])*(q[j]-q[j-1]))*L)
        quadratic=model.a*model.ec*u*u*v;linear=model.b*model.ec*u*v
        activation=quadratic+linear;deactivation=model.c*u
        np.testing.assert_allclose(activation-deactivation,model.reaction(0,q*L,np.vstack([u,v]),np.zeros_like(u))[0],rtol=1e-10,atol=1e-12)
        core_active=integral(u,width,L);half_active=integral(u,L/2,L)
        core_rate=integral(activation,width,L);half_rate=integral(activation,L/2,L)
        total_act=integral(activation,L,L);total_deact=integral(deactivation,L,L)
        steady_error=abs(total_act-total_deact)/total_deact
        assert steady_error<1e-5, 'Checkpoint is not balanced at whole-cell level'
        rows.append(dict(length=L,pole_core_width=width,pole_core_active_material=core_active,
          activation_core=core_rate,activation_core_quadratic=integral(quadratic,width,L),
          activation_core_linear=integral(linear,width,L),deactivation_core=integral(deactivation,width,L),
          activation_per_active_core=core_rate/core_active,
          incumbent_half_active_material=half_active,activation_incumbent_half=half_rate,
          activation_per_active_incumbent_half=half_rate/half_active,
          activation_whole_cell=total_act,deactivation_whole_cell=total_deact,
          whole_cell_relative_balance_error=steady_error,checkpoint=c['checkpoint'],checkpoint_sha256=c['checkpoint_sha256']))
    g.n.csvwrite(DATA/'turnover_vs_length.csv',rows)
    plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    length=[r['length'] for r in rows]
    fig,ax=plt.subplots(figsize=(8.5,5.5))
    ax.plot(length,[r['activation_core'] for r in rows],'o-',color='#0072B2',label='Pole core: half-height region',ms=5,lw=2)
    ax.plot(length,[r['activation_incumbent_half'] for r in rows],'s--',color='#555555',label='Incumbent half: includes pole tails',ms=4,lw=1.7)
    ax.set(xlabel='Cell length L',ylabel='Integrated activation rate\n(concentration × length / native time)',
           title='Goryachev: steady incumbent-pole turnover',xlim=(0,40),ylim=(0,None))
    ax.legend(frameon=False);ax.grid(alpha=.2)
    fig.text(.12,.02,'Exact pre-cue basal states. Gross activation includes both native activation terms; added cue = 0.',fontsize=9,color='#4d5c66')
    fig.tight_layout(rect=[0,.045,1,1])
    for ext in ['png','pdf','svg']:fig.savefig(FIG/('turnover_vs_length.'+ext),dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.5,5.5))
    ax.plot(length,[r['activation_per_active_core'] for r in rows],'o-',color='#0072B2',label='Pole core',ms=5)
    ax.plot(length,[r['activation_per_active_incumbent_half'] for r in rows],'s--',color='#555555',label='Incumbent half',ms=4)
    ax.axhline(model.c,color='#888888',ls=':',label='Native deactivation coefficient c = 0.01733')
    ax.set(xlabel='Cell length L',ylabel='Integrated activation / active material\n(1 / native time)',title='Companion: fractional turnover',xlim=(0,40),ylim=(0,None))
    ax.legend(frameon=False);ax.grid(alpha=.2);fig.tight_layout()
    for ext in ['png','pdf','svg']:fig.savefig(FIG/('fractional_turnover.'+ext),dpi=180)
    plt.close(fig)
    stamp=g.n.now()
    g.n.write(DATA/'provenance.json',dict(generated_at=stamp,producer='finished',cases=len(rows),
          excluded_lengths=[c['length'] for c in manifest['cases'] if c['model']=='goryachev' and c['condition']=='fixed_relative' and not c['eligible']],
          activation_definition='Integral of a*ec*u^2*v + b*ec*u*v, with no added cue; positive activation rather than net reaction.',
          pole_core_definition='Left-boundary-connected region above minimum active concentration plus half the active contrast; interpolated threshold edge.',
          incumbent_half_definition='x from 0 to L/2; includes pole tails. Partial edge cells integrated as native finite-volume averages.',
          maximum_whole_cell_balance_error=max(r['whole_cell_relative_balance_error'] for r in rows),
          model_source_sha256=g.n.sha(g.CODE/'core/models/goryachev.py'),renderer_sha256=g.n.sha(Path(__file__)),
          csv_sha256=g.n.sha(DATA/'turnover_vs_length.csv'),checkpoint_hashes_pass=True))
    (FIG/'index.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>Goryachev incumbent turnover</title></head><body style="font:17px system-ui;max-width:1000px;margin:30px auto"><h1>Steady incumbent-pole activation turnover</h1><p>Generated '+stamp+'; producer finished. Native saved pre-cue basal states. Activation includes both positive native terms, excluding deactivation and additional cue.</p><img style="width:100%" src="turnover_vs_length.png"><p>Blue: half-height pole core. Dashed gray: entire incumbent half, including tails. Finite-volume integrals use fractional edge-cell overlap.</p><p><a href="turnover_vs_length.pdf">Main PDF</a> · <a href="../../data/runs/goryachev_incumbent_turnover_20261008/turnover_vs_length.csv">Data</a></p><img style="width:100%" src="fractional_turnover.png"><p>The companion divides activation rate by active material in the same region. At steady state, local activation need not equal deactivation inside a truncated pole core because diffusion crosses its edge.</p></body></html>')
    print(json.dumps(dict(cases=len(rows),maximum_balance_error=max(r['whole_cell_relative_balance_error'] for r in rows),selected=[r for r in rows if r['length'] in [9,39]]),indent=2))


if __name__=='__main__':main()
