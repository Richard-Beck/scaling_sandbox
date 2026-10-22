"""Canonical Goryachev length mechanism: weaker challenge, stronger incumbent.

prepare freezes compact audited evidence; render needs only those inputs and
the original native checkpoints, allowing exploratory presentations to retire.
"""
from pathlib import Path
import argparse, csv, html, json, shutil, sys
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

CODE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(CODE));import three_test_reversal as t
n=t.n
from models.goryachev import Goryachev
ROOT=n.ROOT
DATA=n.DATA/'runs/goryachev_length_mechanism_20261008'
FIG=ROOT/'figures/Goryachev_length_mechanism'
CAMPAIGN=n.DATA/'runs/three_test_reversal_20261006'
RESTORE=n.DATA/'runs/goryachev_restoring_response_20261008'
ILLUSTRATION=n.DATA/'runs/goryachev_length_illustration_20261008'
COLORS=['#0072B2','#D55E00','#009E73']


def read(path):return json.loads(path.read_text())
def csvread(path):return list(csv.DictReader(path.open()))
def csvwrite(path,rows):n.csvwrite(path,rows)
def now():return datetime.now(timezone.utc).isoformat()
def weights(width,L,N):return np.clip(width-np.arange(N)*L/N,0.,L/N)


def illustrate():
    """Freeze spatial evidence for the revised figure; original arrays stay intact."""
    if ILLUSTRATION.exists():raise ValueError('Illustration evidence already frozen; use render')
    p,states,response=load();ILLUSTRATION.mkdir()
    sources={};arrays={};model=Goryachev()
    for s in states:
        L=s['length'];cp=ROOT/s['checkpoint'];sources[str(cp.relative_to(ROOT))]=n.sha(cp)
        with np.load(cp) as z:arrays[f'state_L{int(L)}']=z['state'].copy()
    manifest=read(RESTORE/'manifest.json')
    for L in [9.,21.,39.]:
        job=next(j for j in manifest['jobs'] if j['case']['length']==L and j['duration']==1.)
        directory=RESTORE/'workers'/f'{job["slot"]:03d}'
        paths=[directory/'epsilon_0.005.npz',directory/'epsilon_0.npz']
        for path in paths:
            assert n.sha(path)==read(path.with_suffix('.json'))['archive_sha256']
            sources[str(path.relative_to(ROOT))]=n.sha(path)
        with np.load(paths[0]) as z,np.load(paths[1]) as sham:
            loss=sham['recovery_core_active'][0]-z['recovery_core_active'][0]
            x=z['x'];w=weights(job['width'],L,len(x))
            net=model.reaction(0.,x,np.stack([z['recovery_active'][0],z['recovery_inactive'][0]]),np.zeros(len(x)))[0]
            reference=model.reaction(0.,x,np.stack([sham['recovery_active'][0],sham['recovery_inactive'][0]]),np.zeros(len(x)))[0]
            excess=net-reference
            np.testing.assert_allclose(excess@w,z['recovery_core_reaction'][0]-sham['recovery_core_reaction'][0],rtol=1e-8,atol=1e-12)
            prefix=f'L{int(L)}_'
            arrays.update({prefix+'x':x,prefix+'width':np.array(job['width']),prefix+'spatial_restoring':excess/loss,
                prefix+'time':np.r_[z['pulse_time']-1.,z['recovery_time'][1:]],
                prefix+'deficit':np.r_[(sham['pulse_core_active']-z['pulse_core_active'])/loss,
                                     (sham['recovery_core_active'][1:]-z['recovery_core_active'][1:])/loss]})
    path=ILLUSTRATION/'spatial_evidence.npz';np.savez_compressed(path,**arrays)
    n.write(ILLUSTRATION/'provenance.json',dict(created_at=now(),producer='finished',sources=sources,
        archive_sha256=n.sha(path),native_model_sha256=n.sha(CODE/'core/models/goryachev.py'),
        definition='Exact steady checkpoints; 0.5% one-unit activation suppression; sham-subtracted native net reaction immediately after release, normalized by fixed-core active material lost. Spatial integral independently agrees with audited core reaction.'))
    render()


def prepare():
    if (DATA/'inputs').exists():raise ValueError('Compact inputs already frozen; use render')
    for run in [CAMPAIGN,RESTORE]:assert read(run/'results/status.json')['complete']
    assert read(RESTORE/'results/independent_audit.json')['passed']
    inp=DATA/'inputs';inp.mkdir(parents=True)
    cases=sorted([c for c in read(CAMPAIGN/'manifest.json')['cases'] if c['model']=='goryachev' and c['condition']=='fixed_relative' and c['eligible'] and c['length']>=9.],key=lambda c:c['length'])
    model=Goryachev();states=[]
    for c in cases:
        cp=ROOT/c['checkpoint'];assert n.sha(cp)==c['checkpoint_sha256']
        with np.load(cp) as z:u,v=z['state'].copy()
        L=c['length'];q=(np.arange(len(u))+.5)/len(u);level=u.min()+.5*np.ptp(u)
        j=int(np.flatnonzero(u<level)[0]);width=float((q[j-1]+(level-u[j-1])/(u[j]-u[j-1])*(q[j]-q[j-1]))*L)
        activation=model.a*model.ec*u*u*v+model.b*model.ec*u*v
        wc=weights(width,L,len(u));wh=weights(L/2,L,len(u));balance=abs(activation.sum()-model.c*u.sum())/(model.c*u.sum())
        assert balance<1e-5
        signal=float(t.response_cue(c)(0.,q*L,L)[-1])
        states.append(dict(length=L,inactive_right=float(v[-1]),active_peak=float(u.max()),pole_width=width,
             right_cue_at_last_cell=signal,right_initial_added_activation=signal*float(v[-1]),
             core_active=float(u@wc),incumbent_active=float(u@wh),core_turnover=float(activation@wc),
             incumbent_turnover=float(activation@wh),fractional_incumbent_turnover=float(activation@wh/(u@wh)),
             whole_cell_balance_error=float(balance),checkpoint=c['checkpoint'],checkpoint_sha256=c['checkpoint_sha256']))
    csvwrite(inp/'starting_states.csv',states)
    responses=[r for r in csvread(CAMPAIGN/'results/endpoints.csv') if r['model']=='goryachev' and float(r['length'])>=9.]
    csvwrite(inp/'native_response.csv',responses)
    shutil.copy2(RESTORE/'results/response_vs_length.csv',inp/'restoring_response.csv')
    snap=n.DATA/'runs/goryachev_native_snapshots_20261008'
    if not snap.exists():snap=n.DATA/'archive/goryachev_exploration_20261008/goryachev_native_snapshots_20261008'
    assert n.sha(snap/'snapshots.csv')==read(snap/'provenance.json')['csv_sha256']
    shutil.copy2(snap/'snapshots.csv',inp/'native_snapshots.csv')
    manifest=read(RESTORE/'manifest.json');curves=[]
    for L in [9.,21.,39.]:
        for D in [1.,10.]:
            job=next(j for j in manifest['jobs'] if j['case']['length']==L and j['duration']==D)
            directory=RESTORE/'workers'/f'{job["slot"]:03d}'
            p=directory/'epsilon_0.005.npz';sham=directory/'epsilon_0.npz'
            assert n.sha(p)==read(p.with_suffix('.json'))['archive_sha256']
            assert n.sha(sham)==read(sham.with_suffix('.json'))['archive_sha256']
            with np.load(p) as z,np.load(sham) as s,np.load(ROOT/job['case']['checkpoint']) as cp:
                np.testing.assert_array_equal(z['originalstate'],cp['state'])
                arrays=dict(time=z['recovery_time'],active=z['recovery_core_active'],sham_active=s['recovery_core_active'],
                    net_reaction=z['recovery_core_reaction'],sham_net_reaction=s['recovery_core_reaction'],
                    material_rate=z['recovery_core_total_rate'],sham_material_rate=s['recovery_core_total_rate'])
            dest=inp/f'recovery_L{int(L)}_D{int(D)}.npz';np.savez_compressed(dest,**arrays)
            curves.append(dict(length=L,duration=D,source=str(p.relative_to(ROOT)),sha256=n.sha(p),sham_source=str(sham.relative_to(ROOT)),sham_sha256=n.sha(sham)))
    n.write(DATA/'provenance.json',dict(created_at=now(),inputs={str(p.relative_to(DATA)):n.sha(p) for p in inp.iterdir()},
       checkpoints=[dict(path=s['checkpoint'],sha256=s['checkpoint_sha256']) for s in states],
       source_model_sha256=n.sha(CODE/'core/models/goryachev.py'),restoring_sources=curves,
       restoration_original_manifest_sha256=n.sha(RESTORE/'manifest.json'),
       restoration_audit_sha256=n.sha(RESTORE/'results/independent_audit.json'),
       restoration_source_generated_at=read(RESTORE/'results/status.json')['generated_at'],
       native_source_generated_at=read(CAMPAIGN/'results/status.json')['generated_at'],
       native_snapshot_provenance=read(snap/'provenance.json'),
       scope='Two supported coupled effects: weaker initial local countercue activation and stronger immediate local incumbent restoring response. Lengths 9–39; no exhaustive attribution of reversal duration.'))
    render()


def load():
    p=read(DATA/'provenance.json')
    for name,digest in p['inputs'].items():assert n.sha(DATA/name)==digest
    for cp in p['checkpoints']:assert n.sha(ROOT/cp['path'])==cp['sha256']
    assert n.sha(CODE/'core/models/goryachev.py')==p['source_model_sha256']
    numeric=lambda rows:[{k:float(v) if k in ['length','duration','epsilon','inactive_right','active_peak','pole_width','right_initial_added_activation','core_active','incumbent_active','core_turnover','incumbent_turnover','fractional_incumbent_turnover','fractional_damage','reaction_gain_per_damage','half_damage_recovery_time','ninety_percent_damage_recovery_time','active_damage_at_release','excess_net_activation_at_release'] and v else v for k,v in r.items()} for r in rows]
    return p,numeric(csvread(DATA/'inputs/starting_states.csv')),numeric(csvread(DATA/'inputs/restoring_response.csv'))


def save(fig,pdf,stem,title,caption):
    fig.suptitle(title,fontsize=16,y=.985)
    fig.tight_layout(rect=[.015,.105,.995,.955])
    import textwrap
    fig.text(.065,.078,'\n'.join(textwrap.wrap(caption,155)),fontsize=9.1,va='top',color='#3d4d58')
    fig.savefig(FIG/(stem+'.png'),dpi=180);pdf.savefig(fig)
    if stem=='Goryachev_length_mechanism':
        fig.savefig(FIG/(stem+'.pdf'));fig.savefig(FIG/(stem+'.svg'))
    plt.close(fig)
    return dict(stem=stem,title=title,caption=caption)


def render():
    p,states,response=load();FIG.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    L=[r['length'] for r in states];native=csvread(DATA/'inputs/native_response.csv')
    pages=[];pdf=PdfPages(FIG/'Goryachev_length_mechanism_report.pdf')
    illustration=read(ILLUSTRATION/'provenance.json')
    assert n.sha(ILLUSTRATION/'spatial_evidence.npz')==illustration['archive_sha256']
    assert n.sha(CODE/'core/models/goryachev.py')==illustration['native_model_sha256']
    for source,digest in illustration['sources'].items():assert n.sha(ROOT/source)==digest
    fig,axes=plt.subplots(2,3,figsize=(16,10.5))
    lengths=[9,15,21,27,33,39]
    palette=['#3b4cc0','#5977d5','#6f9fae','#cf9960','#cc5a47','#a51f35']
    length_colors=dict(zip(lengths,palette))
    inset=axes[0,0].inset_axes([.37,.38,.58,.48])
    with np.load(ILLUSTRATION/'spatial_evidence.npz') as spatial:
        for length,color in zip(lengths,palette):
            u,v=spatial[f'state_L{length}'];x=(np.arange(len(u))+.5)*length/len(u)
            for ax,field in [(axes[0,0],u),(axes[0,1],v),(inset,u)]:
                ax.plot(x,field,color=color,lw=1.9,label=f'L = {length}')
            axes[0,1].plot(x[-1],v[-1],'o',ms=3,color=color)
        for ax in axes[0,:2]:
            ax.set_xlim(0,40);ax.set_xlabel('Position from incumbent end, x');ax.set_ylim(bottom=0)
        axes[0,0].set(title='A  A taller, localized incumbent',ylabel='Steady active concentration')
        axes[0,1].set(title='B  Less inactive material per unit length',ylabel='Steady inactive concentration')
        axes[0,1].legend(frameon=False,ncol=2,loc='upper right',fontsize=9)
        inset.set_xlim(0,2.5);inset.set_ylim(0,125);inset.set_xticks([0,1,2]);inset.set_yticks([0,60,120])
        inset.tick_params(labelsize=8);inset.set_title('Incumbent close-up',fontsize=9)
        inset.set_xlabel('Position x',fontsize=8)
        for length in [9,21,39]:
            prefix=f'L{length}_';color=length_colors[length];x=spatial[prefix+'x']
            reaction=spatial[prefix+'spatial_restoring'];width=float(spatial[prefix+'width'])
            axes[0,2].plot(x,reaction,color=color,lw=2,label=f'L = {length}')
            axes[0,2].fill_between(x,0,reaction,where=x<=width,color=color,alpha=.10)
            axes[1,0].plot(spatial[prefix+'time'],spatial[prefix+'deficit'],color=color,lw=2,label=f'L = {length}')
    axes[0,2].set(title='C  Restoring reaction arises inside the pole',xlabel='Position from incumbent end, x',
                  ylabel='Excess net reaction / material lost\n(1 / length / native time)',xlim=(0,2.5))
    axes[0,2].axhline(0,color='#777777',lw=.7);axes[0,2].legend(frameon=False,fontsize=9)
    axes[1,0].axvspan(-1,0,color='#aab2bb',alpha=.16)
    axes[1,0].axvline(0,color='#777777',lw=.8,ls=':')
    axes[1,0].set(title='D  Small damage, then immediate repair',xlabel='Time relative to restoring native activation',
                  ylabel='Active core deficit / deficit at release',xlim=(-1,2),ylim=(-.02,1.07))
    axes[1,0].text(-.5,.90,'Activation\nsuppressed',ha='center',fontsize=9,color='#53616b')
    axes[1,0].text(1.0,1.015,'Native activation restored',ha='center',fontsize=9,color='#53616b')
    for eps,color in zip([.001,.005,.01],COLORS):
        rs=sorted([r for r in response if r['duration']==1. and r['epsilon']==eps],key=lambda r:r['length'])
        axes[1,1].plot([r['length'] for r in rs],[100*r['fractional_damage'] for r in rs],'o-',color=color,ms=4,label=f'{100*eps:g}% suppression')
        axes[1,2].plot([r['length'] for r in rs],[r['reaction_gain_per_damage'] for r in rs],'o-',color=color,ms=4)
    axes[1,1].set(title='E  Longer cells lose a smaller fraction',ylabel='Active core material lost (%)',xlabel='Cell length L')
    axes[1,1].legend(frameon=False,fontsize=9)
    axes[1,2].set(title='F  Stronger reaction per unit damage',ylabel='Immediate excess net activation / material lost\n(1 / native time)',xlabel='Cell length L')
    for ax in axes.flat:ax.grid(alpha=.15)
    pages.append(save(fig,pdf,'Goryachev_length_mechanism','Goryachev: a depleted inactive pool and a more resistant incumbent',
        'A–B: full steady profiles; colors identify length. Less inactive concentration gives less initial local activation under an equal countercue. C–D: suppress incumbent activation by 0.5% for one time unit, then restore native kinetics. C: immediate excess net reaction per unit lost material; shading marks the half-height core. D: core deficit relative to its own loss at release. Both subtract matched sham controls. E–F: length trends for three small suppression amplitudes; F curves overlap. Stronger immediate restoration does not imply faster late recovery.'))
    # Native field snapshots: directly from the compact, hash-checked native CSV.
    snaps=csvread(DATA/'inputs/native_snapshots.csv');fig,axes=plt.subplots(4,6,figsize=(16,10))
    labels=['Before cue','Early cue','50% of pulse','Withdrawal','50 after withdrawal','Final recovery']
    for block,length in enumerate([9.,39.]):
        for col,label in enumerate(labels):
            rs=[r for r in snaps if float(r['length'])==length and r['stage']==label]
            rs.sort(key=lambda r:float(r['x']));assert len(rs)==201
            for field,offset,color,ylim in [('active',0,COLORS[0],(0,30.) if length==9 else (0,126.)),('inactive',1,COLORS[1],(0,.62))]:
                ax=axes[2*block+offset,col];ax.plot([float(r['x']) for r in rs],[float(r[field]) for r in rs],color=color,lw=1.8)
                if rs[0]['phase']=='pulse' and label!='Before cue':ax.axvspan(.8*length,length,color='#999999',alpha=.15)
                ax.set_xlim(0,length);ax.set_ylim(*ylim);ax.grid(alpha=.15);ax.set_xticks([0,length/2,length]);ax.tick_params(labelsize=8)
                if col==0:ax.set_ylabel(f'L={length:g}\n{field.capitalize()} concentration',fontsize=10)
                if offset==1:ax.set_xlabel('Position x',fontsize=9)
                else:ax.set_title(label+'\nt='+f'{float(rs[0]["time_since_cue_onset"]):g}',fontsize=9)
    pages.append(save(fig,pdf,'native_active_inactive','Native relative-20% reversal: the chemical fields',
        'Blue: active Cdc42. Orange: inactive Cdc42. Gray: right 20% cue support while the cue is on, peak 6.4. Shortest passing tested native pulses: 70.5 at L=9 and 2680 at L=39. Withdrawal snapshots are immediately before removal; native recovery follows. Inactive axes are identical; active scales differ by length. Both retained reversals pass the scalar strength and field-tail checks. These are concentration snapshots, not normalized fields.'))
    fig,axes=plt.subplots(2,2,figsize=(13,10))
    for length,color in zip([9.,21.,39.],COLORS):
        with np.load(DATA/f'inputs/recovery_L{int(length)}_D1.npz') as z:
            loss=z['sham_active'][0]-z['active'][0];time=z['time']
            axes[0,0].plot(time,(z['sham_active']-z['active'])/loss,color=color,label=f'L={length:g}')
            axes[0,1].plot(time,(z['net_reaction']-z['sham_net_reaction'])/loss,color=color)
    for ax in axes[0]:ax.set_xscale('symlog',linthresh=.001);ax.set_xlim(0,1000);ax.set_xlabel('Time after native activation restored');ax.axhline(0,color='#777777',lw=.8)
    axes[0,0].set(title='A  Immediate repair and a slower late tail',ylabel='Remaining active deficit / deficit at release');axes[0,0].legend(frameon=False)
    axes[0,1].set(title='B  The larger initial reaction is brief',ylabel='Excess net activation / initial deficit\n(1 / native time)')
    for D,color in zip([1.,10.],COLORS):
        rs=sorted([r for r in response if r['duration']==D and r['epsilon']==.005],key=lambda r:r['length'])
        axes[1,0].plot([r['length'] for r in rs],[r['reaction_gain_per_damage'] for r in rs],'o-',color=color,label=f'{D:g}-unit pulse',ms=4)
        axes[1,1].plot([r['length'] for r in rs],[r['ninety_percent_damage_recovery_time'] for r in rs],'o-',color=color,ms=4)
    axes[1,0].set(title='C  Stronger initial response for both pulses',ylabel='Immediate excess net activation / material lost\n(1 / native time)');axes[1,0].legend(frameon=False)
    axes[1,1].set(title='D  Late recovery lengthens',ylabel='Time to sustained 90% recovery')
    for ax in axes[1]:ax.set_xlabel('Cell length L')
    for ax in axes.flat:ax.grid(alpha=.2)
    pages.append(save(fig,pdf,'restoring_dynamics','Incumbent resistance is an immediate local response, not uniformly faster recovery',
        'All panels use 0.5% activation suppression in the fixed incumbent core, with native recovery and matched sham subtraction. A–B follow a one-unit pulse. C–D compare one- and ten-unit pulses. The normalized immediate response is larger in longer cells for both durations, but its fast component decays and late recovery is slower. Recovery times are first sustained recorded tolerance crossings through 5000 units. The response is specific to this perturbation, not a global stability eigenvalue.'))
    fig,axes=plt.subplots(1,2,figsize=(13,8))
    axes[0].plot(L,[s['incumbent_turnover'] for s in states],'o-',color=COLORS[0]);axes[0].set(title='Steady activation turnover increases',xlabel='Cell length L',ylabel='Integrated activation in incumbent half\n(concentration × length / native time)')
    axes[1].plot(L,[s['fractional_incumbent_turnover'] for s in states],'o-',color=COLORS[1]);axes[1].axhline(.01733,color='#777777',ls=':');axes[1].set(title='Fractional steady turnover is nearly constant',xlabel='Cell length L',ylabel='Integrated activation / active material\n(1 / native time)',ylim=(0,.025))
    for ax in axes:ax.grid(alpha=.2)
    # Reserve extra caption space for the methods appendix.
    fig.suptitle('Turnover context and evidence checks',fontsize=16,y=.98);fig.tight_layout(rect=[.03,.32,.98,.94])
    note=('Steady turnover alone does not establish resistance: a larger pole activates and deactivates more material, while its fractional turnover stays near 0.01733. The perturbation tests establish the immediate restoring response directly.\n\n'
          'Evidence: exact shared relative-20% native checkpoints at 18 lengths (9–39); both native activation terms; no added cue during incumbent suppression; actual material conserved. All 144 perturbation archives passed independent restart, reaction/diffusion observable and conservation checks. Normalized reaction gains vary by at most 0.045% across tested amplitudes for one-unit pulses, and 0.151% for ten-unit pulses. Compact report inputs and starting checkpoints are hash checked.\n\n'
          'Limits: inactive concentration is not absolute inactive amount. The stronger incumbent response is immediate, local and perturbation-specific. The two mechanisms share conserved material and are not independent causal fractions. Diffusion and later competition still participate in reversal. Original 201-cell grids are retained; no new systematic mesh-convergence study is claimed.')
    import textwrap
    fig.text(.07,.28,'\n\n'.join('\n'.join(textwrap.wrap(para,145)) for para in note.split('\n\n')),fontsize=10,va='top')
    fig.savefig(FIG/'turnover_and_methods.png',dpi=180);pdf.savefig(fig);plt.close(fig)
    pages.append(dict(stem='turnover_and_methods',title='Turnover context and evidence checks',caption=note))
    pdf.close();stamp=now()
    prose=[('Weaker initial activation at the challenger end','Longer cells have lower inactive Cdc42 concentration at the opposite end. Added activation is the local stimulus multiplied by that concentration. The same patch peak therefore produces a smaller initial local activation rate in the longer cell. For a relative-20% patch, the stimulated region also grows with length; integrated cue input additionally depends on that width. A constant-slope gradient has a growing peak and is not an equal-local-stimulus comparison.'),
           ('A stronger immediate restoring response at the incumbent','Briefly suppressing both incumbent activation terms causes less fractional loss in longer cells. Once native kinetics are restored, the increase in net activation per unit material lost is larger. For a 0.5% one-unit suppression, L=39 loses 0.00189% of core active material versus 0.00833% at L=9; its immediate reaction response per unit damage is 4.681 versus 0.371, approximately 12.6 times larger.'),
           ('One coupled length mechanism','These observations support weaker challenger initiation together with stronger local incumbent resistance in longer cells. They arise in the same conserved active/inactive system. The restoring response is brief and the late recovery tail is slower in longer cells; the evidence does not claim uniformly faster recovery, an exhaustive explanation of cue duration, or additive causal percentages.')]
    htmlprose=''.join('<h2>'+html.escape(title)+'</h2><p>'+html.escape(body)+'</p>' for title,body in prose)
    htmlpages=''.join('<section><h2>'+html.escape(page['title'])+'</h2><img src="'+page['stem']+'.png" alt="'+html.escape(page['title'])+'"><p>'+html.escape(page['caption']).replace('\n','<br>')+'</p></section>' for page in pages)
    (FIG/'index.html').write_text('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Goryachev length mechanism</title><style>body{font:17px/1.55 system-ui;max-width:1500px;margin:32px auto;padding:0 24px;color:#20313d}img{width:100%}section{margin-top:40px;border-top:1px solid #ddd}a{color:#12628c}.stamp{color:#52616a;font-size:14px}</style></head><body><h1>Goryachev length mechanism</h1><p>Weaker countercue activation and stronger immediate incumbent resistance.</p><p><a href="Goryachev_length_mechanism.pdf">Main figure PDF</a> · <a href="Goryachev_length_mechanism_report.pdf">Four-page report PDF</a> · <a href="MECHANISM.md">Mechanism and methods</a> · <a href="../../data/runs/goryachev_length_mechanism_20261008/inputs/restoring_response.csv">Perturbation measurements</a></p><p class="stamp">Generated '+stamp+'; all source producers finished. Native campaign data: '+p['native_source_generated_at']+'. Restoring experiment data: '+p['restoration_source_generated_at']+'.</p>'+htmlprose+htmlpages+'</body></html>')
    text='# Goryachev length mechanism\n\nGenerated '+stamp+'; all source producers finished.\n\n'
    text+='\n\n'.join('## '+title+'\n\n'+body for title,body in prose)
    text+='\n\n## Measurements and interpretation limits\n\n'+note+'\n\n'
    text+='The main figure shows full physical-position active/inactive profiles at six lengths (9, 15, 21, 27, 33, 39), with a close-up of the localized active pole. The spatial restoring response and short damage/recovery trace use lengths 9, 21, 39. Their matching colors connect the starting state to the perturbation response. These panels replace the redundant local activation curve and native cue-duration outcome plot. Spatial and pulse evidence is frozen and hash checked in `data/runs/goryachev_length_illustration_20261008/`.\n\n'
    text+='The length trends use all 18 lengths from 9–39. Pole suppression and restoring measurements use the fixed pre-cue half-height core, with fractional edge-cell finite-volume integration. Net activation means gross activation minus deactivation; the plotted excess subtracts the matched sham response. Actual material recovery includes active diffusion and is stored separately.\n\n'
    text+='Rebuild with `code/mechanisms/goryachev_length_mechanism.py render`. New restoring experiments use `code/mechanisms/goryachev_restoring_response.py init --out <fresh-directory>`. Compact inputs, provenance and cleanup inventory live under `data/runs/goryachev_length_mechanism_20261008/`. Earlier exploratory presentations are archived; the native comparison and audited raw evidence remain preserved.\n'
    (FIG/'MECHANISM.md').write_text(text)
    (FIG/'CAPTION.txt').write_text(pages[0]['caption']+'\n')
    n.write(DATA/'render_provenance.json',dict(generated_at=stamp,producer='finished',input_hashes_pass=True,checkpoint_hashes_pass=True,
             illustration_provenance_sha256=n.sha(ILLUSTRATION/'provenance.json'),illustration_archive_sha256=n.sha(ILLUSTRATION/'spatial_evidence.npz'),
             renderer_sha256=n.sha(Path(__file__)),figures={str(f.relative_to(FIG)):n.sha(f) for f in FIG.iterdir() if f.is_file()},pages=len(pages)))
    print(json.dumps(dict(report=str(FIG/'Goryachev_length_mechanism_report.pdf'),main_figure=str(FIG/'Goryachev_length_mechanism.png'),pages=len(pages))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','illustrate','render']);args=p.parse_args()
    prepare() if args.action=='prepare' else illustrate() if args.action=='illustrate' else render()
