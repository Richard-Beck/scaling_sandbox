"""Trial comparison of small incumbent suppression: morphology and restoration."""
from pathlib import Path
import argparse, csv, json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CODE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(CODE))
import three_test_reversal as t
from models.goryachev import Goryachev
from models.wave_pinning import WavePinning
n=t.n;ROOT=n.ROOT
DATA=n.DATA/'runs/goryachev_wave_pinning_trial_20261008'
FIG=ROOT/'figures/Goryachev_wave_pinning_trial'
MODELS={'goryachev':('Goryachev',Goryachev(),'#8c398f'),
        'wave_pinning':('Wave-pinning',WavePinning(),'#168577')}

def read(p):return json.loads(p.read_text())
def width(u,x):
    threshold=u.min()+.5*np.ptp(u)
    j=int(np.flatnonzero(u<threshold)[0]);assert j>0
    return float(x[j-1]+(threshold-u[j-1])/(u[j]-u[j-1])*(x[j]-x[j-1]))

def prepare():
    if DATA.exists():raise ValueError('Evidence already frozen; use render')
    DATA.mkdir();arrays={};sources={};rows=[];morph=[];timestamps={}
    for name,(label,model,color) in MODELS.items():
        run=n.DATA/f'runs/{name}_restoring_response_20261008'
        status=read(run/'results/status.json');assert status['complete'] and status['audit_passed']
        assert read(run/'results/independent_audit.json')['passed']
        timestamps[name]=status['generated_at']
        for file in ['manifest.json','results/response_vs_length.csv','results/independent_audit.json']:
            path=run/file;sources[str(path.relative_to(ROOT))]=n.sha(path)
        summary=list(csv.DictReader((run/'results/response_vs_length.csv').open()))
        manifest=read(run/'manifest.json')
        for job in manifest['jobs']:
            if job['duration']!=1.:continue
            L=job['case']['length'];directory=run/'workers'/f'{job["slot"]:03d}'
            paths=[directory/'epsilon_0.005.npz',directory/'epsilon_0.npz']
            for path in paths:
                assert n.sha(path)==read(path.with_suffix('.json'))['archive_sha256']
                sources[str(path.relative_to(ROOT))]=n.sha(path)
            cp=ROOT/job['case']['checkpoint'];assert n.sha(cp)==job['case']['checkpoint_sha256']
            sources[str(cp.relative_to(ROOT))]=n.sha(cp)
            with np.load(paths[0]) as z,np.load(paths[1]) as sham,np.load(cp) as checkpoint:
                np.testing.assert_array_equal(z['originalstate'],checkpoint['state'])
                x=z['x'];w=np.clip(job['width']-np.arange(len(x))*L/len(x),0.,L/len(x))
                u=z['recovery_active'];v=z['recovery_inactive'];a=sham['recovery_active'];b=sham['recovery_inactive']
                damage=float((a[0]-u[0])@w)
                rate=float((model.reaction(0.,x,np.stack([u[0],v[0]]),np.zeros(len(x)))[0]-
                            model.reaction(0.,x,np.stack([a[0],b[0]]),np.zeros(len(x)))[0])@w)
                row=next(r for r in summary if float(r['length'])==L and float(r['duration'])==1. and float(r['epsilon'])==.005)
                fractional=damage/float(a[0]@w);gain=rate/damage
                np.testing.assert_allclose([fractional,gain],[float(row['fractional_damage']),float(row['reaction_gain_per_damage'])],rtol=1e-8)
                rows.append(dict(model=name,length=L,fractional_loss=fractional,reaction_per_damage=gain))
                times=z['recovery_time'];height=100*(a.max(axis=1)-u.max(axis=1))/a.max(axis=1)
                widths=100*np.array([(width(aa,x)-width(uu,x))/width(aa,x) for aa,uu in zip(a,u)])
                for target in [0.,1.,10.,100.]:
                    i=int(np.argmin(abs(times-target)))
                    morph.append(dict(model=name,length=L,time=float(times[i]),peak_height_loss_percent=float(height[i]),
                                      relative_half_height_width_loss_percent=float(widths[i])))
                if L in [9.,21.,39.]:arrays[f'{name}_L{int(L)}_steady']=z['originalstate']
                if L==21.:
                    prefix=name+'_';i=int(np.argmin(abs(times-10.)))
                    arrays.update({prefix+'x':x,prefix+'time':times,prefix+'height_loss':height,prefix+'width_loss':widths,
                        prefix+'deficit_release':(a[0]-u[0])/damage,prefix+'deficit_later':(a[i]-u[i])/damage,
                        prefix+'later_time':np.array(times[i]),prefix+'core_width':np.array(job['width'])})
    np.savez_compressed(DATA/'evidence.npz',**arrays)
    n.csvwrite(DATA/'comparison.csv',rows);n.csvwrite(DATA/'morphology.csv',morph)
    n.write(DATA/'provenance.json',dict(generated_at=n.now(),producer='finished',sources=sources,source_generation_times=timestamps,
        inputs={f:n.sha(DATA/f) for f in ['evidence.npz','comparison.csv','morphology.csv']},
        native_models={name:n.sha(CODE/f'core/models/{name}.py') for name in MODELS},
        audit='Archive/checkpoint identity; independent fixed-core loss and native net-reaction integrals agree with audited summaries.',
        protocol='0.5% suppression of both basal and feedback activation inside the fixed pre-cue half-height core for one native time unit; native recovery; matched sham subtraction.',
        width_definition='First interpolated crossing of each profile\'s own min + half its contrast. This relative width separates height changes from edge movement; it differs from fixed-threshold sink diagnostics.'))
    render()

def render():
    provenance=read(DATA/'provenance.json')
    for name,digest in provenance['inputs'].items():assert n.sha(DATA/name)==digest
    for name,digest in provenance['native_models'].items():assert n.sha(CODE/f'core/models/{name}.py')==digest
    for name,digest in provenance['sources'].items():assert n.sha(ROOT/name)==digest
    FIG.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axes=plt.subplots(4,2,figsize=(14,16))
    length_colors=['#3b4cc0','#6f9fae','#a51f35'];release='#3268ac';later='#d17b28'
    with np.load(DATA/'evidence.npz') as z:
        for col,(name,(label,model,color)) in enumerate(MODELS.items()):
            ax=axes[0,col]
            for L,c in zip([9,21,39],length_colors):
                u=z[f'{name}_L{L}_steady'][0];x=(np.arange(len(u))+.5)*L/len(u)
                ax.plot(x,u,color=c,lw=2,label=f'L = {L}')
            ax.set(title=f'{"AB"[col]}  {label}: incumbent steady profiles',xlabel='Position from incumbent end, x',ylabel='Active concentration',xlim=(0,40),ylim=(0,None))
            ax.legend(frameon=False,fontsize=9)
            if col==0:
                inset=ax.inset_axes([.36,.32,.59,.49])
                for L,c in zip([9,21,39],length_colors):
                    u=z[f'{name}_L{L}_steady'][0];x=(np.arange(len(u))+.5)*L/len(u);inset.plot(x,u,color=c,lw=1.7)
                inset.set_xlim(0,2.5);inset.set_ylim(0,125);inset.tick_params(labelsize=8);inset.set_title('Pole close-up',fontsize=9)
            ax=axes[1,col];x=z[name+'_x'];steady=z[name+'_L21_steady'][0]
            backdrop=ax.twinx();backdrop.fill_between(x,0,steady/steady.max(),color='#888888',alpha=.12)
            backdrop.set_ylim(0,1.12);backdrop.set_yticks([0,1]);backdrop.set_ylabel('Steady active / peak',color='#777777',fontsize=9);backdrop.tick_params(colors='#777777',labelsize=9)
            backdrop.spines['top'].set_visible(False)
            ax.plot(x,z[name+'_deficit_release'],color=release,lw=2,label='Immediately after release')
            ax.plot(x,z[name+'_deficit_later'],color=later,lw=2,label=f'{float(z[name+"_later_time"]):.2f} after release')
            ax.axvline(float(z[name+'_core_width']),color='#999999',lw=.8,ls=':')
            ax.axhline(0,color='#777777',lw=.6)
            ax.set(title=f'{"CD"[col]}  Where active material is lost (L = 21)',xlabel='Position from incumbent end, x',
                   ylabel='Active deficit / core loss at release\n(1 / native length)',xlim=(0,3 if col==0 else 14))
            ax.legend(frameon=False,fontsize=9)
            ax=axes[2,col];time=z[name+'_time']
            ax.plot(time,z[name+'_height_loss'],color=release,lw=2,label='Peak-height loss')
            ax.plot(time,z[name+'_width_loss'],color=later,lw=2,ls='--',label='Relative half-height width loss')
            ax.axhline(0,color='#777777',lw=.8)
            ax.set(title=f'{"EF"[col]}  Height and width evolve differently (L = 21)',xlabel='Time after native activation restored',
                   ylabel='Loss relative to matched sham (%)',xscale='symlog',xlim=(0,100))
            ax.set_xscale('symlog',linthresh=.1);ax.set_xticks([0,.1,1,10,100]);ax.set_xticklabels(['0','0.1','1','10','100'])
            ax.legend(frameon=False,fontsize=9)
    rows=list(csv.DictReader((DATA/'comparison.csv').open()))
    for name,(label,model,color) in MODELS.items():
        r=sorted([r for r in rows if r['model']==name],key=lambda r:float(r['length']))
        axes[3,0].plot([float(s['length']) for s in r],[100*float(s['fractional_loss']) for s in r],'o-',color=color,label=label,ms=4)
        axes[3,1].plot([float(s['length']) for s in r],[float(s['reaction_per_damage']) for s in r],'o-',color=color,label=label,ms=4)
    axes[3,0].set(title='G  Fractional material loss at release',xlabel='Cell length L',ylabel='Active core material lost (%)',yscale='log')
    axes[3,1].set(title='H  Immediate restoring reaction per unit damage',xlabel='Cell length L',ylabel='Excess net reaction / core material lost\n(1 / native time)')
    for ax in axes[3]:ax.legend(frameon=False)
    for ax in axes.flat:ax.grid(alpha=.15)
    fig.suptitle('Trial contrast: how Goryachev and Wave-pinning incumbents deform and restore',fontsize=17,y=.99)
    fig.tight_layout(rect=[.015,.095,.985,.97],h_pad=2.)
    caption=('Matched 0.5% activation suppression for one native time unit; original half-height core held fixed. C–F use L=21; gray profiles in C–D show the steady incumbent, dotted lines its initial core edge. Deficits and losses subtract matched sham controls. Width follows each profile’s own half-height above its minimum; negative width loss means broadening. Wave-pinning initially loses height, then its residual deficit shifts toward the edge and narrowing. Goryachev remains lower and slightly broader. G–H use all 18 lengths, 9–39; net reaction excludes diffusion. Native times and lengths are not calibrated across models. No change is visually amplified.')
    import textwrap
    fig.text(.07,.074,'\n'.join(textwrap.wrap(caption,150)),va='top',fontsize=9.5,color='#465360')
    for suffix in ['png','pdf','svg']:fig.savefig(FIG/f'goryachev_wave_pinning_trial.{suffix}',dpi=180)
    plt.close(fig);(FIG/'CAPTION.txt').write_text(caption+'\n')
    n.write(DATA/'render_provenance.json',dict(generated_at=n.now(),producer='finished',renderer_sha256=n.sha(Path(__file__)),
        input_hashes_pass=True,source_hashes_pass=True,figures={str(f.name):n.sha(f) for f in FIG.iterdir() if f.is_file()}))
    print(json.dumps(dict(figure=str(FIG/'goryachev_wave_pinning_trial.png'))))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','render']);args=p.parse_args()
    prepare() if args.action=='prepare' else render()
