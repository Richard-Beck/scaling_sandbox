"""Render both mechanism figures using audited inactive removal/return trials."""
from pathlib import Path
import argparse, csv, html, json, sys, textwrap
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

CODE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE))
import three_test_reversal as t
from models.goryachev import Goryachev
from models.wave_pinning import WavePinning
n=t.n;ROOT=n.ROOT
RUN=n.DATA/'runs/inactive_pool_restoring_response_20261008'
DATA=n.DATA/'runs/inactive_pool_figure_exports_20261008'
MAIN=ROOT/'figures/Goryachev_length_mechanism/inactive_pool_perturbation'
COMPARE=ROOT/'figures/Goryachev_wave_pinning_trial/inactive_pool_perturbation'
MODELS={'goryachev':('Goryachev',Goryachev(),'#8c398f'),'wave_pinning':('Wave-pinning',WavePinning(),'#168577')}
LENGTHS=[9,15,21,27,33,39]
PALETTE=['#3b4cc0','#5977d5','#6f9fae','#cf9960','#cc5a47','#a51f35']
LC=dict(zip(LENGTHS,PALETTE));EC=['#0072B2','#D55E00','#009E73']
def read(p):return json.loads(p.read_text())
def rows(p):return list(csv.DictReader(p.open()))
def tag(D):return f'D{D:g}'.replace('.','p')
def width(u,x):
    h=u.min()+.5*np.ptp(u);j=int(np.flatnonzero(u<h)[0]);assert j>0
    return float(x[j-1]+(h-u[j-1])/(u[j]-u[j-1])*(x[j]-x[j-1]))

def prepare():
    if DATA.exists():raise ValueError('Compact evidence already frozen; use render')
    status=read(RUN/'results/status.json');assert status['complete'] and status['audit_passed']
    assert read(RUN/'results/independent_audit.json')['passed']
    manifest=read(RUN/'manifest.json');DATA.mkdir();sources={};morph=[]
    for filename in ['manifest.json','provenance.json','results/status.json','results/independent_audit.json','results/response_vs_length.csv']:
        sources[str((RUN/filename).relative_to(ROOT))]=n.sha(RUN/filename)
    for D in manifest['durations']:
        arrays={}
        for job in manifest['jobs']:
            if job['duration']!=D:continue
            name=job['case']['model'];L=job['case']['length'];model=MODELS[name][1]
            directory=RUN/'workers'/f'{job["slot"]:03d}';path=directory/'epsilon_0.005.npz';sham=directory/'epsilon_0.npz'
            for file in [path,sham]:
                assert n.sha(file)==read(file.with_suffix('.json'))['archive_sha256']
                sources[str(file.relative_to(ROOT))]=n.sha(file)
            with np.load(path) as z,np.load(sham) as s:
                arrays[f'{name}_L{int(L)}_steady']=z['originalstate']
                if L not in [9.,21.,39.]:continue
                x=z['x'];loss=s['recovery_core_active'][0]-z['recovery_core_active'][0]
                u=z['recovery_active'];v=z['recovery_inactive'];a=s['recovery_active'];b=s['recovery_inactive']
                native=model.reaction(0.,x,np.stack([u[0],v[0]]),np.zeros(len(x)))[0]
                reference=model.reaction(0.,x,np.stack([a[0],b[0]]),np.zeros(len(x)))[0]
                w=np.clip(job['width']-np.arange(len(x))*L/len(x),0.,L/len(x))
                np.testing.assert_allclose((native-reference)@w,z['recovery_core_reaction'][0]-s['recovery_core_reaction'][0],rtol=1e-8,atol=1e-12)
                time=z['recovery_time'];i=int(np.argmin(abs(time-10.)))
                height=100*(a.max(axis=1)-u.max(axis=1))/a.max(axis=1)
                widths=100*np.array([(width(aa,x)-width(uu,x))/width(aa,x) for aa,uu in zip(a,u)])
                prefix=f'{name}_L{int(L)}_'
                arrays.update({prefix+'x':x,prefix+'core_width':np.array(job['width']),prefix+'time':time,
                    prefix+'height_loss':height,prefix+'width_loss':widths,prefix+'reaction_density':(native-reference)/loss,
                    prefix+'deficit_release':(a[0]-u[0])/loss,prefix+'deficit_later':(a[i]-u[i])/loss,prefix+'later_time':np.array(time[i]),
                    prefix+'recovery_deficit':(s['recovery_core_active']-z['recovery_core_active'])/loss,
                    prefix+'recovery_reaction':(z['recovery_core_reaction']-s['recovery_core_reaction'])/loss,
                    prefix+'joined_time':np.r_[z['pulse_time']-D,time[1:]],
                    prefix+'joined_deficit':np.r_[(s['pulse_core_active']-z['pulse_core_active'])/loss,(s['recovery_core_active'][1:]-z['recovery_core_active'][1:])/loss]})
                for target in [0.,1.,10.,100.]:
                    k=int(np.argmin(abs(time-target)));morph.append(dict(model=name,length=L,duration=D,time=float(time[k]),peak_height_loss_percent=float(height[k]),relative_half_height_width_loss_percent=float(widths[k])))
        np.savez_compressed(DATA/f'evidence_{tag(D)}.npz',**arrays)
    n.csvwrite(DATA/'morphology.csv',morph)
    (DATA/'response_vs_length.csv').write_bytes((RUN/'results/response_vs_length.csv').read_bytes())
    n.write(DATA/'provenance.json',dict(generated_at=n.now(),producer='finished',durations=manifest['durations'],primary_duration=min(manifest['durations']),
        source_generated_at=status['generated_at'],sources=sources,
        inputs={file.name:n.sha(file) for file in DATA.iterdir() if file.is_file()},
        models={name:n.sha(CODE/f'core/models/{name}.py') for name in MODELS},definition=manifest['definition'],
        audit='All 288 original archives audited; compact spatial reaction densities integrate to audited fixed-core reaction. Matched sham subtraction throughout.',
        width_definition='First interpolated crossing of each profile\'s own minimum plus half its contrast. Negative width loss is broadening.'))
    render()

def load():
    p=read(DATA/'provenance.json')
    for name,digest in p['inputs'].items():assert n.sha(DATA/name)==digest
    for name,digest in p['sources'].items():assert n.sha(ROOT/name)==digest
    for name,digest in p['models'].items():assert n.sha(CODE/f'core/models/{name}.py')==digest
    numeric=rows(DATA/'response_vs_length.csv')
    for r in numeric:
        for key in r:
            if key not in ['model','archive']:r[key]=float(r[key]) if r[key] else np.nan
    return p,numeric

def save(fig,folder,stem,caption,pdf):
    fig.text(.065,.073,'\n'.join(textwrap.wrap(caption,150)),fontsize=9.2,va='top',color='#42525d')
    for suffix in ['png','pdf','svg']:fig.savefig(folder/f'{stem}.{suffix}',dpi=180)
    pdf.savefig(fig);plt.close(fig);(folder/f'{stem}_caption.txt').write_text(caption+'\n')

def main_figure(D,z,data,stem,pdf):
    fig,axes=plt.subplots(2,3,figsize=(16,10.5));inset=axes[0,0].inset_axes([.37,.38,.58,.48])
    for L,color in zip(LENGTHS,PALETTE):
        u,v=z[f'goryachev_L{L}_steady'];x=(np.arange(len(u))+.5)*L/len(u)
        for ax,field in [(axes[0,0],u),(axes[0,1],v),(inset,u)]:ax.plot(x,field,color=color,lw=1.9,label=f'L = {L}')
        axes[0,1].plot(x[-1],v[-1],'o',color=color,ms=3)
    axes[0,0].set(title='A  A taller, localized incumbent',ylabel='Steady active concentration')
    axes[0,1].set(title='B  Lower inactive concentration',ylabel='Steady inactive concentration');axes[0,1].legend(frameon=False,ncol=2,fontsize=9)
    for ax in axes[0,:2]:ax.set(xlim=(0,40),ylim=(0,None),xlabel='Position from incumbent end, x')
    inset.set(xlim=(0,2.5),ylim=(0,125),xticks=[0,1,2],yticks=[0,60,120]);inset.tick_params(labelsize=8);inset.set_title('Incumbent close-up',fontsize=9);inset.set_xlabel('Position x',fontsize=8)
    for L in [9,21,39]:
        prefix=f'goryachev_L{L}_';x=z[prefix+'x'];y=z[prefix+'reaction_density'];color=LC[L]
        axes[0,2].plot(x,y,color=color,lw=2,label=f'L = {L}');axes[0,2].fill_between(x,0,y,where=x<=float(z[prefix+'core_width']),color=color,alpha=.1)
        axes[1,0].plot(z[prefix+'joined_time'],z[prefix+'joined_deficit'],color=color,lw=2)
    axes[0,2].set(title='C  Net reaction when inactive is returned',xlabel='Position from incumbent end, x',ylabel='Excess net reaction / active material lost\n(1 / length / native time)',xlim=(0,2.5));axes[0,2].legend(frameon=False,fontsize=9);axes[0,2].axhline(0,color='#777777',lw=.7)
    axes[1,0].axvspan(-D,0,color='#aab2bb',alpha=.17);axes[1,0].axvline(0,color='#777777',ls=':',lw=.8)
    axes[1,0].set(title='D  Brief pool reduction, then return',xlabel='Time relative to returning inactive material',ylabel='Active core deficit / deficit at return',xlim=(-D,2),ylim=(-.02,1.08))
    axes[1,0].text(.9,1.015,'Removed inactive returned at t = 0',fontsize=9,ha='center',color='#53616b')
    for epsilon,color in zip([.001,.005,.01],EC):
        r=sorted([r for r in data if r['model']=='goryachev' and r['duration']==D and r['epsilon']==epsilon],key=lambda r:r['length'])
        axes[1,1].plot([s['length'] for s in r],[100*s['fractional_damage'] for s in r],'o-',color=color,label=f'{epsilon*100:g}% pool reduction',ms=4)
        axes[1,2].plot([s['length'] for s in r],[s['reaction_gain_per_damage'] for s in r],'o-',color=color,ms=4)
    axes[1,1].set(title='E  Fractional pole loss during pool reduction',xlabel='Cell length L',ylabel='Active core material lost (%)');axes[1,1].legend(frameon=False,fontsize=9)
    axes[1,2].set(title='F  Restoring reaction per unit pole loss',xlabel='Cell length L',ylabel='Immediate excess net reaction / material lost\n(1 / native time)')
    for ax in axes.flat:ax.grid(alpha=.15)
    fig.suptitle(f'Goryachev: incumbent response to a brief inactive-pool perturbation ({D:g} time units)',fontsize=16,y=.985)
    fig.tight_layout(rect=[.015,.12,.995,.955])
    caption=(f'Remove 0.1%, 0.5% or 1% of the inactive profile throughout the cell, evolve for {D:g} native time units, then add back exactly the removed profile. Native kinetics and diffusion never change. A–B: unscaled steady concentrations. C–D: 0.5% removal; matched sham subtraction. C shows net reaction immediately after return; shading marks the fixed original half-height core. D normalizes core deficit to its own loss at return. E–F: all three amplitudes at 18 lengths. Net reaction excludes diffusion; the brief pulse tests a specific restoring response.')
    save(fig,MAIN,stem,caption,pdf)

def comparison_figure(D,z,data,stem,pdf):
    fig,axes=plt.subplots(4,2,figsize=(14,16));release='#3268ac';later='#d17b28'
    for col,(name,(label,model,color)) in enumerate(MODELS.items()):
        ax=axes[0,col]
        for L in [9,21,39]:
            u=z[f'{name}_L{L}_steady'][0];x=(np.arange(len(u))+.5)*L/len(u);ax.plot(x,u,color=LC[L],lw=2,label=f'L = {L}')
        ax.set(title=f'{"AB"[col]}  {label}: incumbent steady profiles',xlabel='Position from incumbent end, x',ylabel='Active concentration',xlim=(0,40),ylim=(0,None));ax.legend(frameon=False,fontsize=9)
        if col==0:
            inset=ax.inset_axes([.36,.32,.59,.49])
            for L in [9,21,39]:
                u=z[f'{name}_L{L}_steady'][0];x=(np.arange(len(u))+.5)*L/len(u);inset.plot(x,u,color=LC[L],lw=1.7)
            inset.set(xlim=(0,2.5),ylim=(0,125));inset.tick_params(labelsize=8);inset.set_title('Pole close-up',fontsize=9)
        prefix=f'{name}_L21_';x=z[prefix+'x'];steady=z[name+'_L21_steady'][0];ax=axes[1,col]
        backdrop=ax.twinx();backdrop.fill_between(x,0,steady/steady.max(),color='#888888',alpha=.12);backdrop.set(ylim=(0,1.12),yticks=[0,1]);backdrop.set_ylabel('Steady active / peak',color='#777777',fontsize=9);backdrop.tick_params(colors='#777777',labelsize=9)
        ax.plot(x,z[prefix+'deficit_release'],color=release,lw=2,label='Immediately after inactive return')
        ax.plot(x,z[prefix+'deficit_later'],color=later,lw=2,label=f'{float(z[prefix+"later_time"]):.2f} after return')
        ax.axvline(float(z[prefix+'core_width']),color='#999999',ls=':',lw=.8);ax.axhline(0,color='#777777',lw=.6)
        ax.set(title=f'{"CD"[col]}  Where active material is lost (L = 21)',xlabel='Position from incumbent end, x',ylabel='Active deficit / core loss at return\n(1 / native length)',xlim=(0,3 if col==0 else 14));ax.legend(frameon=False,fontsize=8.5)
        ax=axes[2,col];time=z[prefix+'time'];ax.plot(time,z[prefix+'height_loss'],color=release,lw=2,label='Peak-height loss');ax.plot(time,z[prefix+'width_loss'],color=later,lw=2,ls='--',label='Relative half-height width loss')
        ax.axhline(0,color='#777777',lw=.8);ax.set(title=f'{"EF"[col]}  Height and width after pool return (L = 21)',xlabel='Time after inactive material returned',ylabel='Loss relative to matched sham (%)',xlim=(0,100));ax.set_xscale('symlog',linthresh=.1);ax.set_xticks([0,.1,1,10,100]);ax.set_xticklabels(['0','0.1','1','10','100']);ax.legend(frameon=False,fontsize=9)
    styles=['--','-',':']
    for name,(label,model,color) in MODELS.items():
        for eps,style in zip([.001,.005,.01],styles):
            r=sorted([r for r in data if r['model']==name and r['duration']==D and r['epsilon']==eps],key=lambda r:r['length'])
            for ax,key,factor in [(axes[3,0],'fractional_damage',100),(axes[3,1],'reaction_gain_per_damage',1)]:
                ax.plot([s['length'] for s in r],[factor*s[key] for s in r],color=color,ls=style,lw=1.8,marker='o' if eps==.005 else None,ms=3)
    handles=[Line2D([],[],color=color,label=label,lw=2) for label,model,color in MODELS.values()]+[Line2D([],[],color='#444444',ls=style,label=f'{eps*100:g}% pool reduction') for eps,style in zip([.001,.005,.01],styles)]
    axes[3,0].set(title='G  Fractional material loss at return',xlabel='Cell length L',ylabel='Active core material lost (%)',yscale='log')
    axes[3,1].set(title='H  Immediate restoring reaction per unit damage',xlabel='Cell length L',ylabel='Excess net reaction / active core material lost\n(1 / native time)')
    axes[3,1].legend(handles=handles,frameon=False,fontsize=8.5,ncol=2,loc='upper left')
    for ax in axes.flat:ax.grid(alpha=.15)
    fig.suptitle(f'Inactive-pool perturbation: Goryachev versus Wave-pinning ({D:g} time units)',fontsize=17,y=.99);fig.tight_layout(rect=[.015,.10,.985,.97],h_pad=2.)
    caption=(f'Remove a fraction of the entire inactive profile for {D:g} native time units, then return that exact material pointwise. Native reactions and diffusion remain unchanged. C–F: 0.5% removal at L=21; gray profiles show the steady incumbent, dotted lines its initial core edge. All deficits subtract matched shams. Width uses each profile’s own half-height above its minimum; negative loss means broadening or a higher peak. G–H: 0.1%, 0.5%, 1% removal across all 18 lengths. G uses a log scale; H curves nearly overlap within each model. Net reaction excludes diffusion. Native lengths and times are not calibrated across models.')
    save(fig,COMPARE,stem,caption,pdf)

def recovery_figure(D,z,data,pdf):
    fig,axes=plt.subplots(2,2,figsize=(13,9))
    for col,(name,(label,model,color)) in enumerate(MODELS.items()):
        for L in [9,21,39]:
            prefix=f'{name}_L{L}_';time=z[prefix+'time'];axes[0,col].plot(time,z[prefix+'recovery_deficit'],color=LC[L],label=f'L = {L}');axes[1,col].plot(time,z[prefix+'recovery_reaction'],color=LC[L])
        axes[0,col].set_title(label);axes[0,col].legend(frameon=False)
        axes[0,col].set_ylabel('Remaining core deficit / deficit at return');axes[1,col].set_ylabel('Excess net reaction / deficit at return\n(1 / native time)')
    for ax in axes.flat:ax.set_xscale('symlog',linthresh=.001);ax.set_xlim(0,1000);ax.set_xlabel('Time after inactive material returned');ax.axhline(0,color='#777777',lw=.6);ax.grid(alpha=.15)
    fig.suptitle(f'Recovery after 0.5% inactive removal for {D:g} native time units',fontsize=16);fig.tight_layout(rect=[.02,.10,.98,.95])
    save(fig,MAIN,f'recovery_{tag(D)}','Matched sham subtraction; same fixed original half-height core throughout. Inactive material is returned at time zero; native reactions and diffusion remain unchanged. Deficits are normalized to each case’s active loss at return. These traces belong to the inactive-pool protocol, not the earlier activation-suppression experiment.',pdf)

def render():
    p,data=load();MAIN.mkdir(parents=True,exist_ok=True);COMPARE.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    primary=p['primary_duration'];mainpdf=PdfPages(MAIN/'inactive_pool_report.pdf');comparepdf=PdfPages(COMPARE/'inactive_pool_comparison_report.pdf')
    mainstems=[];comparestems=[]
    for D in p['durations']:
        mainstem='Goryachev_length_mechanism' if D==primary else f'Goryachev_length_mechanism_{tag(D)}'
        comparestem='goryachev_wave_pinning_trial' if D==primary else f'goryachev_wave_pinning_trial_{tag(D)}'
        with np.load(DATA/f'evidence_{tag(D)}.npz') as z:
            main_figure(D,z,data,mainstem,mainpdf);recovery_figure(D,z,data,mainpdf);comparison_figure(D,z,data,comparestem,comparepdf)
        mainstems.append(mainstem);comparestems.append(comparestem)
    mainpdf.close();comparepdf.close();stamp=n.now()
    protocol=p['definition']
    note='The figures use the brief '+f'{primary:g}'+'-unit removal as primary and the longer pulse as a timing sensitivity check. Activation kinetics are never suppressed. Instantaneous zero-duration removal and return would produce no damage. Width is a relative half-height measure; negative width loss indicates broadening. The material returned is the original removed profile, added to the evolved inactive field.'
    for folder,title,stems,pdfname in [(MAIN,'Goryachev length mechanism: inactive-pool perturbation',mainstems,'inactive_pool_report.pdf'),(COMPARE,'Goryachev and Wave-pinning: inactive-pool perturbation',comparestems,'inactive_pool_comparison_report.pdf')]:
        (folder/'index.html').write_text('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><style>body{font:17px/1.5 system-ui;max-width:1550px;margin:30px auto;padding:0 22px;color:#253440}img{width:100%}.stamp{font-size:14px;color:#52616a}</style></head><body><h1>'+html.escape(title)+'</h1><p>'+html.escape(protocol)+'</p><p>'+html.escape(note)+'</p><p><a href="'+pdfname+'">PDF report</a> · <a href="../../../data/runs/inactive_pool_figure_exports_20261008/README.md">Evidence and methods</a></p><p class="stamp">Figures generated '+stamp+'; producer finished. Source data generated '+p['source_generated_at']+'; all 288 archive audits pass.</p>'+''.join('<section><h2>'+('Primary perturbation' if i==0 else 'Longer-pulse sensitivity check')+'</h2><img src="'+stem+'.png" alt="'+html.escape(title)+'"><p>'+html.escape((folder/f'{stem}_caption.txt').read_text())+'</p></section>' for i,stem in enumerate(stems))+'</body></html>')
    n.write(DATA/'render_provenance.json',dict(generated_at=stamp,producer='finished',renderer_sha256=n.sha(Path(__file__)),input_hashes_pass=True,source_hashes_pass=True,
        figures={str(file.relative_to(ROOT)):n.sha(file) for folder in [MAIN,COMPARE] for file in folder.iterdir() if file.is_file()}))
    print(json.dumps(dict(main_figure=str(MAIN/'Goryachev_length_mechanism.png'),comparison=str(COMPARE/'goryachev_wave_pinning_trial.png'))))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','render']);args=p.parse_args()
    prepare() if args.action=='prepare' else render()
