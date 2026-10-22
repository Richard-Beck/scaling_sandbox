"""Render audited beta=0 results with independent axes for spatial profiles.

The integration campaign and its frozen executable sources remain immutable.
"""
from pathlib import Path
import argparse,csv,shutil,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'code'),str(Path(__file__).parent)]
from debelly_beta0_reversal import n,r,BASE,OUT,TESTS,sources
import length_dev_output_layout as layout
OUT=OUT/'beta0_reversal'
def render(out,rows=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    assert r.read(out/'results/validation.json')['passed']
    m,groups,_,_=r.collect(out)
    if rows is None:rows,_=r.rows_and_records(m,groups,True)
    import csv
    native=list(csv.DictReader((BASE/'results/endpoints.csv').open()))
    fig,axs=plt.subplots(2,3,figsize=(12.7,8.4),layout='constrained',sharex=False)
    colors={'orientation':'#246a9a','full_polarity':'#b14f36'}
    for col,(geometry,title) in enumerate(zip(TESTS,['Fixed patch · 3 µm','Relative patch · 20%','Constant slope · 0.24'])):
        ax=axs[0,col];ax.set_title(title);ax.set_ylabel('Shortest passing cue (s)');ax.set_xlabel('Cell length (µm)');ax.grid(alpha=.15)
        for gate,color in colors.items():
            rs=sorted([q for q in rows if q['condition']==geometry and q['definition']==gate],key=lambda q:q['length'])
            ax.plot([q['length'] for q in rs],[q['endpoint'] if q['status']=='success' else np.nan for q in rs],color=color,lw=1.6,ls='-' if gate=='orientation' else '--')
            for q in rs:
                if q['status']=='success':ax.plot(q['length'],q['endpoint'],'o' if gate=='orientation' else 's',color=color,mfc=color if q['selected_tail_field_pass_0_1_percent'] else 'white',ms=5)
                elif q['status']=='censored':ax.plot(q['length'],q['maximum_tested_duration'],'^',color=color,mfc='white')
                elif gate=='orientation':ax.plot(q['length'],.035,'x',color='#777',transform=ax.get_xaxis_transform(),clip_on=False)
            ns=sorted([q for q in native if q['model']=='debelly' and q['condition']==geometry and q['definition']==gate and q['status']=='success'],key=lambda q:float(q['length']))
            ax.plot([float(q['length']) for q in ns],[float(q['endpoint']) for q in ns],color=color,alpha=.35,lw=1,ls=':')
        ax.set_yscale('log');ax.set_ylim(70,2600);ax.set_xlim(0,40);ax.set_xticks([0,10,20,30,40])
    for ax,(field,label) in zip(axs[1],[('initial_orientation','Basal Rac polarity magnitude'),('initial_contrast','Basal Rac contrast'),('rac','Verified beta=0 Rac profiles')]):
        ax.set_title(label);ax.grid(alpha=.15);ax.set_xlabel('Cell length (µm)' if field!='rac' else 'Position x/L')
        bs=[c for c in m['cases'] if c['condition']=='fixed_relative' and c['checkpoint']]
        if field!='rac':
            ax.set_xlim(17,40);ax.set_xticks([18,23,31,39]);vals=[abs(c[field]) for c in bs];ax.plot([c['length'] for c in bs],vals,'o-',color='#333',ms=4)
            if field=='initial_orientation':ax.axhline(.05,ls=':',color='#777',label='Original cutoff 0.05');ax.legend(frameon=False,fontsize=8)
        else:
            ax.set_xlim(0,1);ax.set_xticks([0,.25,.5,.75,1]);selected=[bs[0],bs[len(bs)//2],bs[-1]]
            for c in selected:
                _,tr=n.read_checkpoint(n.ROOT/c['checkpoint'],'debelly');ax.plot(tr.x/c['length'],tr.fields['rac'][-1],label=f'{c["length"]:g} µm')
            ax.legend(frameon=False,fontsize=8)
    handles=[Line2D([],[],color=c,marker='o' if g=='orientation' else 's',label='β=0 orientation (>0.005)' if g=='orientation' else 'β=0 own-start 1% strength') for g,c in colors.items()]
    handles+=[Line2D([],[],color='#777',ls=':',label='Native feedback (original gates)'),Line2D([],[],color='#777',marker='o',mfc='white',ls='',label='Field stability flag fails'),Line2D([],[],color='#777',marker='^',mfc='white',ls='',label='Censored'),Line2D([],[],color='#777',marker='x',ls='',label='Original preparation ineligible')]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False,fontsize=9)
    fig.suptitle('Reversal with tension→Rho feedback disabled',fontsize=15)
    artifacts={}
    for ext in ['png','pdf','svg']:
        path=layout.dev_output_path(OUT,'beta0_reversal.'+ext);fig.savefig(path,dpi=200);artifacts[layout.dev_artifact_key(OUT,path)]=n.sha(path)
    plt.close(fig)
    caption='Tension-dependent Rho activation is zero throughout new steady-state preparation, finite countercue and recovery. The original homogeneous reference and protrusion threshold, direct tension inhibition and MCA→Rac feedback remain. Each length shares its exact new checkpoint across all three cues. Top: shortest tested successful pulse, with all sampled fail→pass boundaries refined to ≤1%; recovery of 10,000 s is excluded. Blue: orientation; red: 1% replacement of the own-start polarity magnitude. Pale dotted lines show the native manuscript baseline. Open symbols fail the separate field-tail flag. Bottom: new basal asymmetry and selected full Rac profiles. Eligibility policy: '+m['eligibility_policy']+'. Asymmetric policy requires left moment <−0.005 and contrast ≥0.05; reversal requires right moment >0.005, contrast ≥0.05 and final-500-s moment range <0.005, followed by the separate 1% strength requirement. Original 0.05 eligibility and reversal flags are exported. Ten native preparation exclusions are retained; no de novo establishment search was performed.'
    cp=layout.dev_output_path(OUT,'beta0_reversal_caption.txt');cp.write_text(caption+'\n');artifacts[layout.dev_artifact_key(OUT,cp)]=n.sha(cp)
    for name in ['endpoints.csv','independent_audit.json','validation.json']:
        dest=layout.dev_output_path(OUT,'beta0_reversal_'+name);shutil.copy2(out/'results'/name,dest);artifacts[layout.dev_artifact_key(OUT,dest)]=n.sha(dest)
    r.write(layout.dev_output_path(OUT,'beta0_reversal.manifest.json'),dict(generated_at=n.now(),producer='finished',artifact_root='../..',caption=caption,campaign=str(out.relative_to(n.ROOT)),sha256=artifacts,checkpoint_sha256={str(c['length']):c['checkpoint_sha256'] for c in m['cases'] if c['condition']=='fixed_relative' and c['eligible']},sources=sources(),renderer_sha256=n.sha(Path(__file__)),reproduction=f'{sys.executable} {__file__} --out {out}'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=n.DATA/'runs/debelly_beta0_reversal_20261010');a=ap.parse_args();render(a.out.resolve())
