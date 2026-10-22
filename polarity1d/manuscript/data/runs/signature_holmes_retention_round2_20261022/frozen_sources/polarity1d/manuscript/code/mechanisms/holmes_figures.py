"""Report assets for the source-backed Holmes fixed-slope mechanism study."""
from __future__ import annotations
import csv,gzip,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import holmes_mechanism as h

BLUE='#2563eb';ORANGE='#ea7b25';PURPLE='#7c3aed';GREEN='#059669';GRAY='#8b929a';PINK='#db2777'
L_COLORS={10:'#008b8b',28:'#dc7c20',50:'#7c3aed'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,
                    'axes.spines.top':False,'axes.spines.right':False,'legend.fontsize':8,
                    'pdf.fonttype':42,'svg.fonttype':'none'})
ASSETS=[]

def load(L,gate='orientation',result='pass'):
    tag=f'L{L}_{gate}_{result}'
    if L==10 and gate=='full_polarity':tag=f'L10_orientation_{result}'
    return np.load(h.OUT/(tag+'.npz')),json.loads((h.OUT/(tag+'.json')).read_text())

def save(fig,name,claim,caption):
    h.FIG.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf']:
        fig.savefig(h.FIG/(name+'.'+ext),dpi=220,facecolor='white')
    ASSETS.append(dict(name=name,pdf=str(h.FIG/(name+'.pdf')),png=str(h.FIG/(name+'.png')),
                       claim=claim,caption=caption))
    plt.close(fig)

def footer(fig,text):
    fig.text(.035,.025,text,fontsize=9,va='bottom',linespacing=1.35)

def snapshots():
    fig,axs=plt.subplots(3,6,figsize=(11,7),gridspec_kw={'left':.065,'right':.985,'top':.84,'bottom':.175,'wspace':.20,'hspace':.52})
    fig.suptitle('Holmes 3 | Fast right seeding does not remove the old Rac domain',x=.035,y=.98,ha='left',fontsize=16,weight='bold')
    fig.text(.035,.925,'Exact native states; fixed slope g = 0.0667. Each row uses its shortest tested ordinary-orientation pulse.',fontsize=10)
    fig.text(.035,.891,'L = 28 is an ordinary-duration minimum, tied with L = 30; the strength-duration minimum is L = 24.',fontsize=9,color='#505964')
    for row,L in enumerate([10,28,50]):
        z,m=load(L);x=z['x']/L;T=m['duration'];p=z['pulse_state'];pt=z['pulse_time'];r=z['recovery_state'];rt=z['recovery_time']
        frames=[p[0],p[np.abs(pt-T/2).argmin()],p[-1],r[np.abs(rt-10).argmin()],r[np.abs(rt-100).argmin()],r[-1]]
        titles=['Initial',f'On: t = {T/2:.3g}',f'Cutoff: t = {T:.4g}','Recovery +10','Recovery +100','Recovery +10,000']
        for col,(s,title) in enumerate(zip(frames,titles)):
            ax=axs[row,col];ax.plot(x,s[0],color=BLUE,lw=1.8);ax.plot(x,s[2],color=ORANGE,lw=1.7)
            ax.set_xlim(0,1);ax.set_ylim(0,{10:2.1,28:2.65,50:3.55}[L]);ax.set_xticks([0,.5,1]);ax.set_title(title,pad=7)
            if col in [1,2]:ax.set_facecolor('#f5f2ff')
            if col==0:ax.set_ylabel(f'L = {L}\nConcentration')
            if row==2:ax.set_xlabel('Position x/L')
            if col==5:
                ax.text(.98,.97,f'P = {z["recovery_orientation"][-1]:.4f}',transform=ax.transAxes,ha='right',va='top',fontsize=8)
                if L>10:
                    ax.annotate('old-domain\nremnant',xy=({28:.23,50:.33}[L],1.48),xytext=(.02,.67),textcoords='axes fraction',
                                fontsize=7.5,arrowprops=dict(arrowstyle='-',color='#555555',lw=.7))
        if row==0:axs[row,0].legend([Line2D([],[],color=BLUE),Line2D([],[],color=ORANGE)],['Active Rac','Active Rho'],loc='lower left',frameon=False,fontsize=7.5)
    footer(fig,'Rac and Rho suppress each other: the cue raises Rac at the right, lowers Rho there and creates a second Rac-rich domain.\n'
               'At L = 28 and 50 the old domain retreats from the left edge into the interior, yet survives the 10,000-unit recovery.\n'
               'The weak rightward first moment passes even though two Rac-high domains remain. Times are native seconds; cue cutoff is exact.')
    save(fig,'01_spatial_snapshots','Longer cells seed the new right domain quickly but retain a displaced old Rac domain.',
         'Spatial Rac/Rho profiles at three actual sampled lengths. Purple panel backgrounds identify the finite cue. '
         'Preparation is the exact relative_20pct checkpoint and recovery uses native basal input. The final panels are measured endpoints, not independently verified new steady states. '
         'Orientation-selected durations are 8.0625, 3.14375 and 7.25. Native concentration units are micromolar.')

def domains(y,x,threshold):
    mask=y>threshold;d=np.diff(np.r_[False,mask,False].astype(int));a=np.flatnonzero(d==1);b=np.flatnonzero(d==-1)
    return [(float(x[i]),float(x[j-1])) for i,j in zip(a,b)]

def evidence():
    fig,axs=plt.subplots(2,3,figsize=(11,7),gridspec_kw={'left':.065,'right':.985,'top':.835,'bottom':.28,'wspace':.29,'hspace':.60})
    fig.suptitle('Holmes 3 | Pool communication sustains the long-length rise',x=.035,y=.98,ha='left',fontsize=16,weight='bold')
    fig.text(.035,.926,'Stronger stimulation with L accelerates right seeding; speeding both inactive pools removes the rising arm in this diagnostic.',fontsize=10)
    rows=[r for r in csv.DictReader((h.CAMPAIGN/'results/endpoints.csv').open()) if r['model']=='holmes_model3' and r['condition']=='relative_20pct']
    ax=axs[0,0]
    for gate,label,color in [('orientation','Ordinary orientation',BLUE),('full_polarity','1% scalar strength',PURPLE)]:
        rr=sorted([r for r in rows if r['definition']==gate and r['status']=='success'],key=lambda r:float(r['length']))
        xx=np.array([float(r['length']) for r in rr]);yy=np.array([float(r['endpoint']) for r in rr])
        ax.plot(xx,yy,color=color,lw=1.5,label=label)
        for r in rr:
            ax.plot(float(r['length']),float(r['endpoint']),'o',ms=3.5,mfc=color if r['selected_field_tail_pass']=='True' else 'white',mec=color)
    ctrl=[json.loads((h.OUT/'mixing_control'/f'bracket_L{L}.json').read_text()) for L in [10,28,50]]
    xx=[r['length'] for r in ctrl];yy=[r['shortest_tested_full_strength_duration'] for r in ctrl]
    ax.plot(xx,yy,'s--',color=GREEN,label='Fast pools: strength dose',ms=4)
    ax.set_title('A  Native curve and transport control',loc='left');ax.set_xlabel('Length');ax.set_ylabel('Finite cue duration');ax.set_ylim(0,11.8);ax.set_xlim(8,52)
    ax.legend(frameon=False,fontsize=7.2,loc='upper center');ax.axvline(28,color='#dddddd',lw=.8,zorder=0)
    ax=axs[0,1];seed=[];threshold=[];slopepeaks=[]
    for L in [10,28,50]:
        z,m=load(L);pt=z['pulse_time'];ps=z['pulse_state'];oldpeak=z['preparation_state'][0].max()
        seed.append(float(pt[np.flatnonzero(ps[:,0,-1]>=oldpeak)[0]]));threshold.append(m['duration']);slopepeaks.append(m['slope']*L)
    k=np.arange(3);ax.bar(k-.18,seed,.34,color=GREEN,label='Tip reaches old Rac plateau');ax.bar(k+.18,threshold,.34,color=BLUE,label='Required ordinary duration')
    for i,v in enumerate(seed):ax.text(i-.18,v+.12,f'{v:.2f}',ha='center',fontsize=7.5)
    ax.set_xticks(k,[f'L = {L}' for L in [10,28,50]]);ax.set_ylim(0,10.5);ax.set_ylabel('Time');ax.set_title('B  Tip activation gets faster, not slower',loc='left')
    ax.legend(frameon=False,fontsize=7,loc='upper center')
    ax=axs[0,2]
    for L in [10,28,50]:
        z,m=load(L);x=z['x']/L;ps=z['pulse_state'][-1]
        ax.plot(x,ps[1]-ps[1].mean(),color=L_COLORS[L],label=f'Rac inactive, L={L}')
        ax.plot(x,ps[3]-ps[3].mean(),color=L_COLORS[L],ls='--',alpha=.7)
    ax.axhline(0,color='#cccccc',lw=.6);ax.set_xlim(0,1);ax.set_ylim(-1.05,.8);ax.set_xticks([0,.5,1]);ax.set_xlabel('Position x/L');ax.set_ylabel('Deviation from spatial mean')
    ax.set_title('C  Pools develop gradients during cue',loc='left');ax.legend(frameon=False,fontsize=6.8,loc='lower left')
    ax.text(.98,.97,'Solid: inactive Rac\nDashed: inactive Rho',ha='right',va='top',transform=ax.transAxes,fontsize=7.2)
    for col,L in enumerate([10,28,50]):
        ax=axs[1,col];z,m=load(L);q,mm=load(L,'full_polarity');x=z['x']/L
        fast=json.loads((h.OUT/'mixing_control'/f'bracket_L{L}.json').read_text())['shortest_tested_full_strength_duration']
        f=np.load(h.OUT/'mixing_control'/f'L{L}_T{fast:g}.npz')
        ax.plot(x,z['recovery_state'][-1,0],color=BLUE,lw=1.8,label=f'Ordinary: T={m["duration"]:.4g}')
        ax.plot(x,q['recovery_state'][-1,0],color=PURPLE,ls='--',lw=1.8,label=f'Strength: T={mm["duration"]:.4g}')
        ax.plot(f['x']/L,f['final_state'][0],color=GREEN,ls=':',lw=2,label=f'Fast pools: T={fast:.4g}')
        mesh=h.OUT/'mesh_check'/f'L{L}_orientation_N401.npz'
        if mesh.exists():
            mz=np.load(mesh);ax.plot(mz['x'][::12]/L,mz['final_state'][0,::12],'.',color='#444444',ms=2.2,label='401-cell ordinary diagnostic')
        ax.set_ylim(.27,1.63);ax.set_xlim(0,1);ax.set_xticks([0,.5,1]);ax.set_xlabel('Position x/L');ax.set_title(f'{"DEF"[col]}  Final active Rac at L = {L}',loc='left')
        if col==0:ax.set_ylabel('Active Rac')
        ax.legend(frameon=False,fontsize=6.8,loc='upper center',bbox_to_anchor=(.5,-.26),ncol=1)
    footer(fig,'The slope is fixed, so ramp peaks are 0.667, 1.867 and 3.333 at L = 10, 28 and 50. Direct input also stimulates the incumbent domain.\n'
               'Two inactive pools diffuse at 50; membrane-active Rac/Rho diffuse at 0.1. Faster communication (inactive diffusion 5,000) yields\n'
               'one right domain at shorter doses. This changes transport physics; zero controls stay left, with maximum state excursion 0.82%.')
    save(fig,'02_pool_mechanism','Faster inactive-pool coupling removes the long-length upturn, while right-tip activation already accelerates with length.',
         'A shows the native measured curve and a explicitly perturbed transport diagnostic. Green durations are locally refined shortest tested full-strength successes '
         '(7.8125, 2.4375, 1.859375), with brackets below 1%; its coarse scan was 0,1,2,3,4,6,8,10,12. Open native markers preserve failed field-tail checks. '
         'B uses reaching the initial active-Rac maximum at the right boundary as a seeding benchmark, not a commitment clock. C shows mean-centered inactive profiles '
         'at each native ordinary cutoff. D-F compare final shapes at ordinary, stronger and faster-transport doses. All start the same original state; the faster-transport '
         'zero controls quantify a modest shift in the state, so this intervention supports transport-dependent coupling but does not prove a universal diffusion-time formula. '
         'Gray points show a 401-cell diagnostic with conservative-group-preserving interpolation of the original fields; the two-domain ordinary endpoints remain at L28 and50. '
         'Finer-grid zero controls remain left, with maximum frozen-scale excursion0.13%; these bounded examples are not a systematic quantitative mesh-convergence study.')

def recovery():
    fig,axs=plt.subplots(2,3,figsize=(11,7),gridspec_kw={'left':.065,'right':.985,'top':.835,'bottom':.26,'wspace':.29,'hspace':.58})
    fig.suptitle('Holmes 3 | Pulse thresholds select different post-cue fates',x=.035,y=.98,ha='left',fontsize=16,weight='bold')
    fig.text(.035,.927,'Neighboring failed and passing trials restart the same checkpoint. The gate is evaluated after 10,000 native units of recovery.',fontsize=10)
    for col,L in enumerate([10,28,50]):
        ax=axs[0,col];strength=abs(load(L)[1]['initial_orientation'])
        ax.axhspan(.99*strength,1.01*strength,color=PURPLE,alpha=.12);ax.axhline(.05,color=BLUE,ls=':',lw=1)
        for gate,result,color,style in [('orientation','pass',BLUE,'-'),('orientation','fail',GRAY,'-'),('full_polarity','pass',PURPLE,'-'),('full_polarity','fail',PINK,'-')]:
            if L==10 and gate=='full_polarity':continue
            z,m=load(L,gate,result);rt=z['recovery_time'];o=z['recovery_orientation'];ax.plot(rt,o,color=color,lw=1.4,label=f'{gate.replace("full_polarity","Strength").replace("orientation","Ordinary")} {result}: {m["duration"]:.4g}')
            tag=f'L{L}_{gate}_{result}';path=h.OUT/'continuation'/(tag+'.npz')
            if path.exists():
                c=np.load(path);ax.plot(c['time'],c['orientation'],color=color,ls='--',lw=1)
        ax.axvline(10000,color='#999999',ls=':',lw=.8);ax.set_xscale('symlog',linthresh=1);ax.set_xlim(0,100000);ax.set_ylim(-.36,.38)
        ax.set_xticks([0,10,1000,10000,100000],['0','10','1k','10k','100k']);ax.set_xlabel('Time since cue withdrawal');ax.set_title(f'{"ABC"[col]}  L = {L}: scalar recovery',loc='left')
        ax.get_xticklabels()[-1].set_horizontalalignment('right')
        if col==0:ax.set_ylabel('Orientation P')
        ax.legend(loc='upper center',frameon=False,fontsize=6.2)
        ax.text(.01,.045,'dashed = added native continuation',transform=ax.transAxes,fontsize=6.5,color='#505964')
        ax=axs[1,col]
        for result,color in [('pass',PURPLE),('fail',PINK)]:
            z,m=load(L,'full_polarity',result);ax.plot(z['x']/L,z['recovery_state'][-1,0],color=color,lw=1.8,label=f'Strength {result}: T={m["duration"]:.4g}')
        ax.set_xlim(0,1);ax.set_ylim(.25,1.65);ax.set_xticks([0,.5,1]);ax.set_xlabel('Position x/L');ax.set_title(f'{"DEF"[col]}  Adjacent strength doses, same gradient',loc='left',fontsize=9)
        if col==0:ax.set_ylabel('Active Rac at recovery 10k')
        ax.legend(frameon=False,fontsize=7.2,loc='upper center',bbox_to_anchor=(.5,-.23))
    footer(fig,'L = 10: a small dose difference chooses which end survives after competition. L = 28: slow coarsening makes the finite horizon matter.\n'
               'The failed strength dose T = 3.64148 fully reverses by recovery 100k; its ordinary lower neighbor T = 3.125 also later passes.\n'
               'L = 50: the ordinary two-domain state remains essentially unchanged to 100k; a longer dose is needed to remove the interior remnant.')
    save(fig,'03_withdrawal_fates','The finite recovery gate can accept durable two-domain patterns and can reject slow eventual reversals.',
         'Solid histories extend through the retained 10,000-unit recovery; dashed lines add an exact native no-cue continuation to recovery 100,000. '
         'The blue dotted line is the ordinary minimum first moment and the purple band is the 1%-strength target. The retained gates additionally require contrast '
         'and a stable final-500 scalar window. At L10 the same pass/fail pair also brackets strength. The lower panels show native full-strength pass/fail endpoints '
         'at the assay horizon. Added continuations diagnose time-horizon dependence and do not alter any manuscript duration or original success label.')

def audit_and_evidence():
    records=json.load(gzip.open(h.CAMPAIGN/'results/records.json.gz','rt'))['cases'];checks=[];summary=[]
    for task in h.trials():
        case=next(c for c in records if c['case_index']==task['case']['case_index']);orig=next(t for t in case['tests'] if t['duration']==task['duration'])
        z=np.load(h.OUT/(task['tag']+'.npz'));q=np.load(h.CAMPAIGN/orig['trace']);pf=orig['trace_prefix'];m=json.loads((h.OUT/(task['tag']+'.json')).read_text())
        row=dict(tag=task['tag'],orientation_gate_matches=m['orientation_success']==orig['orientation_success'],
                 strength_gate_matches=m['full_polarity_success']==orig['full_polarity_success'],field_flag_matches=m['field_tail_pass']==orig['field_tail_pass'],
                 max_finalstate_error=float(np.max(abs(q[pf+'finalstate']-z['recovery_state'][-1]))),
                 orientation_endpoint_error=abs(orig['final_orientation']-float(z['recovery_orientation'][-1])))
        checks.append(row)
    h.write(h.OUT/'reconstruction_audit.json',dict(generated_at=h.n.now(),checks=checks,maximum_accepted_absolute_field_error=1e-7,
          passed=all(r['orientation_gate_matches'] and r['strength_gate_matches'] and r['field_flag_matches'] and r['max_finalstate_error']<1e-7 for r in checks)))
    for L in [10,28,50]:
        z,m=load(L);s=z['pulse_state'];t=z['pulse_time'];prep=z['preparation_state'];thr=(prep[0].max()+prep[0].min())/2
        seed=float(t[np.flatnonzero(s[:,0,-1]>=prep[0].max())[0]])
        q,mm=load(L,'full_polarity')
        summary.append(dict(length=L,ordinary_duration=m['duration'],strength_duration=mm['duration'],initial_orientation=m['initial_orientation'],
              ordinary_final_orientation=float(z['recovery_orientation'][-1]),strength_final_orientation=float(q['recovery_orientation'][-1]),
              ordinary_domain_intervals_halfheight=domains(z['recovery_state'][-1,0],z['x'],thr),
              strength_domain_intervals_halfheight=domains(q['recovery_state'][-1,0],q['x'],thr),
              tip_reaches_initial_Rac_high=seed,Rac_inactive_cutoff_range=float(np.ptp(s[-1,1])),
              Rho_inactive_cutoff_range=float(np.ptp(s[-1,3])),inactive_first_mode_mixing_scale=L*L/(np.pi*np.pi*50.)))
    h.write(h.OUT/'mechanism_evidence.json',dict(generated_at=h.n.now(),summary=summary,
            native_source_hashes=h.n.sources(),script_hashes={str(p):h.n.sha(p) for p in [Path(__file__),Path(h.__file__)]},
            literature_source='https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002366'))

def main():
    audit_and_evidence();snapshots();evidence();recovery();h.write(h.FIG/'ordered_assets.json',ASSETS)
    print('\n'.join(a['pdf'] for a in ASSETS))

if __name__=='__main__':main()
