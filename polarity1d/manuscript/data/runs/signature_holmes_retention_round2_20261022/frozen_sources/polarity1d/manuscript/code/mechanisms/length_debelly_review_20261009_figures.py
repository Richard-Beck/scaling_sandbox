"""Second review: native deBelly withdrawal, transport and static mechanics figures.

Native integrations are immutable SLURM outputs. Only brief exact one-step budget
checks and static elliptic controls run here. Archived drawings supply style only.
"""
from pathlib import Path
from types import SimpleNamespace
from datetime import datetime,timezone
import json,sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import length_revision_20261009_figures as f
ROOT=f.ROOT
sys.path[:0]=[str(ROOT/'code/mechanisms')]
import debelly_corrected_reversal as corrected
import debelly_corrected_advection as ad
n=corrected.n
RUN=ROOT/'data/runs/length_debelly_review_20261009'
FLOW=ROOT/'data/runs/debelly_corrected_advection_20261008'
REV=ROOT/'data/runs/debelly_corrected_advection_revision_20261008'
ZERO=ROOT/'data/runs/debelly_corrected_zero_crossing_20261008'
f.FIG=ROOT/'figures/length_mechanisms_20261008/review_20261009_debelly'
f.DATA=ROOT/'data/runs/length_mechanisms_20261008/review_20261009_debelly'
f.FIG.mkdir(parents=True,exist_ok=True);f.DATA.mkdir(parents=True,exist_ok=True)
read,load,register,panel,csvwrite=f.read,f.load,f.register,f.panel,f.csvwrite
BLUE,ORANGE,PURPLE,GREY='#0072B2','#D55E00','#762a83','#777777'
COLORS={18:GREY,19:'#CC79A7',21:'#009E73',23:BLUE,25:'#56B4E9',27:'#009E73',29:'#8B3F82',31:'#D55E00',33:'#E69F00',35:'#222222',37:'#CC79A7',39:ORANGE}

def save(fig,stem,caption,methods):
    f.save(fig,stem,caption,methods)
    p=f.FIG/(stem+'.manifest.json');m=json.loads(p.read_text())
    m['renderer_sha256']=f.sha(Path(__file__))
    m['reproduction'][0]=m['reproduction'][0].replace('length_revision_20261009_figures.py',Path(__file__).name)
    m['shared_plot_helpers_sha256']=f.sha(Path(f.__file__))
    p.write_text(json.dumps(m,indent=2)+'\n')

def selected(path,keys):
    register(path)
    with np.load(path) as z:return {k:z[k].copy() for k in keys}
def polarity(y,x,L):return y@(2*x/L-1)/y.sum(axis=-1)
def x50(y,x,L):
    one=y.ndim==1;y=np.atleast_2d(y);a=np.cumsum(y,axis=1);a/=a[:,-1,None]
    edges=np.r_[0,x+L/len(x)/2]
    out=np.array([np.interp(.5,np.r_[0,c],edges) for c in a])
    return out[0] if one else out

def audit():
    m=read(RUN/'manifest.json');assert m['center_boundary_motion'] is False
    assert n.model_for('debelly').parameters.center_boundary_motion is False
    for path,digest in m['sources'].items():assert f.sha(ROOT/path)==digest,path
    records=[]
    for slot,j in enumerate(m['jobs']):
        p=RUN/f'workers/{slot:03d}.npz';meta=read(p.with_suffix('.json'));assert f.sha(p)==meta['archive_sha256'];register(p)
        assert meta['job']==j;c=j['case'];L=c['length'];cp=load(ROOT/c['checkpoint'])
        assert f.sha(ROOT/c['checkpoint'])==c['checkpoint_sha256']
        with np.load(p) as z:
            x=z['x'];D=meta['withdrawal_time'];idx=int(np.argmin(abs(z['pulse_time']-D)))
            fields={};max_start=0.
            for k in n.DB_NAMES:
                np.testing.assert_array_equal(z['originalstate_'+k],cp['state_'+k])
                np.testing.assert_allclose(z['pulse_field_'+k][0],cp['state_'+k],rtol=1e-12,atol=1e-12)
                np.testing.assert_allclose(z['recovery_field_'+k][0],z['pulse_field_'+k][idx],rtol=1e-12,atol=1e-12)
                for phase in ['pulse','recovery']:
                    a=z[phase+'_field_'+k];assert np.isfinite(a).all()
                    if k in ['rac','rho','mca_bound','mca_unbound']:assert a.min()>-1e-10
            for phase in ['pulse','recovery']:
                a=z[phase+'_field_rac'];o=polarity(a,x,L);q=np.ptp(a,axis=1)
                np.testing.assert_allclose(o,z[phase+'_orientation'],atol=2e-13)
                np.testing.assert_allclose(q,z[phase+'_contrast'],atol=2e-13)
            tail=z['recovery_time']>=9500;o=polarity(z['recovery_field_rac'][tail],x,L);q=np.ptp(z['recovery_field_rac'][tail],axis=1)
            error=float(np.max(abs(o/abs(c['initial_orientation'])-1)))
            orientation=bool(min(o)>.05 and min(q)>=.05 and np.ptp(o)<.005)
            maximum=0.
            for k in ['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension','physical_length']:
                a=z['recovery_field_boundary_separation'][tail]+L if k=='physical_length' else z['recovery_field_'+k][tail]
                if k=='displacement':a=a-a.mean(axis=1)[:,None]
                scale=c['fixed_scales'][k]
                maximum=max(maximum,float(np.ptp(a,axis=0).max()/scale),float(abs(a[-1]-a[-2]).max()/2.5*500/scale))
            checks=dict(orientation_success=orientation,full_polarity_success=bool(orientation and error<=.01),field_tail_pass=maximum<=.001)
            for k,v in checks.items():assert bool(meta['metrics'][k])==v,(slot,k)
            assert abs(maximum-meta['metrics']['field_tail_max_scaled'])<1e-10
            assert abs(o[-1]-meta['metrics']['final_orientation'])<1e-12
            if 'expected' in j:
                for k,v in checks.items():assert bool(j['expected'][k])==v
            records.append(dict(slot=slot,length=L,condition=c['condition'],purpose=j['purpose'],duration=D,**checks,
                                final_rac=float(o[-1]),field_tail_max_scaled=maximum,checkpoint_identity=True,phase_continuity=True))
        if slot%10==0:print('Audited',slot,flush=True)
    a=dict(passed=True,generated_at=datetime.now(timezone.utc).isoformat(),center_boundary_motion=False,cases=records,
           checks='Exact checkpoints, frozen sources and output hashes; independent Rac polarity, contrast, orientation/strength gates and full-field tail flag; finite chemistry and phase continuity.')
    (RUN/'independent_audit.json').write_text(json.dumps(a,indent=2)+'\n')
    (RUN/'status.json').write_text(json.dumps(dict(producer='finished',complete=True,jobs=len(records),execution_errors=0,audit_passed=True,generated_at=a['generated_at']),indent=2)+'\n')
    return m,records

def job_index(m,L,geo='fixed_slope',purpose='pass'):
    return next(i for i,j in enumerate(m['jobs']) if j['case']['length']==L and j['case']['condition']==geo and j['purpose']==purpose)
def trajectory(m,L,geo='fixed_slope',purpose='pass'):
    slot=job_index(m,L,geo,purpose);p=RUN/f'workers/{slot:03d}.npz'
    keys=['x']+[s+'_'+k for s in ['pulse','recovery'] for k in ['time','field_rac','field_mca_bound']]
    z=selected(p,keys);meta=read(p.with_suffix('.json'));return z,meta

def continuous(L,m):
    if L<=21:
        z,meta=trajectory(m,L,purpose='continuous_and_zero_withdrawal')
        return z['pulse_time'],z['x'],z['pulse_field_rac'],z['pulse_field_mca_bound']
    p=FLOW/f'L{L:g}_native.npz';summ=read(FLOW/f'L{L:g}_summary.json');assert f.sha(p)==summ['native_sha256']
    z=selected(p,['time','x','field_rac','field_mca_bound']);return z['time'],z['x'],z['field_rac'],z['field_mca_bound']

def measure(m,records):
    result=[]
    for L in [18,19,21,23,25,27,29,31,33,35,37,39]:
        t,x,r,b=continuous(L,m);pr=polarity(r,x,L);pm=polarity(b,x,L)
        ir=np.flatnonzero(pr>=0)[0];im=np.flatnonzero(pm<=0)[0]
        z,meta=trajectory(m,L);D=meta['withdrawal_time'];end=polarity(z['pulse_field_mca_bound'][-1],x,L)
        shift=L/2*(pm[0]-pm[im]);under=L/2*(pm[0]-end)
        zero_rec=next((a for a in records if a['length']==L and a['purpose']=='continuous_and_zero_withdrawal'),None)
        result.append(dict(length=L,rac_zero=t[ir],mca_zero=t[im],minimum_orientation_cue=D,
            initial_mca_polarity=pm[0],mca_polarity_at_cut=end,nominal_centroid_shift_to_zero=shift,
            nominal_centroid_shift_under_minimum_cue=under,mean_net_centroid_speed_to_zero=shift/t[im],
            initial_mca_x50=x50(b[0],x,L),initial_mca_x50_distance_from_midpoint=abs(x50(b[0],x,L)-L/2),
            mca_x50_distance_at_minimum_cut=abs(x50(z['pulse_field_mca_bound'][-1],x,L)-L/2),
            minimum_cue_field_tail_pass=meta['metrics']['field_tail_pass'],
            zero_withdrawal_orientation_success=zero_rec['orientation_success'] if zero_rec else (True if L in [23,39] else None)))
    csvwrite(f.DATA/'memory_measurements.csv',result);return result

def memory_figure(m,rr):
    fig=plt.figure(figsize=(11.5,11.4),layout='constrained');gs=fig.add_gridspec(3,6)
    for col,L in enumerate([18,23,39]):
        for row in [0,1]:
            ax=fig.add_subplot(gs[row,col*2:col*2+2]);panel(ax,chr(65+3*row+col),f'L = {L} µm' if row==0 else None)
            for purpose,color,ls,label in [('pass',BLUE,'-','Shortest tested reversing cue'),('fail',ORANGE,'--','Adjacent failing cue'),('zero',PURPLE,':','Withdrawal at MCA zero')]:
                if purpose=='zero':
                    if L==18:z,meta=trajectory(m,L,purpose='continuous_and_zero_withdrawal');D=meta['withdrawal_time']
                    else:
                        p=ZERO/f'withdrawal_L{L}_offset0.npz';summary=read(ZERO/'withdrawal_summary.json');entry=next(a for a in summary if a['length']==L and a['offset']==0);assert f.sha(p)==entry['archive_sha256'];z=selected(p,['x','duration']+[s+'_'+k for s in ['pulse','recovery'] for k in ['time','field_rac','field_mca_bound']]);D=float(z['duration'])
                else:z,meta=trajectory(m,L,purpose=purpose);D=meta['withdrawal_time']
                field='rac' if row==0 else 'mca_bound';sgn=1 if row==0 else -1
                keep=z['pulse_time']<=D+1e-8
                t=np.r_[z['pulse_time'][keep]-D,z['recovery_time'][1:]]
                values=np.r_[polarity(z['pulse_field_'+field][keep],z['x'],L),polarity(z['recovery_field_'+field][1:],z['x'],L)]*sgn
                ax.plot(t,values,color=color,ls=ls,lw=1.6,label=label)
                if purpose=='pass':ax.plot(t[-1],values[-1],marker='o',ms=5,color=color,mfc=color if meta['metrics']['field_tail_pass'] else 'white',clip_on=False)
            ax.axvspan(-800,0,color='#eeeeee');ax.axvline(0,color=GREY,lw=.9);ax.axhline(0,color=GREY,lw=.7)
            ax.set_xscale('symlog',linthresh=100);ax.set_xlim(-800,10000)
            ax.set_xticks([-500,0,100,1000,10000],['−500','0','100','1,000','10,000'])
            if row==0:ax.axhline(.05,color=GREY,ls=':',lw=.7);ax.set_ylim(-.16,.34)
            else:ax.set_ylim(-.17,.15)
            if col==0:ax.set_ylabel('Rac polarity toward new front' if row==0 else '−MCA polarity (new rear bias)')
            if row==1:ax.set_xlabel('Time relative to withdrawal (s)')
            if row==0 and col==1:ax.legend(frameon=False,fontsize=7,loc='upper left')
    ax=fig.add_subplot(gs[2,:3]);panel(ax,'G','Cue and memory clocks')
    for key,color,label in [('rac_zero',GREY,'Rac zero, continuous cue'),('mca_zero',PURPLE,'MCA zero, continuous cue'),('minimum_orientation_cue',BLUE,'Minimum reversing pulse')]:
        ax.plot([a['length'] for a in rr],[a[key] for a in rr],color=color,marker='s' if key=='minimum_orientation_cue' else 'o',ms=4,label=label)
    ax.set_xlabel('Cell length (µm)');
    for a in rr:ax.scatter(a['length'],a['minimum_orientation_cue'],marker='s',s=22,edgecolors=BLUE,facecolors=BLUE if a['minimum_cue_field_tail_pass'] else 'white',zorder=4)
    ax.set_ylabel('Time (s)');ax.legend(frameon=False)
    ax=fig.add_subplot(gs[2,3:]);panel(ax,'H','MCA half-mass distance')
    for key,color,label in [('initial_mca_x50_distance_from_midpoint',PURPLE,'Initial state'),('mca_x50_distance_at_minimum_cut',BLUE,'At minimum cue withdrawal')]:
        ax.plot([a['length'] for a in rr],[a[key] for a in rr],color=color,marker='s' if key=='minimum_orientation_cue' else 'o',ms=4,label=label)
    ax.set_xlabel('Cell length (µm)');ax.set_ylabel('|MCA x50 − L/2| (µm)');ax.legend(frameon=False)
    save(fig,'debelly_memory',
        'Figure 8. Rac and attachment memory through cue withdrawal. A–C: Rac polarity at L = 18, 23 and 39; D–F: sign-reversed bound-MCA polarity for the same trials, so positive values indicate the attachment bias appropriate to a new right Rac front. Blue: shortest tested reversing gradient pulse; orange: adjacent failing duration; purple: withdrawal at the first continuous-cue MCA zero. Gray time interval precedes withdrawal; time is compressed beyond ±100 s. Recovery continues for 10000 s, with the orientation gate evaluated over its final 500 s. The dotted Rac line is the 0.05 orientation threshold; contrast and tail-range criteria also apply. G: first Rac and MCA zero times during an uninterrupted gradient and minimum tested reversal pulses. G reversal-pulse squares use the same open/filled full-field flag convention. H: absolute nominal MCA half-mass distance from the midpoint, initially and at minimum-cue withdrawal. Terminal blue circles are open when the separate full-field tail check fails and filled when it passes. Gradient slope 0.24; independent endpoint commands. MCA remains rear-enriched, opposite Rac. These distribution landmarks are not tracked material particles or proof of full-field equilibration.',
        'Exact corrected basal checkpoints; frozen native replays; independently recomputed outcome gates. Piecewise uniform finite-volume half-mass interpolation. No trajectory fitting.')

def flow_figure(m,rr):
    fig,axs=plt.subplots(2,3,figsize=(12,8.4),layout='constrained')
    xx=[a['length'] for a in rr]
    for ax,letter in zip(axs.flat,'ABCDEF'):panel(ax,letter)
    ax=axs[0,0];ax.plot(xx,[a['mca_zero'] for a in rr],color=PURPLE,marker='o',label='Continuous-cue MCA zero')
    ax.plot(xx,[a['minimum_orientation_cue'] for a in rr],color=BLUE,marker='s',label='Minimum reversing pulse')
    ax.scatter([18],[rr[0]['mca_zero']],marker='x',s=70,color=ORANGE,zorder=5,label='Zero-withdrawal fails reversal')
    
    for a in rr:ax.scatter(a['length'],a['minimum_orientation_cue'],marker='s',s=22,edgecolors=BLUE,facecolors=BLUE if a['minimum_cue_field_tail_pass'] else 'white',zorder=4)
    ax.set_ylabel('Time (s)');ax.legend(frameon=False,fontsize=7)
    ax=axs[0,1]
    ax.plot(xx,[a['nominal_centroid_shift_to_zero'] for a in rr],color=PURPLE,marker='o',label='Shift to MCA zero')
    ax.plot(xx,[a['nominal_centroid_shift_under_minimum_cue'] for a in rr],color=BLUE,marker='s',label='Shift during reversing pulse')
    ax.set_ylabel('Nominal MCA centroid displacement (µm)');ax.legend(frameon=False,fontsize=7)
    ax=axs[0,2];ax.plot(xx,np.array([a['mean_net_centroid_speed_to_zero'] for a in rr])*1000,color='#009E73',marker='o');ax.set_ylabel('Mean speed to zero (10⁻³ µm/s)')
    for ax in axs[0]:ax.set_xlabel('Cell length (µm)');ax.set_xlim(17,40)
    p=n.model_for('debelly').parameters;b0=.078125;xi=p.linker_friction*b0
    z=load(FLOW/'L23_profiles.npz');refx=z['zero_x_fraction']*23;refrho=z['zero_rho_concentration']
    mesh=[]
    for L in [23,39,60,100]:
        vv=[]
        for N in [2001,4001]:
            dx=L/N;x=(np.arange(N)+.5)*dx;rho=np.interp(x,refx,refrho)
            force=p.active_stress*p.rho_to_rac_coupling*n.ds._gradient(rho,dx)
            v=-n.ds._dirichlet_elliptic(np.full(N,xi),p.cortical_viscosity,force,dx,0.,0.)
            vv.append(np.interp(5,x,v))
            if N==2001:axs[1,0].plot(np.r_[0,x,L],np.r_[0,v,0],label=f'L = {L} µm')
        rel=abs(vv[0]/vv[1]-1);assert rel<.001;mesh.append(dict(length=L,relative_mesh_change=rel))
    axs[1,0].axvline(5,color=GREY,ls=':',lw=.7);axs[1,0].set_xlabel('Distance from original left end (µm)');axs[1,0].set_ylabel('Rho-only leftward speed (µm/s)');axs[1,0].legend(frameon=False,fontsize=7)
    controls=f.rows(REV/'fixed_reference_rho_extension.csv')
    for friction,color,ls in [('True',BLUE,'-'),('False',GREY,'--')]:
        a=[r for r in controls if r['friction']==friction]
        axs[1,1].plot([float(r['length']) for r in a],[float(r['speed_at_x5']) for r in a],color=color,ls=ls,label='MCA friction' if friction=='True' else 'Zero friction')
    axs[1,1].set_xlabel('Extended domain length (µm)');axs[1,1].set_ylabel('Rho-only speed at x = 5 µm (µm/s)');axs[1,1].legend(frameon=False,fontsize=7)
    budget=[]
    for L in xx:
        if L<=21:
            slot=job_index(m,L,purpose='continuous_and_zero_withdrawal');a=load(RUN/f'workers/{slot:03d}.npz')
            tr=SimpleNamespace(x=a['x'],time=a['pulse_time'],fields={k:a['pulse_field_'+k] for k in n.DB_NAMES})
            im=int(np.flatnonzero(polarity(tr.fields['mca_bound'],tr.x,L)<=0)[0]);row,prof=ad.decompose(tr,im,L,corrected.response_cue(m['jobs'][slot]['case']));row['stage']='zero';budget.append(row)
        else:budget.append(next(r for r in read(FLOW/f'L{L}_summary.json')['rows'] if r['stage']=='zero'))
    for key,color,ls,label in [('rate_boundary',ORANGE,'-','Endpoint advection'),('rate_rho',BLUE,'-','Rho advection'),('other',GREY,'-','Other terms'),('reorientation_rate','#222222','--','Net rate')]:
        y=[sum(r['rate_'+k] for k in ['membrane','diffusion','exchange','replenishment']) if key=='other' else r[key] for r in budget]
        axs[1,2].plot(xx,np.array(y)*1e4,color=color,ls=ls,label=label)
    axs[1,2].axhline(0,color=GREY,lw=.7);axs[1,2].set_xlabel('Cell length (µm)');axs[1,2].set_ylabel('−dP(MCA)/dt (10⁻⁴ s⁻¹)');axs[1,2].legend(frameon=False,fontsize=7)
    assert max(abs(r['budget_residual']) for r in budget)<1e-11
    csvwrite(f.DATA/'zero_rate_budget.csv',budget)
    save(fig,'debelly_flow',
        'Figure 9. Length trends and distinct mechanical controls. A–C include every qualifying deBelly length from 18 to 39 µm. A compares uninterrupted-cue MCA zero with the shortest reversing pulse; the orange cross marks the failed MCA-zero withdrawal at L = 18. B shows nominal bound-MCA centroid displacement to zero and during the actual shortest reversing pulse. Minimum-pulse squares are open if the separate full-field tail check fails. Every eligible length reverses at some tested dose; no eligible “no reversal at any dose” points exist. Initial preparations at L ≤ 17 are ineligible and are excluded. C is displacement to zero divided by time to zero. D holds the local L = 23 Rho profile at MCA zero fixed in physical reference coordinates and extends its low-Rho tail; both endpoint velocities and membrane forcing are zero and basal attachment friction is uniform. E samples that static control at x = 5 µm with and without friction. D–E isolate geometry and screening, rather than evolving a new cell. F is the exact instantaneous native bound-MCA polarity-rate budget at the first MCA zero: positive terms erase the original attachment bias. Endpoint, Rho, and remaining diffusion, membrane, exchange and replenishment terms sum to the net changing-normalization rate. F is an additive identity at each native state, not a separately evolved intervention. All native trajectories use gradient slope 0.24 and independent endpoint commands.',
        dict(control_mesh_checks=mesh,native_budget_max_residual=max(abs(r['budget_residual']) for r in budget),nominal_coordinates=True))
    return dict(control_mesh=mesh,native_budget_max_residual=max(abs(r['budget_residual']) for r in budget))

# The archived dashpot generator supplies the visual grammar only.
def spring(ax,x0,x1,y,color='#4d4d4d'):
    xx=np.linspace(x0,x1,15);yy=np.full(15,y);yy[2:-2]+=np.resize([.045,-.045],11)
    ax.plot(xx,yy,color=color,lw=1.1)
def dashpot(ax,x0,x1,y,color='#4d4d4d',lw=1):
    d=x1-x0;ax.plot([x0,x0+.32*d],[y,y],color=color,lw=lw);ax.plot([x0+.62*d,x1],[y,y],color=color,lw=lw)
    ax.add_patch(Rectangle((x0+.32*d,y-.04),.30*d,.08,fill=False,ec=color,lw=lw));ax.plot([x0+.53*d,x0+.53*d],[y-.045,y+.045],color=color,lw=lw)
def vertical_dashpot(ax,x):
    ax.plot([x,x],[.13,.3],color='#d89400',lw=1);ax.add_patch(Rectangle((x-.4,.3),.8,.12,fill=False,ec='#d89400',lw=1));ax.plot([x,x],[.36,.63],color='#d89400',lw=1);ax.plot([x-.45,x+.45],[.36,.36],color='#d89400',lw=1)
def schematic(ax,L,kind,ell):
    xx=np.linspace(0,L,16)
    for a,b in zip(xx[:-1],xx[1:]):spring(ax,a,b,.66);dashpot(ax,a,b,.1)
    ax.plot(xx,np.full(16,.66),'o',ms=3,color='#4d4d4d');ax.plot(xx,np.full(16,.1),'s',ms=3,color='#2166ac')
    for x in xx:vertical_dashpot(ax,x)
    if kind=='rac':
        response=np.sinh(xx/ell)/np.sinh(L/ell)
        for x,v in zip(xx,response):
            if v>1e-3:ax.arrow(x,-.08,-3*v,0,width=.01,head_width=.07,head_length=min(.65,1.5*v),color=PURPLE,length_includes_head=True)
        ax.scatter([L],[.1],s=55,color='#b2182b',zorder=5)
    else:
        z=load(FLOW/'L23_profiles.npz');rx=z['zero_x_fraction']*23;rh=z['zero_rho_concentration']
        N=2001;dx=L/N;xf=(np.arange(N)+.5)*dx;p=n.model_for('debelly').parameters
        force=p.active_stress*p.rho_to_rac_coupling*n.ds._gradient(np.interp(xf,rx,rh),dx)
        velocity=n.ds._dirichlet_elliptic(np.full(N,p.linker_friction*.078125),p.cortical_viscosity,force,dx,0.,0.)
        for x,v in zip(xx,np.interp(xx,np.r_[0,xf,L],np.r_[0,velocity,0])):
            if abs(v)>1e-6:ax.arrow(x,-.08,200*v,0,width=.01,head_width=.07,head_length=min(.45,abs(200*v)/2),color=PURPLE,length_includes_head=True)
    ax.set_xlim(-4,73);ax.set_ylim(-.3,.86);ax.set_yticks([]);ax.set_xlabel('Nominal distance (µm)');ax.grid(False);ax.spines['left'].set_visible(False)
def mechanics_figure():
    p=n.model_for('debelly').parameters;xi=p.linker_friction*.078125;ell=np.sqrt(p.cortical_viscosity/xi)
    fig,axs=plt.subplots(3,2,figsize=(11.5,9.6),layout='constrained')
    for i,ax in enumerate(axs.flat):panel(ax,chr(65+i))
    schematic(axs[0,0],23,'rac',ell);axs[0,0].set_title('Rac end driving · L = 23 µm')
    schematic(axs[0,1],69,'rac',ell);axs[0,1].set_title('Rac end driving · L = 69 µm')
    schematic(axs[1,0],23,'rho',ell);axs[1,0].set_title('Local Rho forcing · L = 23 µm')
    schematic(axs[1,1],69,'rho',ell);axs[1,1].set_title('Local Rho forcing · L = 69 µm')
    for L in [23,39,69,100]:
        x=np.linspace(0,L,1000);distance=L-x;v=np.sinh(x/ell)/np.sinh(L/ell)
        axs[2,0].plot(distance,v,label=f'L = {L} µm')
    distance=np.linspace(0,100,1000);axs[2,0].plot(distance,np.exp(-distance/ell),color=GREY,ls='--',label='Long-domain exp(−d/ℓ)')
    axs[2,0].axvline(ell,color=GREY,ls=':',lw=.8);axs[2,0].set_xlabel('Distance from driven Rac end (µm)');axs[2,0].set_ylabel('|v| / |v at driven end|');axs[2,0].legend(frameon=False,fontsize=7)
    z=load(FLOW/'L23_profiles.npz');rx=z['zero_x_fraction']*23;rh=z['zero_rho_concentration']
    for L in [23,39,69,100]:
        N=2001;dx=L/N;x=(np.arange(N)+.5)*dx;rho=np.interp(x,rx,rh);force=p.active_stress*p.rho_to_rac_coupling*n.ds._gradient(rho,dx)
        v=-n.ds._dirichlet_elliptic(np.full(N,xi),p.cortical_viscosity,force,dx,0.,0.)
        axs[2,1].plot(np.r_[0,x,L],np.r_[0,v,0],label=f'L = {L} µm')
    axs[2,1].set_xlabel('Distance from local Rho forcing (µm)');axs[2,1].set_ylabel('Rho-only leftward speed (µm/s)');axs[2,1].legend(frameon=False,fontsize=7)
    handles=[Line2D([],[],color='#4d4d4d',label='Elastic membrane springs'),Line2D([],[],color='#2166ac',marker='s',label='Viscous cortex (axial dashpots)'),Line2D([],[],color='#d89400',label='MCA friction (vertical dashpots)'),Line2D([],[],color=PURPLE,label='Cortical velocity arrows'),Line2D([],[],marker='o',color='#b2182b',ls='',label='Rac-driven boundary')]
    fig.legend(handles=handles,loc='outside upper center',ncol=3,frameon=False,fontsize=8)
    save(fig,'debelly_screening',
        'Figure 10. Static spring/dashpot illustration of mechanical screening. A–D use 16 nodes per chain and a common physical horizontal scale for short (23 µm) and long (69 µm) cells. Upper springs represent membrane elasticity; lower axial dashpots represent cortical viscosity; vertical gold dashpots represent MCA friction. Node count is illustrative and is not a simulated coarse model: continuum coefficients are unchanged when spacing changes. A–B: a Rac-driven right endpoint moves the cortex leftward; the opposite endpoint is held at zero velocity. Arrow lengths follow the analytical uniform-friction response. C–D: arrows show the static cortical velocity produced by an identical local Rho stress-gradient forcing, with a common arrow scale between those two panels. Arrow scales for endpoint driving and Rho forcing are separate. E: the resulting Rac-boundary response is sinh(x/ℓ)/sinh(L/ℓ); d = L − x. The dashed curve is its long-domain exponential envelope. F: the same local Rho profile and low-tail extension as Figure 9D, with zero endpoint velocities. Uniform basal MCA gives ℓ = sqrt(cortical viscosity / attachment friction) = 35.8 µm. Membrane velocity is frozen at zero for these isolated controls. The upper membrane chain identifies the elastic element of the coupled native system; it is held during this screening calculation. Drawing style follows the archived deBelly dashpot figures; all calculations use the corrected current model. No animations.',
        dict(screening_length=ell,basal_bound_mca=.078125,viscosity=p.cortical_viscosity,friction=xi,analytical_endpoint_solution=True,schematic_nodes=16))
    return ell

def patch_figure(m):
    fig,axs=plt.subplots(3,2,figsize=(11.5,11.3),layout='constrained');measure=[]
    for col,geo in enumerate(['fixed_physical','fixed_relative']):
        rows=[]
        for L in [25,27,29,31,33,35]:
            z,meta=trajectory(m,L,geo);t=z['pulse_time'];x=z['x'];L=float(L);pm=polarity(z['pulse_field_mca_bound'],x,L);post=polarity(z['recovery_field_mca_bound'],x,L);D=meta['withdrawal_time'];i=np.searchsorted(t,200)
            # All selected pulses last beyond 200 s. Common-time central derivative.
            assert 0<i<len(t)-1
            speed=-L/2*(pm[i+1]-pm[i-1])/(t[i+1]-t[i-1])
            row=dict(condition=geo,length=L,duration=D,mca_at_cut=pm[-1],mca_at_200=pm[i],nominal_centroid_speed_at_200=speed,
                     nominal_centroid_shift_first_100_postcue=L/2*(post[0]-np.interp(100,z['recovery_time'],post)),
                     final_rac=meta['metrics']['final_orientation'],field_tail_pass=meta['metrics']['field_tail_pass'])
            rows.append(row);measure.append(row)
        xx=[r['length'] for r in rows]
        ax=axs[0,col];panel(ax,chr(65+col),'Fixed physical width (3 µm)' if col==0 else 'Fixed relative width (20%)')
        ax.plot(xx,[r['duration'] for r in rows],color=BLUE,label='Shortest reversing duration')
        for row in rows:ax.scatter(row['length'],row['duration'],s=25,edgecolors=BLUE,facecolors=BLUE if row['field_tail_pass'] else 'white',zorder=4)
        for L in xx:
            slot=job_index(m,L,geo,'fail');meta=read(RUN/f'workers/{slot:03d}.json');ax.scatter([L],[meta['withdrawal_time']],color=ORANGE,marker='x',s=30)
        ax.plot([],[],color=ORANGE,ls='',marker='x',label='Adjacent failing duration');ax.set_ylabel('Pulse duration (s)');ax.legend(frameon=False,fontsize=7)
        ax=axs[1,col];panel(ax,chr(67+col));ax.plot(xx,[r['mca_at_cut'] for r in rows],color=PURPLE,marker='o',label='At minimum cue withdrawal')
        ax.plot(xx,[r['mca_at_200'] for r in rows],color=GREY,marker='s',label='At common cue time 200 s');ax.axhline(0,color=GREY,lw=.7);ax.set_ylabel('MCA polarity (old rear bias)');ax.legend(frameon=False,fontsize=7)
        ax=axs[2,col];panel(ax,chr(69+col));ax.plot(xx,np.array([r['nominal_centroid_speed_at_200'] for r in rows])*1000,color=BLUE,marker='o',label='MCA centroid speed at 200 s')
        twin=ax.twinx();twin.spines['right'].set_visible(True);twin.plot(xx,[r['nominal_centroid_shift_first_100_postcue'] for r in rows],color=ORANGE,ls='--',marker='s',label='MCA shift in first 100 s after withdrawal')
        ax.set_ylabel('Centroid speed (10⁻³ µm/s)',color=BLUE);twin.set_ylabel('Post-cue MCA centroid shift (µm)',color=ORANGE)
        ax.legend(handles=ax.get_lines()+twin.get_lines(),labels=[a.get_label() for a in ax.get_lines()+twin.get_lines()],frameon=False,fontsize=6.5,loc='upper left')
        for ax in axs[:,col]:ax.set_xlabel('Cell length (µm)');ax.set_xticks(xx)
    csvwrite(f.DATA/'patch_nonmonotonicity.csv',measure)
    save(fig,'debelly_patch_turnovers',
        'Figure 11. Diagnostics of the deBelly patch-duration dips near L = 25–35 µm. Left column: fixed 3 µm patch; right: 20% relative patch. A–B show shortest tested reversing durations (open circles fail the separate full-field tail flag) and adjacent failures from the corrected campaign, exactly replayed from the same checkpoints. C–D compare MCA polarity at a common cue time (200 s) and at the shortest successful withdrawal. E–F show nominal MCA centroid speed at 200 s (left blue axis) and additional centroid displacement in the first 100 s after withdrawal (right orange axis). Positive speed and displacement erase the initial rear attachment bias. Common-time progress becomes faster even as the initial MCA bias grows. The successful cutoff is length-dependent and remains before MCA zero throughout these panels. Post-withdrawal progress also varies, with a local increase in the fixed-width case near its duration dip. These diagnostics identify competing trends and a changing post-cue commitment state; they do not establish a unique causal explanation of every local turnover. Durations are discrete tested minima, not continuous bifurcation thresholds. Recovery gates use the final 500 of 10000 s; field convergence is separately reported in the attached data.',
        'Native corrected pulse/recovery replays; central difference over one second at common cue time; centroid displacement computed from changing-normalization MCA polarity, not tracked particle transport.')
    return measure

def main():
    for p in [Path(__file__),Path(f.__file__),Path(ad.__file__),Path(corrected.__file__),Path(n.__file__),Path(n.ds.__file__)]:register(p)
    style=ROOT.parents[1]/'literature_models/debelly_2026_mechanochemical/scripts/generate_debelly_dashpot_visuals.py'
    # External-to-manuscript style is documented, not a scientific data dependency.
    style_hash=f.sha(style) if style.exists() else None
    for p in [FLOW/'independent_audit.json',REV/'controlled_calculation_audit.json']:
        a=read(p);assert a['passed']
    m,records=audit();rr=measure(m,records);memory_figure(m,rr);control=flow_figure(m,rr);ell=mechanics_figure();patch=patch_figure(m)
    stamp=datetime.now(timezone.utc).isoformat()
    result=dict(generated_at=stamp,producer='finished',passed=True,sources=f.sources,
        renderer_sha256=f.sha(Path(__file__)),native_trials=42,corrected_independent_endpoints_only=True,
        mechanics=control,screening_length=ell,archived_style_sha256=style_hash,
        figure_sha256={p.name:f.sha(p) for p in f.FIG.iterdir() if p.suffix in ['.png','.pdf','.svg']},
        table_sha256={p.name:f.sha(p) for p in f.DATA.glob('*.csv')},
        caveats=['MCA zero is a distribution landmark, not a universal successful withdrawal threshold.',
                 'All eligible lengths have a reversing tested cue. Ineligible starts are not reversal failures.',
                 'Static screening controls freeze membrane velocity; drawings do not alter model parameters.',
                 'Patch dip diagnostics are descriptive competing trends, not unique causal attribution.'])
    (f.DATA/'provenance.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f.FIG,flush=True)

if __name__=='__main__':main()
