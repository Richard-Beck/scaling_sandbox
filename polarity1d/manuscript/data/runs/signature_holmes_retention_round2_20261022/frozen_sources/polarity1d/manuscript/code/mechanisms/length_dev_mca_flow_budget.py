"""Exact MCA transport/moment budgets from retained corrected cue histories.

No native trajectories are integrated. Independent banded transport operators
recover the two pre-closure MCA pools from consecutive retained native states.
Every reconstructed pool and local/moment balance is checked against the data.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import json
import shutil

import numpy as np
from scipy.linalg import solve_banded
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from length_dev_basal_front_tension_decomposition import solve
from length_dev_output_layout import dev_output_path

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/'figures/length_mechanisms_20261008'
PARENT = ROOT/'data/runs/length_narrative_memory_20261009'
BASE = ROOT/'data/runs/debelly_corrected_three_test_20261008'
RUN = ROOT/'data/runs/debelly_mca_flow_budget_20261010'
OUT = REPORT/'dev/mca_flow_budget'
LENGTHS = [23, 25, 27, 29, 31]
STEM = 'MCA_flow_redistribution_budget'
SOURCES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source(path):
    path = Path(path)
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    return path


def read(path):
    return json.loads(source(path).read_text())


def csvwrite(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def transport(old, velocity, edge, D, decay, production, dt, dx):
    """Independently assembled implicit FV upwind/diffusion/exchange solve."""
    n = len(old)
    vf = .5*(velocity[:-1]+velocity[1:])
    pos, neg = np.maximum(vf, 0)*dt/dx, np.minimum(vf, 0)*dt/dx
    r = dt*D/dx**2
    matrix = np.zeros((3, n))
    matrix[1] = 1+dt*decay+2*r
    matrix[1, [0, -1]] -= r
    matrix[0, 1:] = -r+neg
    matrix[2, :-1] = -r-pos
    matrix[1, :-1] += pos
    matrix[1, 1:] -= neg
    matrix[1, 0] -= dt*edge[0]/dx
    matrix[1, -1] += dt*edge[1]/dx
    return solve_banded((1, 1), matrix, old+dt*production)


def flux_parts(y, velocity, edge, D, dx):
    vf = .5*(velocity[:-1]+velocity[1:])
    interior = np.r_[0., vf*np.where(vf>=0, y[:-1], y[1:]), 0.]
    left = np.zeros(len(y)+1)
    right = left.copy()
    left[0] = edge[0]*y[0]
    right[-1] = edge[1]*y[-1]
    diffusion = np.r_[0., -D*np.diff(y)/dx, 0.]
    return dict(interior_advection=interior, old_end_exit=left,
                new_end_exit=right, diffusion=diffusion)


def commands(f, p, threshold):
    switch = .5*(1+np.tanh((f['rac'][[0, -1]]-threshold)/p['protrusion_switch_scale']))
    edge = p['rac_to_rho_coupling']*p['protrusion_speed']*switch/(1+f['tension'][[0, -1]]**2/p['protrusion_tension_scale'])
    edge[1] *= -1
    return edge


def collect():
    manifest = read(PARENT/'manifest.json')
    assert manifest['center_boundary_motion'] is False
    assert read(PARENT/'status.json')['audit_passed']
    assert read(BASE/'results/independent_audit.json')['passed']
    parameter_meta = read(REPORT/'dev/resting_state/supporting/metadata/basal_endpoint_metrics_vs_length.manifest.json')
    p = parameter_meta['parameters']
    assert p['center_boundary_motion'] is False
    endpoint = source(REPORT/'dev/deBelly_redraft/supporting/data/Figure_alpha_steady_state_mechanism_endpoint_metrics.csv')
    with endpoint.open() as stream:
        thresholds = {int(r['length_um']):float(r['rac_switch_threshold']) for r in csv.DictReader(stream)}
    with gzip.open(source(BASE/'results/records.json.gz'), 'rt') as stream:
        records = {r['case_index']:r for r in json.load(stream)['cases']}
    rows, audits, snapshots, trajectories = [], [], {}, {}
    fb = p['mca_on']/(p['mca_on']+p['mca_off'])
    for L in LENGTHS:
        job = next(j for j in manifest['jobs'] if j['case']['length']==L and j['case']['condition']=='fixed_relative')
        case, D = job['case'], float(job['duration'])
        record = records[case['case_index']]
        assert D == record['result']['full_polarity']['endpoint']
        assert D == min(t['duration'] for t in record['tests'] if t['full_polarity_success'])
        assert case['response_amplitude']==3.6 and case['width_mode']=='relative'
        path = ROOT/job['source'] if job['source'] else PARENT/f"workers/{manifest['missing_jobs'].index(job):03d}.npz"
        meta = read(path.with_suffix('.json'))
        assert sha(source(path)) == meta['archive_sha256']
        checkpoint = source(ROOT/case['checkpoint'])
        assert sha(checkpoint)==case['checkpoint_sha256']
        with np.load(path) as z:
            x, t = z['x'].copy(), z['pulse_time'].copy()
            fields = {k[len('pulse_field_'):]:z[k].copy() for k in z.files if k.startswith('pulse_field_')}
            with np.load(checkpoint) as cp:
                keys = [k[6:] for k in cp.files if k.startswith('state_')]
                assert len(keys)==11
                for k in keys:
                    np.testing.assert_array_equal(fields[k][0], cp['state_'+k])
                    np.testing.assert_array_equal(fields[k][-1], z['recovery_field_'+k][0])
            tail = z['recovery_time']>=9500
            rac = z['recovery_field_rac'][tail]
            orientation = rac@(2*x/L-1)/rac.sum(axis=1)
            strength_error = float(np.max(abs(orientation/abs(case['initial_orientation'])-1)))
            assert orientation.min()>.05 and np.ptp(rac,axis=1).min()>=.05 and np.ptp(orientation)<.005
            assert strength_error<=.01 and meta['metrics']['full_polarity_success']
        assert t[0]==0 and abs(t[-1]-D)<1e-9 and np.all((np.diff(t)>0)&(np.diff(t)<=.5+1e-9))
        np.testing.assert_allclose(x, (np.arange(len(x))+.5)*L/len(x), atol=1e-12)
        dx = L/len(x)
        c = fields['mca_bound']@x/fields['mca_bound'].sum(axis=1)
        front = x>=.8*L
        total_mass = float(np.sum(fields['mca_bound'][0]+fields['mca_unbound'][0]))
        errors = dict(pool_reconstruction=0.,local_bound_budget=0.,local_total_budget=0.,centroid_budget=0.,
                      total_mass=0.,flow_superposition=0.,uniform_source=0.,flux_telescoping=0.,mechanical_transport_budget=0.,
                      independent_moment_formula=0.)
        first = {k:v[0] for k,v in fields.items()}
        initial_edge = commands(first,p,thresholds[L])
        initial_friction = p['linker_friction']*np.maximum(first['mca_bound'],1e-10)
        old_component = dict(left_drive=solve(initial_friction,p['cortical_viscosity'],np.zeros_like(x),dx,initial_edge[0],0.),
                             right_drive=solve(initial_friction,p['cortical_viscosity'],np.zeros_like(x),dx,0.,initial_edge[1]),
                             rho=solve(initial_friction,p['cortical_viscosity'],p['active_stress']*p['rho_to_rac_coupling']*np.gradient(first['rho'],dx),dx),
                             membrane=np.zeros_like(x))
        initial_remainder = first['velocity']-sum(old_component.values())
        assert abs(initial_remainder).max()<1e-5
        old_component['initial_remainder'] = initial_remainder
        for i in range(1, len(t)):
            dt = float(t[i]-t[i-1])
            old = {k:v[i-1] for k,v in fields.items()}
            new = {k:v[i] for k,v in fields.items()}
            edge = commands(old,p,thresholds[L])
            vm = (new['displacement']-old['displacement'])/dt
            me = [(new[k][0]-old[k][0])/dt for k in ['left_boundary_displacement','right_boundary_displacement']]
            bpre = transport(old['mca_bound'],old['velocity'],edge,p['mca_bound_diffusion'],p['mca_off'],p['mca_on']*old['mca_unbound'],dt,dx)
            upre = transport(old['mca_unbound'],vm,me,p['mca_unbound_diffusion'],p['mca_on'],p['mca_off']*bpre,dt,dx)
            correction = (total_mass-(bpre+upre).sum())/len(x)
            bfinal, ufinal = bpre+fb*correction, upre+(1-fb)*correction
            errors['pool_reconstruction'] = max(errors['pool_reconstruction'],float(abs(bfinal-new['mca_bound']).max()),float(abs(ufinal-new['mca_unbound']).max()))
            flux = flux_parts(bpre,old['velocity'],edge,p['mca_bound_diffusion'],dx)
            rate = {k:-np.diff(v)/dx for k,v in flux.items()}
            rate['exchange'] = p['mca_on']*old['mca_unbound']-p['mca_off']*bpre
            rate['uniform_replenishment'] = np.full(len(x),fb*correction/dt)
            actual = (new['mca_bound']-old['mca_bound'])/dt
            errors['local_bound_budget'] = max(errors['local_bound_budget'],float(abs(sum(rate.values())-actual).max()))
            fu = flux_parts(upre,vm,me,p['mca_unbound_diffusion'],dx)
            unbound_rate = -np.diff(sum(fu.values()))/dx+p['mca_off']*bpre-p['mca_on']*upre+(1-fb)*correction/dt
            observed_total = (new['mca_bound']+new['mca_unbound']-old['mca_bound']-old['mca_unbound'])/dt
            errors['local_total_budget'] = max(errors['local_total_budget'],float(abs(sum(rate.values())+unbound_rate-observed_total).max()))
            errors['total_mass'] = max(errors['total_mass'],float(abs((new['mca_bound']+new['mca_unbound']).sum()-total_mass)*dx))
            errors['uniform_source'] = max(errors['uniform_source'],float(np.ptp((new['mca_bound']-bpre)/dt)))
            total_flux = sum(flux.values())+sum(fu.values())
            errors['flux_telescoping'] = max(errors['flux_telescoping'],float(abs((-np.diff(total_flux)).sum()-(total_flux[0]-total_flux[-1]))))
            # Exact finite-step changing-mass centroid identity. Positive speed
            # means the bound-MCA centroid moves left, toward the new rear.
            weights = -(x-c[i-1])/new['mca_bound'].sum()
            speeds = {k:float(v@weights) for k,v in rate.items()}
            alternate = dict(interior_advection=-float(flux['interior_advection'][1:-1].sum()/new['mca_bound'].sum()),
                             old_end_exit=-float((x[0]-c[i-1])*flux['old_end_exit'][0]/(dx*new['mca_bound'].sum())),
                             new_end_exit=float((x[-1]-c[i-1])*flux['new_end_exit'][-1]/(dx*new['mca_bound'].sum())))
            errors['independent_moment_formula'] = max(errors['independent_moment_formula'],max(abs(alternate[k]-speeds[k]) for k in alternate))
            # Freeze donor selection on the actual full transport velocity.
            # Interior velocity is staggered (previous cortical solve); endpoint
            # transport commands use the current step's Rac/tension inputs.
            actual_vf = .5*(old['velocity'][:-1]+old['velocity'][1:])
            donors = np.where(actual_vf>=0,bpre[:-1],bpre[1:])
            mechanical_rate = {}
            for key,v in old_component.items():
                vf = .5*(v[:-1]+v[1:])
                ff = np.r_[edge[0]*bpre[0] if key=='left_drive' else 0.,vf*donors,
                           edge[1]*bpre[-1] if key=='right_drive' else 0.]
                mechanical_rate[key] = -np.diff(ff)/dx
            errors['mechanical_transport_budget'] = max(errors['mechanical_transport_budget'],
                float(abs(sum(mechanical_rate.values())-sum(rate[k] for k in ['interior_advection','old_end_exit','new_end_exit'])).max()))
            observed_speed = float((c[i-1]-c[i])/dt)
            errors['centroid_budget'] = max(errors['centroid_budget'],abs(sum(speeds.values())-observed_speed))
            friction = p['linker_friction']*np.maximum(new['mca_bound'],1e-10)
            current_split = dict(left_drive=solve(friction,p['cortical_viscosity'],np.zeros_like(x),dx,edge[0],0.),
                                 right_drive=solve(friction,p['cortical_viscosity'],np.zeros_like(x),dx,0.,edge[1]))
            component = dict(endpoint=current_split['left_drive']+current_split['right_drive'],
                             rho=solve(friction,p['cortical_viscosity'],p['active_stress']*p['rho_to_rac_coupling']*np.gradient(old['rho'],dx),dx),
                             membrane=solve(friction,p['cortical_viscosity'],friction*vm,dx))
            errors['flow_superposition'] = max(errors['flow_superposition'],float(abs(sum(component.values())-new['velocity']).max()))
            common = dict(length_um=L,time_start_s=float(t[i-1]),time_end_s=float(t[i]),step_s=dt,
                          cue_duration_s=D,centroid_um=float(c[i]),centroid_speed_um_s=observed_speed,
                          new_front_bound_MCA=float(new['mca_bound'][front].mean()),
                          new_endpoint_bound_MCA=float(new['mca_bound'][-1]),
                          mean_cortical_velocity_um_s=float(new['velocity'].mean()),
                          mean_transport_cortical_velocity_um_s=float(old['velocity'].mean()),
                          left_command_um_s=float(edge[0]),right_command_um_s=float(edge[1]),
                          right_bound_outflux=float(flux['new_end_exit'][-1]),
                          left_bound_outflux=float(-flux['old_end_exit'][0]),
                          bound_mass=float(new['mca_bound'].sum()*dx),
                          bound_uniform_source_conc_s=float(fb*correction/dt),
                          native_field_tail_pass=bool(meta['metrics']['field_tail_pass']))
            common.update({'flow_'+k+'_um_s':float(v.mean()) for k,v in component.items()})
            common.update({'speed_'+k+'_um_s':v for k,v in speeds.items()})
            common.update({'front_depletion_'+k+'_conc_s':-float(v[front].mean()) for k,v in rate.items()})
            common.update({'mechanical_speed_'+k+'_um_s':float(v@weights) for k,v in mechanical_rate.items()})
            common.update({'mechanical_front_depletion_'+k+'_conc_s':-float(v[front].mean()) for k,v in mechanical_rate.items()})
            common.update({'flow_'+k+'_um_s':float(v.mean()) for k,v in current_split.items()})
            rows.append(common)
            old_component = dict(**current_split,rho=component['rho'],membrane=component['membrane'],initial_remainder=np.zeros_like(x))
            if abs(t[i]-60)<1e-9 or abs(t[i]-180)<1e-9:
                label = f'L{L}_t{int(t[i])}'
                snapshots.update({label+'_x':x,label+'_bound':new['mca_bound'],label+'_bound_before':old['mca_bound'],
                                  label+'_velocity':new['velocity'],label+'_transport_velocity':old['velocity'],
                                  label+'_edge':edge,label+'_right_tension':np.array(old['tension'][-1])})
                snapshots.update({label+'_flow_'+k:v for k,v in component.items()})
                snapshots.update({label+'_rate_'+k:v for k,v in rate.items()})
        assert errors['pool_reconstruction']<1e-12,errors
        assert errors['flow_superposition']<1e-11,errors
        assert max(errors[k] for k in errors if k!='flow_superposition')<1e-10,errors
        trajectories.update({f'L{L}_time':t,f'L{L}_centroid':c,f'L{L}_bound':fields['mca_bound'],
                             f'L{L}_velocity':fields['velocity'],f'L{L}_tension':fields['tension'],
                             f'L{L}_x':x,f'L{L}_rac':fields['rac'],f'L{L}_rho':fields['rho']})
        audits.append(dict(length_um=L,archive=str(path.relative_to(ROOT)),checkpoint=str(checkpoint.relative_to(ROOT)),
                           cue_duration_s=D,strength_error=strength_error,native_field_tail_pass=bool(meta['metrics']['field_tail_pass']),
                           intervals=len(t)-1,errors=errors))
        print('Verified',L,errors,flush=True)
    return rows,audits,snapshots,trajectories,p


def summarize(rows):
    summary=[]
    for L in LENGTHS:
        rr=[r for r in rows if r['length_um']==L]
        for name,start,end in [('early',0.,60.),('later',60.,180.),('common',0.,180.),('whole',0.,rr[-1]['cue_duration_s'])]:
            w=np.array([max(0.,min(r['time_end_s'],end)-max(r['time_start_s'],start)) for r in rr])
            assert abs(w.sum()-(end-start))<1e-9
            row=dict(length_um=L,window=name,start_s=start,end_s=end)
            for key in rr[0]:
                if key.startswith(('speed_','flow_','front_depletion_','mechanical_')) or key in ['centroid_speed_um_s','mean_cortical_velocity_um_s','mean_transport_cortical_velocity_um_s','right_command_um_s','left_command_um_s']:
                    row[key]=float(w@np.array([r[key] for r in rr])/w.sum())
            summary.append(row)
    return summary


COLORS = dict(endpoint='#D55E00',rho='#0072B2',membrane='#999999',total='#202020',
              L23='#2675B5',L31='#DA6B20',source='#9C6E20',exchange='#009E73')


def panel(ax, letter, title):
    ax.set_title(title,fontsize=11.5,pad=12,loc='left')
    ax.text(-.16,1.04,letter,transform=ax.transAxes,fontsize=15,weight='bold')
    ax.axhline(0,color='#AAAAAA',lw=.65,zorder=0)
    ax.grid(axis='y',color='#DDDDDD',lw=.55,alpha=.65,zorder=0)
    ax.tick_params(labelsize=9)


def render(rows, summary, trajectories):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none',
                         'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7})
    fig,axes=plt.subplots(2,3,figsize=(16.4,9.6))
    fig.subplots_adjust(left=.07,right=.985,bottom=.125,top=.865,hspace=.54,wspace=.34)
    fig.text(.07,.968,'Cortical flow and MCA redistribution during reversal',fontsize=19,weight='bold',va='top')
    fig.text(.07,.923,'Corrected deBelly · 23–31 µm · relative-20% right cue · exact budgets on retained native histories',fontsize=11,color='#444444')
    xs=np.array(LENGTHS)
    def values(window,key):
        return np.array([next(r[key] for r in summary if r['length_um']==L and r['window']==window) for L in LENGTHS])
    def length_axis(ax):
        ax.set(xlabel='Cell length (µm)',xticks=LENGTHS,xlim=(22.5,31.5))
    ax=axes[0,0]
    panel(ax,'A','Endpoint commands through the early transition')
    for L in [23,31]:
        rr=[r for r in rows if r['length_um']==L and r['time_end_s']<=180]
        for key,ls in [('right_command_um_s','-'),('left_command_um_s','--')]:
            ax.plot([r['time_start_s'] for r in rr],[r[key] for r in rr],color=COLORS[f'L{L}'],ls=ls,lw=1.8)
    ax.axvspan(0,60,color='#EEEEEE',zorder=-2)
    ax.axvline(60,color='#777777',lw=.7,ls=':')
    ax.set(xlim=(0,180),xticks=[0,60,120,180],xlabel='Time after cue onset (s)',ylabel='Signed command (µm/s)')
    handles=[Line2D([],[],color=COLORS[f'L{L}'],lw=2,label=f'{L} µm') for L in [23,31]]
    handles += [Line2D([],[],color='#444444',ls=ls,label=label) for ls,label in [('-','New right end'),('--','Old left end')]]
    ax.legend(handles=handles,fontsize=8.5,frameon=False,ncol=2,loc='lower right')
    ax=axes[0,1]
    panel(ax,'B','Mean cortical flow: 60–180 s')
    for key,label,color,ls in [('flow_endpoint_um_s','Endpoint forcing',COLORS['endpoint'],'-'),
                              ('flow_rho_um_s','Rho stress gradient',COLORS['rho'],'-'),
                              ('flow_membrane_um_s','Membrane motion',COLORS['membrane'],'-'),
                              ('mean_cortical_velocity_um_s','Total',COLORS['total'],'--')]:
        ax.plot(xs,1000*values('later',key),color=color,ls=ls,marker='o',ms=3.5,lw=1.8,label=label)
    length_axis(ax)
    ax.set_ylabel('Mean velocity (nm/s; rightward +)')
    ax.legend(fontsize=8.5,frameon=False,loc='upper left')
    ax=axes[0,2]
    panel(ax,'C','Bound-MCA centroid redistribution speed')
    for window,label,color,marker in [('early','First 60 s','#009E73','^'),('later','60–180 s','#8054AA','s'),('whole','Full selected cue',COLORS['L23'],'o')]:
        ax.plot(xs,1000*values(window,'centroid_speed_um_s'),color=color,marker=marker,ms=4,lw=1.8,label=label)
    base=values('whole','centroid_speed_um_s')[0]
    ax.plot(xs,1000*base*xs/23,color='#888888',ls=':',lw=1.5,label='Proportional to length¹')
    length_axis(ax)
    ax.set_ylabel('Centroid speed (nm/s; leftward +)')
    ax.set_ylim(0,7.5)
    ax.legend(fontsize=8.5,frameon=False,loc='lower left')
    ax=axes[1,0]
    panel(ax,'D','Exact centroid budget: 60–180 s')
    keys=['new_end_exit','interior_advection','old_end_exit','uniform_replenishment','exchange','diffusion']
    labels=['Right\nexit','Interior\nadvection','Left\nexit','Uniform\nsource','Binding\nexchange','Diffusion','Net']
    xx=np.arange(7)
    for offset,L in [(-.18,23),(.18,31)]:
        r=next(r for r in summary if r['length_um']==L and r['window']=='later')
        y=np.array([r['speed_'+key+'_um_s'] for key in keys]+[r['centroid_speed_um_s']])*1000
        ax.bar(xx+offset,y,.34,color=COLORS[f'L{L}'],label=f'{L} µm',zorder=3)
    ax.set(xticks=xx,xticklabels=labels,ylabel='Centroid-speed contribution (nm/s)',ylim=(-56,64))
    ax.tick_params(axis='x',labelsize=8,length=0,pad=5)
    ax.legend(fontsize=8.5,frameon=False,loc='upper right')
    ax=axes[1,1]
    panel(ax,'E','Mechanical contributions to MCA advection²')
    for key,label,color,ls in [('right_drive','New-end forcing',COLORS['endpoint'],'-'),
                              ('rho','Rho stress gradient',COLORS['rho'],'-'),
                              ('left_drive','Old-end forcing','#8054AA',':'),
                              ('membrane','Membrane motion',COLORS['membrane'],'-')]:
        ax.plot(xs,1000*values('later','mechanical_speed_'+key+'_um_s'),color=color,ls=ls,marker='o',ms=3.5,lw=1.8,label=label)
    adv=sum(values('later','mechanical_speed_'+key+'_um_s') for key in ['right_drive','rho','left_drive','membrane','initial_remainder'])
    ax.plot(xs,1000*adv,color=COLORS['total'],ls='--',marker='o',ms=3.5,lw=1.8,label='Total advection')
    length_axis(ax)
    ax.set_ylabel('Centroid-speed contribution (nm/s)')
    ax.legend(fontsize=8.3,frameon=False,loc='center left',bbox_to_anchor=(.015,.43))
    ax.set_ylim(-4.5,13.7)
    ax=axes[1,2]
    panel(ax,'F','Bound-MCA depletion near the new front³')
    transport=sum(values('later','front_depletion_'+key+'_conc_s') for key in ['new_end_exit','interior_advection','old_end_exit','diffusion'])
    binding=values('later','front_depletion_exchange_conc_s')
    uniform=values('later','front_depletion_uniform_replenishment_conc_s')
    for y,label,color,ls in [(transport,'Net transport',COLORS['endpoint'],'-'),
                             (binding,'Binding exchange',COLORS['exchange'],'-'),
                             (uniform,'Uniform replacement',COLORS['source'],'-'),
                             (transport+binding+uniform,'Observed depletion',COLORS['total'],'--')]:
        ax.plot(xs,y*1e4,color=color,ls=ls,marker='o',ms=3.5,lw=1.8,label=label)
    length_axis(ax)
    ax.set_ylabel('Mean depletion rate (10⁻⁴ conc./s)')
    ax.legend(fontsize=8,frameon=False,loc='upper center',ncol=2,columnspacing=1.)
    ax.set_ylim(-3.2,7.4)
    fig.text(.07,.065,'¹ Anchored to the full-cue speed at 23 µm.  ² 60–180 s; source superposition at each native state, with full-flow upwind donors.',fontsize=9,color='#444444')
    fig.text(.07,.041,'³ Mean bound MCA in the rightmost ≈20% (40/201 cells), 60–180 s. Positive contributions deplete this region; negative contributions refill it.',fontsize=9,color='#444444')
    fig.text(.07,.017,'All cells remain under cue through 180 s. Full-cue averages use each length’s own shortest tested strength-passing duration.',fontsize=9,color='#444444')
    for ext in ['.png','.pdf','.svg']:
        fig.savefig(dev_output_path(OUT,STEM+ext),dpi=210,facecolor='white')
    plt.close(fig)


def verify_windows(summary,trajectories):
    checks=[]
    zpath=source(REPORT/'dev/deBelly_redraft/supporting/data/Figure_Z_MCA_distance_speed.csv')
    zmeta=read(REPORT/'dev/deBelly_redraft/supporting/metadata/Figure_Z_MCA_distance_speed.manifest.json')
    assert sha(zpath)==zmeta['sha256']['supporting/data/'+zpath.name]
    with zpath.open() as stream:
        zrows={int(r['length_um']):r for r in csv.DictReader(stream) if r['condition']=='fixed_relative'}
    for r in summary:
        L=r['length_um'];t=trajectories[f'L{L}_time'];c=trajectories[f'L{L}_centroid']
        b=trajectories[f'L{L}_bound'];x=trajectories[f'L{L}_x'];front=b[:,x>=.8*L].mean(axis=1)
        start,end=r['start_s'],r['end_s']
        actual=(float(np.interp(start,t,c))-float(np.interp(end,t,c)))/(end-start)
        actual_front=(float(np.interp(start,t,front))-float(np.interp(end,t,front)))/(end-start)
        centroid_error=abs(actual-r['centroid_speed_um_s'])
        front_error=abs(actual_front-sum(r[k] for k in r if k.startswith('front_depletion_')))
        assert centroid_error<1e-12 and front_error<1e-12
        if r['window'] in ['early','whole']:
            k='early_mean_speed_um_s' if r['window']=='early' else 'cue_average_speed_um_s'
            assert abs(actual-float(zrows[L][k]))<1e-12
        checks.append(dict(length_um=L,window=r['window'],integrated_centroid_speed_error_um_s=centroid_error,
                           integrated_front_concentration_rate_error=front_error))
    return checks


def publish(rows,summary,audits,trajectories,generated,audit):
    render(rows,summary,trajectories)
    for name in ['interval_budgets.csv','window_summary.csv','snapshots.npz','trajectories.npz']:
        shutil.copy2(RUN/name,dev_output_path(OUT,STEM+'_'+name))
    caption = '''Cortical flow and bound-MCA redistribution during corrected deBelly reversal at L = 23, 25, 27, 29, 31 µm. Each relative-20% right patch has peak 3.6 and starts at the exact original qualifying basal checkpoint, with independent endpoint commands. Selected shortest tested canonical scalar-1%-strength-passing cue durations are 232.244898, 250.000000, 270.127551, 282.755102, 285.280612 s respectively.

A: actual signed native endpoint commands for L = 23 and 31 µm, shown at the beginning of each native step; gray shading marks the first 60 s. B: time-weighted spatial-mean cortical-flow components over the common 60–180-s cue window. Endpoint forcing includes both ends. C: exact bound-MCA centroid displacement divided by elapsed time over the first 60 s, common 60–180-s window, and each selected full cue. The dotted reference is proportional to length, anchored to the full-cue speed at 23 µm; it is not a fit. B measures rightward flow; C–E use positive for leftward centroid movement toward the new rear.

D: exact changing-bound-mass centroid budget, averaged over 60–180 s, at 23 and 31 µm. Interior rightward advection moves the centroid right, whereas right-boundary loss moves it left; these large contributions cancel substantially. Uniform replenishment, binding exchange and diffusion complete the identity. The centroid remains right of the midpoint throughout these selected cues, so spatially uniform addition shifts it left toward the midpoint even while refilling the new-front region in F. The tiny old-end loss is retained.

E: partition of the advective contribution into old-end forcing, new-end forcing, Rho-gradient forcing and membrane-motion drag, using independent linear mechanical solves at each recorded state. Transport uses the previous cortical solution with the current step's prescribed endpoint velocities, preserving native staggering; all component fluxes share the donor chosen by the full native velocity. These are instantaneous flux allocations, not independently evolved pathway removals. Initial stationary reconstruction leaves an explicitly exported negligible remainder; noninitial mechanical superposition is exact. Rho's negative interior velocity contributes positively to leftward centroid movement, while opposing mean rightward flow in B.

F: exact mean bound-MCA depletion budget in cells whose reference positions are at least 0.8L (40 of 201 cells; actual width 40L/201), averaged over 60–180 s. Net transport includes flux through that region's inner face and outer boundary plus bound diffusion. Positive rates deplete the region; negative rates refill it. Binding exchange is the physical bound-pool on/off contribution; it differs from the numerical split-exchange residual in the total-MCA balance. Uniform replacement is the original implemented mass correction partitioned into the two pools; it is not a newly inferred production pathway. Its instantaneous source is spatially uniform.

B, D, E, F use a common elapsed window; the full-cue curve in C uses different selected durations. All calculations use nominal reference coordinates. Scalar Rac success and native-field tail stability remain separate; native-field flags for 23, 25, 27, 29, 31 µm are false, false, false, false, true. No new native trajectories or parameter changes. Independent implicit-transport reconstruction reproduces every saved bound and unbound update; local bound/total balances, conserved total mass, independently evaluated moment identities, full mechanical superposition, exact checkpoints, phase continuity and original scalar gates are audited. The main report is unchanged.
'''
    dev_output_path(OUT,STEM+'_caption.txt').write_text(caption)
    dev_output_path(OUT,STEM+'_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    results={}
    for window in ['early','later','common','whole']:
        a=next(r for r in summary if r['length_um']==23 and r['window']==window)
        b=next(r for r in summary if r['length_um']==31 and r['window']==window)
        results[window]=dict(mean_flow_ratio_31_over_23=b['mean_cortical_velocity_um_s']/a['mean_cortical_velocity_um_s'],
                             centroid_speed_ratio_31_over_23=b['centroid_speed_um_s']/a['centroid_speed_um_s'],
                             front_depletion_ratio_31_over_23=sum(b[k] for k in b if k.startswith('front_depletion_'))/sum(a[k] for k in a if k.startswith('front_depletion_')))
    later={r['length_um']:r for r in summary if r['window']=='later'}
    adv={L:sum(later[L]['mechanical_speed_'+key+'_um_s'] for key in ['right_drive','rho','left_drive','membrane','initial_remainder']) for L in [23,31]}
    finding = f'''# Flow versus MCA redistribution over 23–31 µm

Generated {generated}; producer finished. Native histories and the main report were preserved.

The first 60 s and later cue response differ. Across 23 to 31 µm, early mean rightward flow falls to {results['early']['mean_flow_ratio_31_over_23']:.3f} times its original value, while early bound-MCA centroid speed falls to {results['early']['centroid_speed_ratio_31_over_23']:.3f}. The longer cell retains stronger old-end forcing and initially weaker new-end forcing. These actual commands coevolve with Rac and tension; the figure does not assign independent causal shares to those feedbacks.

During the common 60–180-s window, mean rightward flow rises {results['later']['mean_flow_ratio_31_over_23']:.3f}-fold, but centroid speed rises only {results['later']['centroid_speed_ratio_31_over_23']:.3f}-fold. Rho-gradient forcing opposes the mean rightward velocity yet contributes positively to leftward MCA-centroid movement. Interior rightward advection and right-end sweep-out largely cancel, so stronger mean cortical flow does not imply proportionately faster centroid movement. The centroid tracks distribution changes, including boundary loss, replacement and binding, rather than molecular paths. Binding and uniform replacement further distinguish local depletion from centroid movement.

The advective centroid contribution itself changes only from {adv[23]*1000:.3f} to {adv[31]*1000:.3f} nm/s over that common window. Rho-associated advection supplies a positive contribution of {later[23]['mechanical_speed_rho_um_s']*1000:.3f} and {later[31]['mechanical_speed_rho_um_s']*1000:.3f} nm/s; the new-end forcing contribution is {later[23]['mechanical_speed_right_drive_um_s']*1000:.3f} and {later[31]['mechanical_speed_right_drive_um_s']*1000:.3f} nm/s. Uniform addition helps move the still-right-biased centroid toward the midpoint, while adding bound MCA back into the challenged-front region. Thus the same source helps centroid reorientation and opposes local front depletion.

Mean bound-MCA depletion in the new-front region during 60–180 s falls to {results['later']['front_depletion_ratio_31_over_23']:.3f} times its 23-µm value despite increased net transport depletion: stronger uniform refill and binding oppose that transport. The full selected-cue centroid speed rises {results['whole']['centroid_speed_ratio_31_over_23']:.3f}-fold for a length increase of {31/23:.3f}-fold. This descriptive average also includes length-dependent cue durations.

All exact native cue intervals at five lengths were reconstructed independently. The maximum bound/unbound pool reconstruction discrepancy is {max(a['errors']['pool_reconstruction'] for a in audits):.3g}; the maximum centroid-speed budget residual is {max(a['errors']['centroid_budget'] for a in audits):.3g} µm/s. No new trajectories were integrated. The mechanical partition is source superposition on the actual coupled trajectory, not an independently evolved ablation.
'''
    dev_output_path(OUT,STEM+'_findings.md').write_text(finding)
    viewer = ('<!doctype html><html lang="en"><meta charset="utf-8"><title>MCA flow and redistribution budget</title>'
              '<style>body{font:16px system-ui;max-width:1600px;margin:24px auto;padding:0 20px}img{width:100%}p{line-height:1.5}</style>'
              '<h1>Cortical flow and MCA redistribution</h1><p>Generated '+generated+' · Producer finished · Exact native budgets verified.</p>'
              '<p><a href="../vectors/'+STEM+'.pdf">Vector PDF</a> · <a href="../metadata/'+STEM+'_findings.md">Findings</a> · '
              '<a href="../metadata/'+STEM+'_audit.json">Audit</a></p><img src="../../'+STEM+'.png" alt="Six panels distinguish cortical flow, MCA centroid movement and local MCA depletion"><p>'+caption+'</p></html>')
    dev_output_path(OUT,STEM+'.html').write_text(viewer)
    (OUT/'README.md').write_text('# MCA flow and redistribution budget\n\n[Development experiments](../README.md)\n\n'
         '- [Figure]('+STEM+'.png) · [PDF](supporting/vectors/'+STEM+'.pdf) · [SVG](supporting/vectors/'+STEM+'.svg)\n'
         '- [Findings](supporting/metadata/'+STEM+'_findings.md) · [Caption](supporting/metadata/'+STEM+'_caption.txt)\n'
         '- [Reader](supporting/viewers/'+STEM+'.html) · [Audit](supporting/metadata/'+STEM+'_audit.json)\n\n'
         'Reproduce with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python -B polarity1d/manuscript/code/mechanisms/length_dev_mca_flow_budget.py` from the repository root.\n')
    artifacts = [p for p in OUT.rglob('*') if p.is_file() and not p.name.endswith('.manifest.json')]
    manifest=dict(generated_at=generated,producer='finished',passed=True,artifact_root='../..',
                  sha256={str(p.relative_to(OUT)):sha(p) for p in artifacts},source_sha256=SOURCES,
                  results=results,length_ratio_31_over_23=31/23,new_native_integrations=0,current_report_changed=False)
    dev_output_path(OUT,STEM+'.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


def main():
    generated=datetime.now(timezone.utc).isoformat()
    source(Path(__file__))
    source(Path(__file__).with_name('length_dev_basal_front_tension_decomposition.py'))
    source(Path(__file__).with_name('length_dev_output_layout.py'))
    report_hash={k:sha(REPORT/k) for k in ['index.html','length_mechanisms.pdf']}
    rows,audits,snapshots,trajectories,p=collect()
    summary=summarize(rows)
    window_checks=verify_windows(summary,trajectories)
    RUN.mkdir(parents=True,exist_ok=True)
    csvwrite(RUN/'interval_budgets.csv',rows)
    csvwrite(RUN/'window_summary.csv',summary)
    np.savez_compressed(RUN/'snapshots.npz',**snapshots)
    np.savez_compressed(RUN/'trajectories.npz',**trajectories)
    audit=dict(generated_at=generated,producer='finished',passed=True,lengths_um=LENGTHS,
               native_trajectory_integrations=0,intervals=len(rows),jobs=audits,parameters=p,
               integrated_window_checks=window_checks,figure_Z_early_and_full_cue_speeds_reproduced=True,
               current_report_unchanged=all(sha(REPORT/k)==v for k,v in report_hash.items()),
               methods='Independent implicit FV MCA reconstruction on every original cue interval; exact changing-bound-mass centroid identity; original staggered mechanical solve reconstruction. No new trajectories, interventions or fitted scaling laws.',
               source_sha256=SOURCES)
    assert audit['current_report_unchanged']
    (RUN/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    publish(rows,summary,audits,trajectories,generated,audit)
    for r in summary:
        print({k:round(v,7) if isinstance(v,float) else v for k,v in r.items() if not k.startswith('front_')},flush=True)


if __name__=='__main__':
    main()
