"""Figure 9: exact bound-MCA balance on the relative-20% L31 trajectory.

Read-only analysis of retained corrected-model evidence. A 60.5-s native
recovery replay and one basal diagnostic step validate the reconstruction.
The upwind face product rule is an allocation of the native numerical flux,
not a separately evolved ablation or a central-difference continuum estimate.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import html
import json
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from length_dev_mca_flow_budget import transport, commands
from length_dev_basal_front_tension_decomposition import solve

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'figures/length_mechanisms_20261008'
OUT = REPORT / 'revision_20261013/debelly_figure9'
PARENT = ROOT / 'data/runs/length_narrative_memory_20261009'
BASE = ROOT / 'data/runs/debelly_corrected_three_test_20261008'
STEM = 'Figure_9_bound_MCA_balance'
LENGTH = 31
SOURCES = {}
TERMS = ['profile_advection', 'flow_gradient', 'bound_diffusion',
         'net_binding', 'uniform_replacement', 'net']
COLORS = dict(profile_advection='#0072B2', flow_gradient='#D55E00',
              bound_diffusion='#7F7F7F', net_binding='#009E73',
              uniform_replacement='#CC79A7', net='#202020',
              endpoint='#D55E00', rho='#0072B2', membrane='#7F7F7F')
LABELS = dict(profile_advection='Profile advection', flow_gradient='Flow gradient',
              bound_diffusion='Bound diffusion', net_binding='Net binding',
              uniform_replacement='Uniform replacement', net='Net change')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source(path):
    path = Path(path)
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    return path


def read(path):
    return json.loads(source(path).read_text())


def csvwrite(name, records):
    path = OUT / name
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    return path


def exact_split(y, velocity, edges, diffusion, dx):
    """Exact product rule using adjacent FV faces and full-flow donors.

    Endpoint faces use the first/last cell as in the native transport solve.
    Thus the profile term vanishes in an endpoint cell whenever both its
    faces select that same cell; the allocation does not invent a boundary
    concentration gradient. Numerical upwind diffusion remains embedded in
    the donor-based terms and is distinct from physical bound diffusion.
    """
    vf = np.r_[edges[0], .5 * (velocity[:-1] + velocity[1:]), edges[1]]
    donor = np.r_[y[0], np.where(vf[1:-1] >= 0, y[:-1], y[1:]), y[-1]]
    vbar = .5 * (vf[:-1] + vf[1:])
    bbar = .5 * (donor[:-1] + donor[1:])
    profile = -vbar * np.diff(donor) / dx
    gradient = -bbar * np.diff(vf) / dx
    advective_flux = vf * donor
    diffusive_flux = np.r_[0., -diffusion * np.diff(y) / dx, 0.]
    actual_advection = -np.diff(advective_flux) / dx
    error = float(np.max(abs(profile + gradient - actual_advection)))
    assert error < 1e-13, error
    return dict(profile_advection=profile, flow_gradient=gradient,
                bound_diffusion=-np.diff(diffusive_flux) / dx,
                advective_flux=advective_flux, diffusive_flux=diffusive_flux,
                face_velocity=vf, face_donor=donor, face_mean_velocity=vbar,
                face_mean_bound=bbar, product_rule_error=error)


def interval(old, new, dt, p, threshold, x):
    dx = LENGTH / len(x)
    edge = commands(old, p, threshold)
    vm = (new['displacement'] - old['displacement']) / dt
    me = [(new[k][0] - old[k][0]) / dt for k in
          ['left_boundary_displacement', 'right_boundary_displacement']]
    bpre = transport(old['mca_bound'], old['velocity'], edge,
                     p['mca_bound_diffusion'], p['mca_off'],
                     p['mca_on'] * old['mca_unbound'], dt, dx)
    upre = transport(old['mca_unbound'], vm, me,
                     p['mca_unbound_diffusion'], p['mca_on'],
                     p['mca_off'] * bpre, dt, dx)
    total_mass = (old['mca_bound'] + old['mca_unbound']).sum()
    correction = (total_mass - (bpre + upre).sum()) / len(x)
    fb = p['mca_on'] / (p['mca_on'] + p['mca_off'])
    split = exact_split(bpre, old['velocity'], edge,
                        p['mca_bound_diffusion'], dx)
    split.update(net_binding=p['mca_on'] * old['mca_unbound'] - p['mca_off'] * bpre,
                 uniform_replacement=np.full(len(x), fb * correction / dt),
                 net=(new['mca_bound'] - old['mca_bound']) / dt,
                 preclosure_bound=bpre, preclosure_unbound=upre)
    uf = exact_split(upre, vm, me, p['mca_unbound_diffusion'], dx)
    total_flux = split['advective_flux'] + split['diffusive_flux'] + uf['advective_flux'] + uf['diffusive_flux']
    # The physical transfer rates do not cancel at different split-step times.
    # Their sum is k_on*(U_old - U_pre), a numerical residual in the total pool.
    split_residual = p['mca_on'] * (old['mca_unbound'] - upre)
    total_rate = (new['mca_bound'] + new['mca_unbound'] - old['mca_bound'] - old['mca_unbound']) / dt
    uniform_total = correction / dt
    errors = dict(
        preclosure_bound_reconstruction=float(np.max(abs(bpre + fb * correction - new['mca_bound']))),
        preclosure_unbound_reconstruction=float(np.max(abs(upre + (1-fb) * correction - new['mca_unbound']))),
        bound_balance=float(np.max(abs(sum(split[k] for k in TERMS[:-1]) - split['net']))),
        product_rule=split['product_rule_error'],
        total_balance=float(np.max(abs(-np.diff(total_flux)/dx + split_residual + uniform_total - total_rate))),
        conserved_total_mass=float(abs(total_rate.sum()*dx)),
        global_source_identity=float(abs(uniform_total*LENGTH + split_residual.sum()*dx - (total_flux[-1]-total_flux[0]))),
    )
    assert max(errors.values()) < 1e-11, errors
    split.update(total_flux=total_flux, unbound_advective_flux=uf['advective_flux'],
                 unbound_diffusive_flux=uf['diffusive_flux'],
                 total_boundary_left_outward=float(-total_flux[0]),
                 total_boundary_right_outward=float(total_flux[-1]),
                 total_uniform_replacement=float(uniform_total*LENGTH),
                 total_split_residual=float(split_residual.sum()*dx), errors=errors)
    return split


def components(state, previous, previous_dt, p, threshold, x):
    dx = LENGTH / len(x)
    if previous is None:
        force = state
        vm = np.zeros_like(x)
    else:
        force = previous
        vm = (state['displacement'] - previous['displacement']) / previous_dt
    edge = commands(force, p, threshold)
    friction = p['linker_friction'] * np.maximum(state['mca_bound'], 1e-10)
    result = dict(endpoint=solve(friction, p['cortical_viscosity'], np.zeros_like(x), dx, *edge),
                  rho=solve(friction, p['cortical_viscosity'],
                            p['active_stress'] * p['rho_to_rac_coupling'] * np.gradient(force['rho'], dx), dx),
                  membrane=solve(friction, p['cortical_viscosity'], friction * vm, dx))
    error = float(np.max(abs(sum(result.values()) - state['velocity'])))
    assert error < (1e-5 if previous is None else 1e-11), error
    return result, error


def load():
    sys.path[:0] = [str(ROOT/'code'), str(ROOT/'code/core')]
    import debelly_corrected_reversal as corrected
    n = corrected.n
    for path in [Path(__file__), Path(corrected.__file__), Path(n.__file__), Path(n.ds.__file__),
                 Path(corrected.e.__file__), ROOT/'code/core/src/polarity1d/contract.py',
                 Path(__file__).with_name('length_dev_mca_flow_budget.py'),
                 Path(__file__).with_name('length_dev_basal_front_tension_decomposition.py')]:
        source(path)
    manifest = read(PARENT/'manifest.json')
    assert manifest['center_boundary_motion'] is False
    assert read(PARENT/'status.json')['audit_passed']
    assert read(BASE/'results/independent_audit.json')['passed']
    pm = read(REPORT/'dev/resting_state/supporting/metadata/basal_endpoint_metrics_vs_length.manifest.json')
    p = pm['parameters']
    assert p['center_boundary_motion'] is False
    for key,value in p.items():
        assert getattr(n.model_for('debelly').parameters,key)==value,(key,value)
    job = next(j for j in manifest['jobs'] if j['case']['length']==LENGTH and j['case']['condition']=='fixed_relative')
    case = job['case']
    D = float(job['duration'])
    assert case['response_amplitude'] == 3.6 and case['width_mode'] == 'relative'
    with gzip.open(source(BASE/'results/records.json.gz'), 'rt') as stream:
        record = next(r for r in json.load(stream)['cases'] if r['case_index']==case['case_index'])
    assert D == record['result']['full_polarity']['endpoint']
    assert D == min(r['duration'] for r in record['tests'] if r['full_polarity_success'])
    archive = source(ROOT/job['source'] if job['source'] else PARENT/f"workers/{manifest['missing_jobs'].index(job):03d}.npz")
    am = read(archive.with_suffix('.json'))
    assert sha(archive) == am['archive_sha256']
    checkpoint = source(ROOT/case['checkpoint'])
    assert sha(checkpoint) == case['checkpoint_sha256']
    with np.load(archive) as z:
        x, tp, tr = z['x'].copy(), z['pulse_time'].copy(), z['recovery_time'].copy()
        pulse = {k:z['pulse_field_'+k].copy() for k in n.DB_NAMES}
        recovery = {k:z['recovery_field_'+k].copy() for k in n.DB_NAMES}
    with np.load(checkpoint) as z:
        for k in n.DB_NAMES:
            np.testing.assert_array_equal(pulse[k][0], z['state_'+k])
            np.testing.assert_array_equal(pulse[k][-1], recovery[k][0])
    assert tp[-1] == D and tr[-1] == 10000
    np.testing.assert_allclose(x, (np.arange(len(x))+.5)*LENGTH/len(x), atol=1e-12)
    # The recorded native switch threshold is independently retained in the
    # original basal export; do not infer it from an arbitrary snapshot.
    with source(REPORT/'dev/deBelly_redraft/supporting/data/Figure_alpha_steady_state_mechanism_endpoint_metrics.csv').open() as stream:
        threshold = float(next(r for r in csv.DictReader(stream) if int(r['length_um'])==LENGTH)['rac_switch_threshold'])
    state_at = lambda ff, i:{k:ff[k][i].copy() for k in n.DB_NAMES}
    basal = state_at(pulse, 0)
    basal_trial = n.run('debelly', LENGTH, basal, n.ConstantInput(0), .5, interval=.5)
    dense = n.run('debelly', LENGTH, state_at(recovery,0), n.ConstantInput(0), 60.5, interval=.5)
    replay_error = 0.
    for j in np.flatnonzero(tr <= 60.5):
        k = int(np.argmin(abs(dense.time-tr[j])))
        assert abs(dense.time[k]-tr[j]) < 1e-12
        for name in n.DB_NAMES:
            replay_error = max(replay_error, float(np.max(abs(dense.fields[name][k]-recovery[name][j]))))
    assert replay_error < 1e-12, replay_error
    specs = [dict(key='basal', label='Basal state', time=0., state=basal,
                  previous=None, prevdt=None, new=basal_trial.state, dt=.5, phase='basal'),
             dict(key='early', label='30 s after cue onset', time=30.,
                  state=state_at(pulse,60), previous=state_at(pulse,59), prevdt=.5,
                  new=state_at(pulse,61), dt=.5, phase='cue on'),
             dict(key='later', label='180 s after cue onset', time=180.,
                  state=state_at(pulse,360), previous=state_at(pulse,359), prevdt=.5,
                  new=state_at(pulse,361), dt=.5, phase='cue on'),
             dict(key='withdrawal', label=f'At withdrawal: {D:.2f} s', time=D,
                  state=state_at(pulse,-1), previous=state_at(pulse,-2), prevdt=float(tp[-1]-tp[-2]),
                  new=state_at(dense.fields,1), dt=.5, phase='cue off'),
             dict(key='recovery', label='60 s after withdrawal', time=D+60.,
                  state=state_at(dense.fields,120), previous=state_at(dense.fields,119), prevdt=.5,
                  new=state_at(dense.fields,121), dt=.5, phase='cue off')]
    assert tp[60] == 30. and tp[61] == 30.5
    diagnostics = []
    for spec in specs:
        spec['budget'] = interval(spec['state'], spec['new'], spec['dt'], p, threshold, x)
        spec['components'], spec['flow_error'] = components(spec['state'],spec['previous'],spec['prevdt'],p,threshold,x)
        b = spec['budget']
        gradient_rms = float(np.sqrt(np.mean(b['flow_gradient']**2)))
        diagnostics.append(dict(snapshot=spec['key'],label=spec['label'],time_after_onset_s=spec['time'],
            phase_for_next_step=spec['phase'],step_s=spec['dt'], flow_reconstruction_error=spec['flow_error'],
            **b['errors'], **{k+'_rms':float(np.sqrt(np.mean(b[k]**2))) for k in TERMS},
            **{k+'_relative_to_gradient_rms':float(np.sqrt(np.mean(b[k]**2))/gradient_rms) for k in TERMS},
            total_left_outflow=b['total_boundary_left_outward'],total_right_outflow=b['total_boundary_right_outward'],
            total_uniform_replacement=b['total_uniform_replacement'],total_split_residual=b['total_split_residual']))
    trace=[]
    maximum_errors = {k:0. for k in specs[0]['budget']['errors']}
    for i in range(1,len(tp)):
        dt = float(tp[i]-tp[i-1])
        b = interval(state_at(pulse,i-1),state_at(pulse,i),dt,p,threshold,x)
        maximum_errors = {k:max(v,b['errors'][k]) for k,v in maximum_errors.items()}
        trace.append(dict(time_start_s=float(tp[i-1]), time_end_s=float(tp[i]), step_s=dt,
                          total_left_outflow=b['total_boundary_left_outward'],
                          total_right_outflow=b['total_boundary_right_outward'],
                          total_boundary_loss=b['total_boundary_left_outward']+b['total_boundary_right_outward'],
                          total_uniform_replacement=b['total_uniform_replacement'],
                          total_split_residual=b['total_split_residual'],
                          bound_uniform_rate=float(b['uniform_replacement'][0])))
    # Compare reconstructed implicit pools and process terms with the actual
    # solver calls, separately from agreement with the saved final states.
    native_errors=[]
    original = n.ds._transport
    for spec in specs:
        calls=[]
        def traced(values, velocity, diffusion, step, dx, **kwargs):
            result=original(values,velocity,diffusion,step,dx,**kwargs)
            calls.append(result.copy())
            return result
        n.ds._transport=traced
        try:
            cue = corrected.response_cue(case) if spec['phase']=='cue on' else n.ConstantInput(0)
            trial = n.run('debelly',LENGTH,spec['state'],cue,.5,interval=.5)
        finally:
            n.ds._transport=original
        assert len(calls)==2
        error=max(float(np.max(abs(calls[0]-spec['budget']['preclosure_bound']))),
                  float(np.max(abs(calls[1]-spec['budget']['preclosure_unbound']))),
                  max(float(np.max(abs(trial.state[k]-spec['new'][k]))) for k in n.DB_NAMES))
        assert error < 1e-12,error
        native_errors.append(dict(snapshot=spec['key'],maximum_native_traced_pool_and_all_field_error=error))
    audit=dict(passed=True,length_um=LENGTH,condition='fixed_relative',width_fraction=.2,peak_extra_Rac_input=3.6,
        cue_duration_s=D,checkpoint=str(checkpoint.relative_to(ROOT)),checkpoint_sha256=sha(checkpoint),
        archive=str(archive.relative_to(ROOT)),archive_sha256=sha(archive),
        exact_initial_checkpoint_and_phase_continuity_fields=n.DB_NAMES,
        shortest_tested_strength_selection_verified=True,native_field_tail_pass=am['metrics']['field_tail_pass'],
        native_cue_intervals=len(trace),cue_maximum_errors=maximum_errors,
        recovery_replay_s=60.5,recovery_replay_all_eleven_fields_error=replay_error,
        native_traced_checks=native_errors,snapshot_diagnostics=diagnostics,
        new_long_integrations=0,diagnostic_native_steps=127)
    return x,p,specs,trace,audit


def panel(ax, letter):
    offset = -.045 if ax.get_position().width > .5 else -.14
    ax.text(offset,1.04,letter,transform=ax.transAxes,fontsize=13,weight='bold',va='bottom')
    ax.grid(axis='y',alpha=.15,lw=.65)
    ax.tick_params(labelsize=10)


def save(fig, name):
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    outside=[]
    for obj in fig.findobj(matplotlib.text.Text):
        if obj.get_visible() and obj.get_text():
            bb=obj.get_window_extent(renderer)
            if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > fig.bbox.width+1 or bb.y1 > fig.bbox.height+1:
                outside.append(obj.get_text())
    assert not outside,outside
    for ext in ['png','pdf','svg']:
        fig.savefig(OUT/(name+'.'+ext),dpi=220,facecolor='white')
    plt.close(fig)


def render(x, snapshots, trace):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                        'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(13.1,13.7))
    gs=fig.add_gridspec(4,3,left=.082,right=.985,bottom=.062,top=.935,
                         height_ratios=[1.15,.95,1.4,.95],hspace=.58,wspace=.30)
    active=snapshots[:3]
    fig.legend(handles=[Line2D([],[],color=COLORS[k],ls='--' if k=='net' else '-',lw=2,label=LABELS[k]) for k in TERMS],
               loc='upper center',bbox_to_anchor=(.535,.995),frameon=False,ncol=3,fontsize=11,columnspacing=2.2)
    bound_max=max(float(s['state']['mca_bound'].max()) for s in active)*1.08
    for j,s in enumerate(active):
        ax=fig.add_subplot(gs[0,j]);panel(ax,'ABC'[j])
        ax.plot(x,s['state']['mca_bound'],color='#222222',lw=2.4)
        ax.set(xlim=(0,LENGTH),ylim=(0,bound_max),xticks=[0,15.5,31],title=s['label'])
        ax.set_xlabel('Reference position (µm)')
        if j==0:ax.set_ylabel('Bound MCA\n(native conc.)')
        ax.axvspan(.8*LENGTH,LENGTH,color='#EEEEEE',zorder=-2)
        ax=fig.add_subplot(gs[1,j]);panel(ax,'DEF'[j]);ax.axhline(0,color='#999999',lw=.7)
        for k,label in [('endpoint','Endpoint forcing'),('rho','Rho stress gradient'),('membrane','Membrane motion')]:
            ax.plot(x,1000*s['components'][k],color=COLORS[k],lw=1.9,label=label)
        ax.plot(x,1000*s['state']['velocity'],color='#222222',ls='--',lw=2.2,label='Cortical flow')
        ax.set(xlim=(0,LENGTH),xticks=[0,15.5,31],xlabel='Reference position (µm)')
        if j==0:ax.set_ylabel('Cortical velocity\n(nm/s; rightward +)')
        ax=fig.add_subplot(gs[2,j]);panel(ax,'GHI'[j]);ax.axhline(0,color='#999999',lw=.7)
        b=s['budget']
        for k in TERMS:
            ax.plot(x,1e4*b[k],color=COLORS[k],ls='--' if k=='net' else '-',lw=2.3 if k=='net' else 1.65,zorder=5 if k=='net' else 3)
        ax.set(xlim=(0,LENGTH),xticks=[0,15.5,31],xlabel='Reference position (µm)')
        if j==0:ax.set_ylabel('Bound-MCA rate\n(10⁻⁴ native conc./s)\nsymlog outside ±1',fontsize=10)
        ax.axvspan(.8*LENGTH,LENGTH,color='#EEEEEE',zorder=-2)
        # A short endpoint segment remains visible without suppressing native
        # boundary-cell values or hiding their larger gradient allocation.
        limit=max(float(abs(b[k]).max()) for k in TERMS)*1e4*1.10
        ax.set_yscale('symlog',linthresh=1.,linscale=1.)
        ax.set_ylim(-limit,limit)
        ax.set_yticks([v for v in [-100,-10,-1,0,1,10,100] if abs(v)<limit])
        ax.set_yticklabels([str(v) for v in [-100,-10,-1,0,1,10,100] if abs(v)<limit])
    fig.legend(handles=[Line2D([],[],color=COLORS[k],lw=2,label=label) for k,label in
                        [('endpoint','Endpoint forcing'),('rho','Rho stress gradient'),('membrane','Membrane motion')]] +
                       [Line2D([],[],color='#222222',ls='--',lw=2,label='Cortical flow')],
               loc='upper center',bbox_to_anchor=(.535,gs[1,0].get_position(fig).y1+.050),
               frameon=False,ncol=4,fontsize=10,columnspacing=2.)
    ax=fig.add_subplot(gs[3,:]);panel(ax,'J')
    tt=np.array([r['time_start_s'] for r in trace])
    for k,label,color,ls in [('total_left_outflow','Left signed outward flux','#0072B2','-'),
        ('total_right_outflow','Right signed outward flux','#D55E00','-'),
        ('total_boundary_loss','Sum of both boundaries','#222222','--'),
        ('total_uniform_replacement','Uniform replacement','#CC79A7','-'),
        ('total_split_residual','Sequential split residual','#999999','--')]:
        ax.plot(tt,1e3*np.array([r[k] for r in trace]),color=color,lw=2.,ls=ls,label=label)
    ax.set(xlim=(0,trace[-1]['time_end_s']),xlabel='Time after cue onset (s)',
           ylabel='Total-pool rate\n(10⁻³ native conc.·µm/s)')
    ax.set_xticks([0,60,120,180,240,trace[-1]['time_end_s']],['0','60','120','180','240','285'])
    ax.get_xticklabels()[-1].set_ha('right')
    ax.legend(frameon=False,fontsize=9.7,ncol=3,loc='lower center',bbox_to_anchor=(.5,1.02))
    save(fig,STEM)


def render_companion(x, snapshots, trace):
    fig=plt.figure(figsize=(13.1,13.4))
    gs=fig.add_gridspec(4,3,left=.085,right=.98,bottom=.07,top=.94,hspace=.66,wspace=.36)
    fig.legend(handles=[Line2D([],[],color=COLORS[k],ls='--' if k=='net' else '-',lw=2,label=LABELS[k]) for k in TERMS],
               loc='upper center',bbox_to_anchor=(.535,.995),frameon=False,ncol=3,fontsize=11,columnspacing=2.2)
    for j,s in enumerate([snapshots[1],snapshots[3],snapshots[4]]):
        b=s['budget'];dx=LENGTH/len(x)
        face_x=np.linspace(0,LENGTH,len(x)+1)
        ax=fig.add_subplot(gs[0,j]);panel(ax,'ABC'[j])
        ax.plot(face_x,1e3*b['advective_flux'],color='#D55E00',lw=2.,label='Bound advection')
        ax.plot(face_x,1e3*b['diffusive_flux'],color='#7F7F7F',lw=1.6,label='Bound diffusion')
        ax.axhline(0,color='#AAAAAA',lw=.7)
        ax.set(title=s['label'],xlabel='Reference face position (µm)',xlim=(0,LENGTH),xticks=[0,15.5,31],ylabel='Bound-MCA flux\n(10⁻³ native conc.·µm/s)')
        if j==0:ax.legend(frameon=False,fontsize=9)
        ax=fig.add_subplot(gs[1,j]);panel(ax,'DEF'[j])
        for k in TERMS:
            ax.plot(x,1e4*b[k],color=COLORS[k],ls='--' if k=='net' else '-',lw=1.7)
        ax.set(xlim=(0,LENGTH),xticks=[0,15.5,31],xlabel='Reference position (µm)',ylabel='Bound-MCA rate\n(10⁻⁴ native conc./s)')
        ax.axhline(0,color='#AAAAAA',lw=.7)
        if s['key']=='withdrawal':
            inset=ax.inset_axes([.08,.10,.69,.50])
            for k in TERMS:
                inset.plot(x,1e4*b[k],color=COLORS[k],ls='--' if k=='net' else '-',lw=1.)
            selected=x<=30.
            limit=1.1e4*max(float(abs(b[k][selected]).max()) for k in TERMS)
            inset.set(xlim=(0,30),ylim=(-limit,limit),xticks=[0,15,30])
            inset.set_title('Bulk view: x ≤ 30 µm',fontsize=8,pad=2)
            inset.axhline(0,color='#AAAAAA',lw=.5)
            inset.tick_params(labelsize=7,pad=1)
        ax=fig.add_subplot(gs[2,j]);panel(ax,'GHI'[j])
        positions=[0,1,int(.1*len(x)),int(.9*len(x)),len(x)-2,len(x)-1]
        inflow=np.array([b['advective_flux'][k]/dx for k in positions])*1e4
        outflow=np.array([b['advective_flux'][k+1]/dx for k in positions])*1e4
        xx=np.arange(len(positions))
        ax.bar(xx-.18,inflow,.34,color='#0072B2',label='Left-face signed influx')
        ax.bar(xx+.18,-outflow,.34,color='#D55E00',label='Right-face signed influx')
        ax.plot(xx,inflow-outflow,'ko',ms=4,label='Net advection')
        ax.axhline(0,color='#AAAAAA',lw=.7)
        ax.set(xticks=xx,xticklabels=[f'{x[k]:.2f}' for k in positions],
               xlabel='Cell-center position (µm)',ylabel='Cell advective rate\n(10⁻⁴ native conc./s)')
        ax.tick_params(axis='x',labelrotation=30,labelsize=9)
    fig.legend(handles=[Line2D([],[],color='#0072B2',lw=6,label='Left-face signed influx'),
                        Line2D([],[],color='#D55E00',lw=6,label='Right-face signed influx'),
                        Line2D([],[],color='black',marker='o',ls='',label='Net advection')],
               loc='upper center',bbox_to_anchor=(.535,gs[2,0].get_position(fig).y1+.055),
               frameon=False,ncol=3,fontsize=10,columnspacing=2.)
    ax=fig.add_subplot(gs[3,:]);panel(ax,'J')
    tt=np.array([r['time_start_s'] for r in trace])
    for k,label,color,ls in [('total_left_outflow','Left signed outward flux','#0072B2','-'),
                           ('total_right_outflow','Right signed outward flux','#D55E00','-'),
                           ('total_uniform_replacement','Uniform replacement','#CC79A7','-'),
                           ('total_split_residual','Sequential split residual','#999999','--')]:
        ax.plot(tt,1e3*np.array([r[k] for r in trace]),color=color,lw=2.,ls=ls,label=label)
    ax.set(xlim=(0,trace[-1]['time_end_s']),xlabel='Time after cue onset (s)',
           ylabel='Total-pool rate\n(10⁻³ native conc.·µm/s)')
    ax.set_xticks([0,60,120,180,240,trace[-1]['time_end_s']],['0','60','120','180','240','285'])
    ax.get_xticklabels()[-1].set_ha('right')
    ax.legend(frameon=False,fontsize=10,ncol=2,loc='lower center',bbox_to_anchor=(.5,1.02))
    save(fig,STEM+'_companion')


def write_package(x, snapshots, trace, audit):
    early=snapshots[1]['budget'];basal=snapshots[0]['budget'];later=snapshots[2]['budget']
    source_fold=early['total_uniform_replacement']/basal['total_uniform_replacement']
    front=x>=.8*LENGTH
    gradient=-float(early['flow_gradient'][front].mean())
    binding=float(early['net_binding'][front].mean())
    uniform=float(early['uniform_replacement'][front].mean())
    caption=(
        '*Figure 9* Flow-gradient dilution lowers bound MCA near the new front, opposed by binding and uniform replacement. '
        'The corrected 31-µm cell receives the same relative-20% right cue as Figure 8 (extra Rac peak 3.6; shortest tested '
        '1%-strength-passing duration 285.280612 s). Columns show the basal state, 30 s and 180 s after onset. '
        'A–C) Bound-MCA profiles; shading marks the rightmost 40/201 cells. '
        'D–F) Cortical velocity and its endpoint, Rho-gradient and membrane-motion components. '
        'G–I) Exact next-0.5-s bound-MCA rate: profile advection, flow-gradient dilution/compression, physical bound diffusion, '
        'net binding, uniform replacement and their observed sum. These panels use signed symlog axes, linear between ±1 on '
        'the displayed rate scale; all native cells and values are retained. The two advection contributions use adjacent-face '
        'averages and the native upwind donor concentrations, so they sum exactly to the native advective flux divergence. '
        'J) Total-pool signed outward flux at each boundary, their sum, uniform replacement and the numerical sequential-split '
        'residual throughout the cue. Replacement plus this residual equals summed boundary loss. These global curves include '
        'both bound and unbound MCA; G–I concern bound MCA only. The final cue interval is shortened to 0.280612 s, and its '
        'rate change precedes withdrawal. The companion retains withdrawal and 60-s recovery, including the native endpoint-cell '
        'timing spike, and compares face flux with its local divergence. Cortical flow uses the preceding mechanical solve; '
        'budgets use the next native update. The mechanical curves are fixed-state source allocations, not evolved ablations. '
        'The separate native-field tail flag is true.')
    paragraphs=[
        'Before stimulation, the protrusion commands drive cortical flow toward the old left front, while the Rho stress '
        'gradient drives an opposing flow (Figure 9D). Their spatial means are −18.95 and +5.93 nm/s, leaving −13.02 nm/s '
        'net cortical motion. After each transport and exchange step, the solver restores the global MCA deficit by adding '
        'the same concentration increment at every position, partitioned between bound and unbound pools in their resting '
        'exchange proportions. This is uniform replacement. The steady bound-MCA profile is maintained during cortical '
        'motion because binding and uniform replacement offset advective loss (Figure 9G). Under the relative-20% cue, most cortical flow points right '
        'by 30 s, while bound MCA still resembles its initial right-high profile (Figure 9B,E).',
        f'The new front loses bound MCA because its transport loss exceeds its refill (Figure 9H,I). A spatial increase in '
        f'cortical velocity spreads the bound material over more cortex, giving the orange flow-gradient dilution term. At 30 s, in the rightmost '
        f'40/201 cells, flow-gradient dilution removes {gradient*1e4:.3f} × 10^-4 native concentration units/s. Net binding returns '
        f'{binding*1e4:.3f} × 10^-4 and uniform replacement returns {uniform*1e4:.3f} × 10^-4 concentration units/s, offsetting {100*binding/gradient:.0f}% and '
        f'{100*uniform/gradient:.0f}% of that loss. Profile advection and physical bound diffusion are small in this particular '
        'region at this time. The remaining net depletion is 0.586 × 10^-4 concentration units/s. The flow gradient therefore '
        'provides the principal local depletion term here; binding and replacement determine how much depletion survives.',
        'Profile advection matters elsewhere. Across the whole cell its spatial RMS is 33%, 22% and 31% of the flow-gradient '
        'term at the three displayed stages; net binding is 39%, 53% and 43%. Bound diffusion is below 0.5% in each. Those '
        'whole-cell magnitudes do not describe every location. At 180 s the MCA profile declines toward the right endpoint; '
        'rightward flow brings higher MCA into that region, so the blue profile-advection term partly restores it. An outward '
        'endpoint flux also does not specify local depletion. A cell can receive '
        'almost the same flux through its other face, making the difference small. The companion shows both face contributions '
        'and their local sum, including the endpoint cells.',
        f'Uniform replacement is {source_fold:.2f} times its basal rate at 30 s (Figure 9J). The total-pool rate is 0.004577 '
        'native concentration·µm/s, while combined boundary loss is 0.007075. The remaining 0.002499 is the sequential-exchange '
        'residual, 35.3% of boundary loss: bound and unbound exchange are evaluated at different stages of the native update. '
        'This residual is numerical bookkeeping, not an additional biological source. The global relation therefore includes '
        'it explicitly. Physical net binding remains a substantial term in the bound-MCA balance and cannot be inferred from '
        'a bound-plus-unbound budget, where physical exchange cancels. The cue raises boundary throughput and uniform refill '
        'together; local MCA depletion reflects their spatially unequal effects.'
    ]
    methods='''Figure 9: exact bound-MCA process allocation, 2026-10-13

Evidence and trajectory. The original corrected independent-endpoint L=31 µm relative-20% basal checkpoint and its retained shortest tested scalar-1%-strength-passing pulse are unchanged. The pulse has native extra Rac peak 3.6 and duration 285.28061224489795 s. All eleven initial fields match the checkpoint exactly and all eleven fields are continuous at withdrawal. Canonical selection is independently checked against the corrected campaign records. The native-field tail flag is separately true. Nominal undeformed coordinates use 201 cells, dx=31/201 µm. All concentrations and concentration-integrated amounts retain native units; no experimental concentration calibration is implied.

Native bound update. Reconstruct the pre-correction bound field B_pre with the independently assembled implicit finite-volume operator, using the previous stored cortical velocity in the interior, endpoint commands from the current Rac/tension state, implicit physical bound diffusion and unbinding, and binding source k_on*U_old. Reconstruct U_pre with the newly calculated membrane velocity, implicit unbound diffusion/unbinding, and source k_off*B_pre. If c is the per-cell total-pool deficit, the spatially uniform correction is f_bound*c for bound and (1-f_bound)*c for unbound, where f_bound=k_on/(k_on+k_off).

Exact advection allocation. At internal faces v_face=(v_i+v_i+1)/2 and B_face is B_pre from the upwind donor selected by the full face velocity; endpoint faces use the first/last B_pre as in the native operator. For each cell, let v_bar and B_bar be arithmetic averages of its left/right face values. Profile advection is -v_bar*(B_face_right-B_face_left)/dx. Flow-gradient dilution/compression is -B_bar*(v_face_right-v_face_left)/dx. The exact two-face product rule gives their sum as -(v_face_right*B_face_right-v_face_left*B_face_left)/dx. These are the discrete analogues of -v*dB/dx and -B*dv/dx, not central-difference estimates. Numerical upwind smoothing remains embedded in these donor-based terms; the separate bound-diffusion curve is only the explicit physical diffusion coefficient. The split is a native-flux allocation, not a proof that the two contributions can be causally perturbed independently.

Other bound contributions. Bound diffusion is the divergence of the no-flux diffusive face flux using B_pre. Net binding is k_on*U_old-k_off*B_pre, which is the physical exchange contribution in the native bound solve. Uniform replacement is f_bound*c/dt at every cell. Their sum with both advective terms reproduces (B_new-B_old)/dt. Gross binding and unbinding are each larger than their difference; the displayed net binding is the rate that changes bound MCA.

Mechanical timing. The displayed cortical velocity belongs to the completed preceding mechanical solve. Its exact source decomposition uses the completed-step bound-MCA drag, the preceding Rho profile, that step's membrane velocity and its endpoint commands. At the basal checkpoint the stationary reconstruction uses zero membrane velocity. Bound transport in the next step uses this stored cortical velocity in the interior, but the next step's endpoint commands from the displayed current Rac/tension. The plotted source decomposition and the bound flux therefore have the intentional native staggered timing. They are not independent pathway ablations. Endpoint motion of the membrane is distinct from imposed cortical endpoint velocity.

Flux and divergence. Positive face flux is rightward. The advective contribution to cell i is (left-face flux-right-face flux)/dx. A negative left endpoint flux is outward loss, but the adjacent internal face can replenish nearly the same amount from the next cell. The companion plots the two signed cell-face contributions and their sum at six specified native cell centers. It does not infer local transport depletion from endpoint flux magnitude alone.

Total-pool identity. Let Q_L=-F_total_left, Q_R=F_total_right, S=L*c/dt, and E=dx*sum(k_on*(U_old-U_pre)). Total includes both pools and both advection/diffusion fluxes. S+E=Q_L+Q_R exactly, because final total MCA is restored. E is the numerical sequential-exchange residual: the physical binding transfers would cancel if evaluated at the same time, but the staggered bound and unbound solves use different unbound fields. S alone is consequently smaller than boundary loss in these snapshots. The individual outward boundary fluxes are signed: the small negative late left flux is unbound inward transport and is retained, not clamped. Bound-only replacement in G-I is f_bound*c/dt, not the integrated total S in J.

Sampling and endpoint cells. Main snapshots are basal, 30 s and 180 s after onset; all plotted process rates are exact next-0.5-s increments. Signed symlog rate axes are linear for magnitude <=1 in units of 10^-4 native concentration/s and logarithmic outside that band. Axis ranges vary by panel; no boundary cell is omitted or clipped. The companion gives 30 s, withdrawal and 60 s after withdrawal, retaining the full native endpoint values on linear axes. Withdrawal fields are the inherited cue-on endpoint, while its rate panel uses the next cue-off 0.5-s update. At withdrawal, the next right cortical endpoint command is 0.124384 µm/s while the retained adjacent interior-face speed is 0.092406 µm/s. This native timing difference yields a large last-cell flow-gradient rate; it is not evidence for an additional localized biochemical mechanism. The last cue step is shortened from 0.5 to 0.280612 s, and the original endpoint relaxation is applied once per step. The sharp final rate change in the cue trace occurs during that shortened cue-on interval, before withdrawal. No timestep or mesh convergence claim is made from this accounting. The main 180-s comparison avoids making that endpoint timing transient the principal mechanism example; the companion preserves it for inspection.

Verification. All 571 retained cue intervals reconstruct both pre-correction MCA fields, both final pools, the six-term bound identity and the total-pool/global identities. Five independent one-step native traces capture the actual two transport solves and match the reconstruction and all eleven final fields to roundoff. A 60.5-s native recovery replay matches all eleven original fields at every retained 2.5-s recovery sample. No new long integration, altered parameter, source input, or original artifact is used. Spatial RMS ratios are unweighted over the 201 equal-volume cells; right-front averages use the last 40 cells, width 40L/201. Whole-cell RMS ratios are descriptive magnitudes and are not local cancellation fractions.
'''
    (OUT/(STEM+'_caption.txt')).write_text(caption+'\n')
    (OUT/(STEM+'_methods.txt')).write_text(methods)
    companion_caption=(
        'Figure 9 companion. A–C) Exact next-step bound-MCA advective and physical diffusive face fluxes at 30 s, withdrawal '
        'and 60 s after withdrawal. D–F) The complete bound-MCA process balance on linear axes; panel E preserves the native '
        'rightmost-cell spike at withdrawal, with an inset resolving the bulk region x ≤ 30 µm. G–I) Signed contributions from each cell face at six explicit cell centers, and their '
        'sum. The left-face contribution is +F_left/dx and the right-face contribution is −F_right/dx; they can be large and '
        'opposing while net advection is small. J) Signed total-pool boundary fluxes, uniform replacement and split residual '
        'during the cue. See the main caption and methods for donor convention, physical versus numerical exchange, the shortened '
        'last cue interval and endpoint/interior velocity timing. None of the endpoint cells or flux signs is suppressed.')
    (OUT/(STEM+'_companion_caption.txt')).write_text(companion_caption+'\n')
    chapter=dict(id='figure9',title='Flow-gradient dilution competes with binding and uniform replacement',
        image=str((OUT/(STEM+'.png')).relative_to(REPORT)),paragraphs=paragraphs,caption=caption,
        manifest=str((OUT/(STEM+'.manifest.json')).relative_to(REPORT)),sources=[
            [str((OUT/(STEM+'_methods.txt')).relative_to(REPORT)),'Exact bound-MCA methods'],
            [str((OUT/(STEM+'_companion.pdf')).relative_to(REPORT)),'Flux, endpoint cells and recovery companion'],
            [str((OUT/(STEM+'_audit.json')).relative_to(REPORT)),'Native-step audit']])
    (OUT/'chapter.json').write_text(json.dumps(chapter,indent=2,ensure_ascii=False)+'\n')
    generated=datetime.now(timezone.utc).isoformat()
    body='<!doctype html><meta charset="utf-8"><title>Bound-MCA balance</title><style>body{font:18px/1.5 sans-serif;max-width:1150px;margin:2rem auto;padding:0 1rem}img{width:100%}a{color:#17609a}</style>'
    body+='<h1>Bound-MCA balance in the corrected 31-µm cell</h1><p>Generated '+generated+'; producer finished.</p>'
    body+=''.join('<p>'+html.escape(p)+'</p>' for p in paragraphs)
    body+='<img src="'+STEM+'.png" alt="Bound MCA, cortical flow, six local process terms and global pool accounting"><p>'+html.escape(caption).replace('*Figure 9*','<em>Figure 9</em>')+'</p>'
    body+='<p><a href="'+STEM+'_methods.txt">Exact methods</a> · <a href="'+STEM+'.pdf">Main vector PDF</a> · <a href="'+STEM+'_audit.json">Native-step audit</a></p>'
    body+='<h2>Flux, native endpoint cells and recovery</h2><img src="'+STEM+'_companion.png" alt="Native face flux, linear local process balance and opposing cell-face rates"><p>'+html.escape(companion_caption)+'</p>'
    (OUT/'index.html').write_text(body+'\n')
    artifacts={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name != STEM+'.manifest.json'}
    manifest=dict(generated_at=generated,producer='finished',passed=True,
        renderer=str(Path(__file__).relative_to(ROOT)),renderer_sha256=sha(__file__),
        artifacts=artifacts,sha256=artifacts,sources=SOURCES,source_sha256=SOURCES,
        caption=caption,audit=audit,new_long_integrations=0,
        reproduction='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python -B polarity1d/manuscript/code/mechanisms/length_revision_20261013_figure9.py')
    (OUT/(STEM+'.manifest.json')).write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    x,p,snapshots,trace,audit=load()
    profiles=[];faces=[];regions=[]
    for s in snapshots:
        b=s['budget']
        for i,xx in enumerate(x):
            profiles.append(dict(snapshot=s['key'],label=s['label'],time_after_onset_s=s['time'],
                reference_position_um=float(xx),bound_MCA=float(s['state']['mca_bound'][i]),
                unbound_MCA_before=float(s['state']['mca_unbound'][i]),
                bound_MCA_after=float(s['new']['mca_bound'][i]),step_s=s['dt'],
                bound_preclosure=float(b['preclosure_bound'][i]),
                unbound_preclosure=float(b['preclosure_unbound'][i]),
                **{k:float(b[k][i]) for k in TERMS},
                **{k+'_velocity_um_s':float(v[i]) for k,v in s['components'].items()},
                cortical_velocity_um_s=float(s['state']['velocity'][i]),
                face_mean_velocity_um_s=float(b['face_mean_velocity'][i]),face_mean_bound_MCA=float(b['face_mean_bound'][i])))
        for i,xx in enumerate(np.linspace(0,LENGTH,len(x)+1)):
            faces.append(dict(snapshot=s['key'],reference_face_position_um=float(xx),
                face_velocity_um_s=float(b['face_velocity'][i]),upwind_bound_donor=float(b['face_donor'][i]),
                bound_advective_flux=float(b['advective_flux'][i]),bound_diffusive_flux=float(b['diffusive_flux'][i]),
                unbound_advective_flux=float(b['unbound_advective_flux'][i]),unbound_diffusive_flux=float(b['unbound_diffusive_flux'][i]),
                total_pool_flux=float(b['total_flux'][i])))
        for name,mask in [('old_front_20pct',x<=.2*LENGTH),('new_front_20pct',x>=.8*LENGTH),
                          ('first_cell',np.arange(len(x))==0),('last_cell',np.arange(len(x))==len(x)-1)]:
            regions.append(dict(snapshot=s['key'],region=name,cells=int(mask.sum()),
                                **{k:float(b[k][mask].mean()) for k in TERMS}))
    csvwrite(STEM+'_profiles_and_rates.csv',profiles)
    csvwrite(STEM+'_face_fluxes.csv',faces)
    csvwrite(STEM+'_regional_rates.csv',regions)
    csvwrite(STEM+'_global_pool_rates.csv',trace)
    csvwrite(STEM+'_snapshot_diagnostics.csv',audit['snapshot_diagnostics'])
    render(x,snapshots,trace)
    render_companion(x,snapshots,trace)
    # Numerical findings are exported before drafting the final figure prose.
    (OUT/(STEM+'_audit.json')).write_text(json.dumps(audit,indent=2)+'\n')
    (OUT/(STEM+'_parameters.json')).write_text(json.dumps(p,indent=2)+'\n')
    (OUT/(STEM+'_sources.json')).write_text(json.dumps(SOURCES,indent=2)+'\n')
    write_package(x,snapshots,trace,audit)
    print(json.dumps(audit,indent=2))


if __name__=='__main__':
    main()
