"""Figure Y: four matched mechanical/flow/total-MCA-flux snapshots.

Retained strength-selected pulse plus short, independently checked dense
recovery replay; original report and prior Figure Y drafts are preserved.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import FancyBboxPatch

from length_dev_basal_front_tension_decomposition import solve, gradient
import length_dev_output_layout as layout

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'figures/length_mechanisms_20261008'
OUT = REPORT / 'dev/resting_state'
BASE = ROOT / 'data/runs/debelly_corrected_three_test_20261008'
LENGTH = 31
SLIP_WINDOW = 10.
EXAGGERATION = 2.5
COLORS = dict(endpoint='#D55E00', rho='#0072B2', membrane='#929292', total='#222222',
              mca='#009E73', cortex='#762A83', rac='#C13C25', tension='#9D245D', diffusion='#009E73',
              unbound_advection='#CC79A7', unbound_diffusion='#0072B2')
FIELDS = ['rac', 'rho', 'mca_bound', 'mca_unbound', 'displacement', 'velocity', 'tension',
          'raw_left_displacement', 'raw_right_displacement', 'left_boundary_displacement', 'right_boundary_displacement']
SOURCES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def register(path):
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    return path


def read(path):
    return json.loads(register(path).read_text())


def output(name):
    return layout.dev_output_path(OUT, name)


def csvwrite(name, records):
    with output(name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def commands(fields, parameters):
    rac, tension = fields['rac'][..., [0, -1]], fields['tension'][..., [0, -1]]
    switch = .5*(1+np.tanh((rac-.899061648345644)/parameters['protrusion_switch_scale']))
    result = parameters['rac_to_rho_coupling']*parameters['protrusion_speed']*switch/(1+tension**2/parameters['protrusion_tension_scale'])
    result[..., 1] *= -1
    return result


def window_integral(values, times, start, stop):
    """Trapezoid quadrature with interpolated exact window endpoints."""
    selected = (times > start) & (times < stop)
    interp = lambda when: np.array([np.interp(when, times, values[:, j]) for j in range(values.shape[1])])
    yy = np.vstack([interp(start), values[selected], interp(stop)])
    tt = np.r_[start, times[selected], stop]
    return np.trapezoid(yy, tt, axis=0), interp(start)


def snapshot(fields, times, i, x, parameters, label, when, basal=False):
    f = {key: value[i].copy() for key, value in fields.items()}
    dx = LENGTH/len(x)
    friction = parameters['linker_friction']*f['mca_bound']
    if basal:
        force_rho = f['rho']
        boundary = commands(f, parameters)
        membrane_v = np.zeros_like(x)
        membrane_boundary_v = np.zeros(2)
        slip = SLIP_WINDOW*f['velocity']
        endpoint_slip = SLIP_WINDOW*boundary
    else:
        dt = times[i]-times[i-1]
        assert 0 < dt <= .5+1e-9
        start = times[i]-SLIP_WINDOW
        assert start >= times[1]
        force_rho = fields['rho'][i-1]
        boundary = commands({key: value[i-1] for key, value in fields.items()}, parameters)
        membrane_v = (f['displacement']-fields['displacement'][i-1])/dt
        integral_v, _ = window_integral(fields['velocity'][:i+1], times[:i+1], start, times[i])
        _, u_start = window_integral(fields['displacement'][:i+1], times[:i+1], start, times[i])
        slip = integral_v-(f['displacement']-u_start)
        boundary_history = commands({key: value[:i] for key, value in fields.items()}, parameters)
        integral_boundary, _ = window_integral(boundary_history, times[1:i+1], start, times[i])
        endpoint_u = np.column_stack([fields[key][:i+1,0] for key in ['left_boundary_displacement','right_boundary_displacement']])
        _, endpoint_start = window_integral(endpoint_u, times[:i+1], start, times[i])
        endpoint_slip = integral_boundary-(endpoint_u[-1]-endpoint_start)
        membrane_boundary_v = (endpoint_u[-1]-endpoint_u[-2])/dt
    components = dict(
        endpoint=solve(friction, parameters['cortical_viscosity'], np.zeros_like(x), dx, *boundary),
        rho=solve(friction, parameters['cortical_viscosity'], parameters['active_stress']*parameters['rho_to_rac_coupling']*np.gradient(force_rho, dx), dx),
        membrane=solve(friction, parameters['cortical_viscosity'], friction*membrane_v, dx))
    residual = float(np.max(abs(sum(components.values())-f['velocity'])))
    assert residual < (1e-5 if basal else 1e-11), (label, residual)
    tension_error = float(np.max(abs(parameters['membrane_stiffness']*gradient(f['displacement'], dx, f['left_boundary_displacement'][0], f['right_boundary_displacement'][0])-f['tension'])))
    assert tension_error < 1e-9
    # Instantaneous native face flux: positive is rightward. These displayed
    # state fluxes are not the pre-closure split-step source budget.
    face_v = .5*(f['velocity'][:-1]+f['velocity'][1:])
    donors = np.where(face_v >= 0, f['mca_bound'][:-1], f['mca_bound'][1:])
    advective = np.r_[boundary[0]*f['mca_bound'][0], face_v*donors, boundary[1]*f['mca_bound'][-1]]
    diffusive = np.r_[0., -parameters['mca_bound_diffusion']*np.diff(f['mca_bound'])/dx, 0.]
    unbound_face_v = .5*(membrane_v[:-1]+membrane_v[1:])
    unbound_donors = np.where(unbound_face_v >= 0, f['mca_unbound'][:-1], f['mca_unbound'][1:])
    unbound_advective = np.r_[membrane_boundary_v[0]*f['mca_unbound'][0], unbound_face_v*unbound_donors, membrane_boundary_v[1]*f['mca_unbound'][-1]]
    unbound_diffusive = np.r_[0., -parameters['mca_unbound_diffusion']*np.diff(f['mca_unbound'])/dx, 0.]
    total_flux = advective+diffusive+unbound_advective+unbound_diffusive
    divergence = -np.diff(total_flux)/dx
    conservation_residual = float(abs(divergence.sum()*dx-(total_flux[0]-total_flux[-1])))
    assert conservation_residual < 1e-12
    return dict(label=label, time=when, f=f, boundary=boundary, slip=slip, endpoint_slip=endpoint_slip,
        components=components, membrane_v=membrane_v, membrane_boundary_v=membrane_boundary_v,
        face_x=np.linspace(0, LENGTH, len(x)+1), advective=advective, diffusive=diffusive,
        unbound_advective=unbound_advective, unbound_diffusive=unbound_diffusive, total_flux=total_flux,
        audit=dict(label=label, velocity_reconstruction_residual_um_per_s=residual, native_tension_identity_residual_pN=tension_error, flux_conservation_residual=conservation_residual,
                   max_abs_advective_flux=float(abs(advective).max()), max_abs_diffusive_flux=float(abs(diffusive).max()),
                   max_abs_unbound_advective_flux=float(abs(unbound_advective).max()), max_abs_unbound_diffusive_flux=float(abs(unbound_diffusive).max())))


def draw_dashpot(ax, start, end, color, width=1., halfwidth_pt=2.5):
    """Draw barrel and piston normal to the actual rendered connector."""
    start, end = np.asarray(start), np.asarray(end)
    ds, de = ax.transData.transform([start, end])
    delta = de-ds
    normal = np.array([-delta[1], delta[0]])/np.linalg.norm(delta)*halfwidth_pt*ax.figure.dpi/72
    a, b, piston = ds+.32*delta, ds+.69*delta, ds+.50*delta
    inverse = ax.transData.inverted()
    lines = [np.array([ds, a]), np.array([a+normal, b+normal, b-normal, a-normal]),
             np.array([piston-normal, piston+normal]), np.array([piston, de])]
    for line in lines:
        ax.plot(*inverse.transform(line).T, color=color, lw=width, solid_capstyle='round', zorder=3)


def draw_spring(ax, left, right, y, color):
    xx = np.r_[left, np.linspace(left+.13*(right-left), right-.13*(right-left), 13), right]
    yy = np.full_like(xx, y)
    yy[2:-2] += .075*np.resize([1, -1], len(xx)-4)
    ax.plot(xx, yy, color=color, lw=1.6, zorder=3)


def title(ax, letter, text):
    ax.set_title(letter+'  '+text, loc='left', fontsize=12, weight='bold', pad=10)


def main():
    generated = datetime.now(timezone.utc).isoformat()
    stem = 'figure_Y_cortical_motion_mca_12panel'
    report_hashes = {name: sha(REPORT/name) for name in ['index.html', 'length_mechanisms.pdf']}
    for path in [Path(__file__), Path(layout.__file__), Path(__file__).with_name('length_dev_basal_front_tension_decomposition.py')]:
        register(path)
    meta = read(OUT/'supporting/metadata/basal_endpoint_metrics_vs_length.manifest.json')
    p = meta['parameters']
    assert not p['center_boundary_motion']
    base = read(BASE/'manifest.json')
    assert not base['center_boundary_motion'] and read(BASE/'results/independent_audit.json')['passed']
    case = next(c for c in base['cases'] if c['model']=='debelly' and c['length']==LENGTH and c['condition']=='fixed_slope')
    checkpoint = register(ROOT/case['checkpoint'])
    assert sha(checkpoint) == case['checkpoint_sha256']
    with np.load(checkpoint) as z:
        x = z['x'].copy()
        initial = {key: z['state_'+key][None, :].copy() for key in FIELDS}
    snapshots = [snapshot(initial, np.array([0.]), 0, x, p, 'Basal steady state', 0., basal=True)]
    sys.path[:0] = [str(ROOT/'code'), str(ROOT/'code/core')]
    import debelly_corrected_reversal as corrected
    n = corrected.n
    for path in [Path(corrected.__file__), Path(n.__file__), Path(n.ds.__file__)]:
        register(path)
    assert not n.model_for('debelly').parameters.center_boundary_motion
    run = ROOT/'data/runs/length_narrative_memory_20261009'
    selection = read(run/'manifest.json')
    job = next(j for j in selection['jobs'] if j['case']['length']==LENGTH and j['case']['condition']=='fixed_slope')
    assert job['case']['checkpoint_sha256'] == sha(checkpoint)
    assert job['expected']['full_polarity_success']
    source = register(ROOT/job['source'] if job['source'] else run/f"workers/{selection['missing_jobs'].index(job):03d}.npz")
    source_meta = read(source.with_suffix('.json'))
    assert sha(source)==source_meta['archive_sha256']
    with np.load(source) as z:
        t = z['pulse_time'].copy()
        fields = {key: z['pulse_field_'+key].copy() for key in FIELDS}
        rt = z['recovery_time'].copy()
        recovery = {key: z['recovery_field_'+key].copy() for key in n.DB_NAMES}
    for key in FIELDS:
        np.testing.assert_array_equal(fields[key][0], initial[key][0])
        np.testing.assert_array_equal(fields[key][-1], recovery[key][0])
    assert t[-1] == job['duration']
    i = int(np.argmin(abs(t-30.)))
    assert t[i] == 30.
    snapshots.append(snapshot(fields, t, i, x, p, 'Early cue: 30 s after onset', 30.))
    snapshots.append(snapshot(fields, t, len(t)-1, x, p, f'Cue off: {t[-1]:.2f} s after onset', t[-1]))

    # Only 80 s / 160 native steps: audit a dense local recovery replay against
    # every retained 2.5-s sample, then refine the requested event-based time.
    start, finish = 40., 120.
    start_i = int(np.argmin(abs(rt-start)))
    assert rt[start_i] == start
    state = {key: value[start_i].copy() for key, value in recovery.items()}
    replay = n.run('debelly', LENGTH, state, n.ConstantInput(0), finish-start, interval=.5)
    native_rt = start+replay.time
    restart_error = 0.
    for j in np.where((rt>=start)&(rt<=finish))[0]:
        k = int(np.argmin(abs(native_rt-rt[j])))
        assert abs(native_rt[k]-rt[j]) < 1e-9
        restart_error = max(restart_error, max(float(np.max(abs(replay.fields[key][k]-recovery[key][j]))) for key in n.DB_NAMES))
    assert restart_error < 1e-9
    threshold = .899061648345644
    incumbent = replay.fields['rac'][:,0]
    crosses = np.where((incumbent[:-1]>=threshold)&(incumbent[1:]<threshold))[0]
    assert len(crosses)==1
    k = int(crosses[0])
    assert np.all(recovery['rac'][rt<=start,0] >= threshold)
    crossing = native_rt[k]+(native_rt[k+1]-native_rt[k])*(incumbent[k]-threshold)/(incumbent[k]-incumbent[k+1])
    target = float(crossing+60.)
    assert target < finish
    exact = n.run('debelly', LENGTH, state, n.ConstantInput(0), target-start, interval=.5)
    local = {key: exact.fields[key].copy() for key in FIELDS}
    np.testing.assert_array_equal(exact.time[:-1], replay.time[:len(exact.time)-1])
    for key in n.DB_NAMES:
        np.testing.assert_array_equal(exact.fields[key][:-1], replay.fields[key][:len(exact.time)-1])
    snapshots.append(snapshot(local, start+exact.time, len(exact.time)-1, x, p, '60 s after old-front Rac switch-off', job['duration']+target))
    replay_audits = [dict(replay_start_after_withdrawal_s=start, replay_end_after_withdrawal_s=finish,
                         maximum_all_retained_state_replay_error=restart_error, final_snapshot_short_step_s=float(np.diff(exact.time)[-1]))]
    event = dict(endpoint='leftmost native Rac sample (incumbent front)', threshold=threshold,
                 crossing_bracket_after_withdrawal_s=[float(native_rt[k]),float(native_rt[k+1])],
                 interpolated_crossing_after_withdrawal_s=float(crossing), final_snapshot_after_withdrawal_s=target,
                 final_snapshot_after_onset_s=job['duration']+target, time_after_crossing_s=60.)
    csvwrite(stem+'_event_trace.csv', [dict(time_after_withdrawal_s=float(tt), incumbent_rac=float(rr), threshold=threshold) for tt,rr in zip(native_rt,incumbent)])
    clock_description = f'Finite constant-slope cue: off at {job["duration"]:.2f} s; final row at {target:.2f} s after withdrawal.'

    nodes = np.linspace(0, LENGTH, 17)
    data_rows, flux_rows, node_rows = [], [], []
    for s in snapshots:
        f = s['f']
        left, right = f['left_boundary_displacement'][0], f['right_boundary_displacement'][0]
        u = np.interp(nodes, np.r_[0, x, LENGTH], np.r_[left, f['displacement'], right])
        slip = np.interp(nodes, np.r_[0, x, LENGTH], np.r_[s['endpoint_slip'][0], s['slip'], s['endpoint_slip'][1]])
        s['nodes'], s['membrane_x'], s['cortex_x'] = nodes, nodes+EXAGGERATION*u, nodes+EXAGGERATION*(u+slip)
        assert np.all(np.diff(s['membrane_x'])>0) and np.all(np.diff(s['cortex_x'])>0)
        for j, xx in enumerate(x):
            data_rows.append(dict(snapshot=s['label'], time_after_onset_s=s['time'], reference_position_um=float(xx), **{key: float(f[key][j]) for key in ['rac', 'rho', 'mca_bound', 'mca_unbound', 'tension', 'displacement', 'velocity']}, native_membrane_velocity_um_per_s=float(s['membrane_v'][j]),
                local_10s_slip_um=float(s['slip'][j]), **{key+'_velocity_um_per_s': float(value[j]) for key, value in s['components'].items()}))
        for j, xx in enumerate(s['face_x']):
            flux_rows.append(dict(snapshot=s['label'], time_after_onset_s=s['time'], reference_face_position_um=float(xx), bound_advective_flux=float(s['advective'][j]), bound_diffusive_flux=float(s['diffusive'][j]), unbound_advective_flux=float(s['unbound_advective'][j]), unbound_diffusive_flux=float(s['unbound_diffusive'][j]), total_mca_spatial_flux=float(s['total_flux'][j])))
        for j, xx in enumerate(nodes):
            node_rows.append(dict(snapshot=s['label'], reference_position_um=float(xx), native_membrane_displacement_um=float(u[j]), local_10s_slip_um=float(slip[j]), plotted_membrane_position_um=float(s['membrane_x'][j]), plotted_cortex_marker_position_um=float(s['cortex_x'][j]), exaggeration=EXAGGERATION))
    csvwrite(stem+'_profiles.csv', data_rows)
    csvwrite(stem+'_fluxes.csv', flux_rows)
    csvwrite(stem+'_drawing_nodes.csv', node_rows)

    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False, 'pdf.fonttype':42, 'svg.fonttype':'none'})
    fig = plt.figure(figsize=(17.5, 21.5))
    grid = fig.add_gridspec(4, 3, left=.105, right=.975, bottom=.060, top=.850, width_ratios=[1.7, 1, 1], hspace=.33, wspace=.30)
    drawing_positions = []
    for s in snapshots:
        node_velocity = np.interp(nodes, np.r_[0,x,LENGTH], np.r_[s['boundary'][0],s['f']['velocity'],s['boundary'][1]])
        drawing_positions.extend([s['membrane_x'],s['cortex_x'],s['cortex_x']+40*node_velocity,
                                  s['cortex_x'][[0,-1]]+40*s['boundary']])
    low, high = float(np.concatenate(drawing_positions).min())-1, float(np.concatenate(drawing_positions).max())+1
    velocity_values = np.concatenate([np.r_[s['f']['velocity'], *s['components'].values()] for s in snapshots])
    vlo, vhi = float(velocity_values.min()), float(velocity_values.max())
    vpad = .12*max(vhi-vlo, .01)
    flux_values = np.concatenate([np.r_[s['advective'], s['diffusive'], s['unbound_advective'], s['unbound_diffusive'],s['total_flux']] for s in snapshots])*1000
    flo, fhi = float(flux_values.min()), float(flux_values.max())
    fpad = .12*max(fhi-flo, .01)
    diffusion_limit = max(float(max(abs(s['diffusive']).max(),abs(s['unbound_diffusive']).max())) for s in snapshots)*1.18e6
    max_tension = max(5., max(float(abs(s['f']['tension']).max()) for s in snapshots))
    tnorm = TwoSlopeNorm(vmin=-max_tension, vcenter=0, vmax=max_tension)
    bnorm = Normalize(0, 12)
    tmap, bmap = plt.get_cmap('coolwarm'), plt.get_cmap('YlGn')
    mechanical_axes = []
    for row, s in enumerate(snapshots):
        inner = grid[row, 0].subgridspec(5, 1, height_ratios=[3.1, 1, 1, 1, 1], hspace=.10)
        mech = fig.add_subplot(inner[0])
        title(mech, 'ADGJ'[row], s['label'])
        mech.set(xlim=(low, high), ylim=(-.58, 2.0), xticks=[], yticks=[1.7, 1.15, .05], yticklabels=['Reference', 'Membrane', 'Cortex'])
        mech.tick_params(axis='y', length=0, labelsize=9)
        mech.spines[['left', 'bottom']].set_visible(False)
        mechanical_axes.append(mech)
        profiles = [('rac','Rac',COLORS['rac']), ('rho','Rho',COLORS['rho']), ('mca_bound','Bound MCA',COLORS['mca']), ('tension','Tension',COLORS['tension'])]
        for j, (key, label, color) in enumerate(profiles):
            ax = fig.add_subplot(inner[j+1])
            values = s['f'][key]
            ax.plot(x, values, color=color, lw=1.65)
            ax.fill_between(x, 0, values, color=color, alpha=.09)
            if key=='rac':
                ax.set_ylim(0, max(1., float(values.max())*1.08))
                ax.axhline(.899061648345644, color='#888888', lw=.6, ls=':')
            elif key=='rho': ax.set_ylim(0, 1)
            elif key=='mca_bound': ax.set_ylim(0, .11)
            else: ax.set_ylim(min(-1, -max_tension*.04), max_tension*1.08)
            ax.yaxis.set_major_locator(MaxNLocator(2, min_n_ticks=2))
            ax.set_xlim(0, LENGTH)
            ax.set_xticks([0, 10, 20, 30])
            ax.tick_params(labelsize=8, length=2)
            ax.text(-.10, .5, label+('\n(pN)' if key=='tension' else '\n(conc.)'), transform=ax.transAxes, ha='right', va='center', color=color, fontsize=9)
            ax.grid(alpha=.12, axis='y')
            if j<3: ax.tick_params(axis='x', labelbottom=False)
            else: ax.set_xlabel('Reference position (µm)', fontsize=9, labelpad=1)

        ax = fig.add_subplot(grid[row, 1])
        title(ax, 'BEHK'[row], 'Cortical-flow components')
        for key, label in [('endpoint','Endpoint forcing'), ('rho','Rho stress gradient'), ('membrane','Membrane motion'), ('total','Total cortical flow')]:
            values = s['f']['velocity'] if key=='total' else s['components'][key]
            ax.plot(x, values, color=COLORS[key], ls='--' if key=='total' else '-', lw=2 if key=='total' else 1.7, label=label)
        ax.axhline(0, color='#aaaaaa', lw=.8)
        ax.grid(alpha=.13)
        ax.set(xlim=(0,LENGTH), ylim=(vlo-vpad, vhi+vpad), xlabel='Reference position (µm)', ylabel='Cortical velocity (µm/s)', xticks=[0,10,20,30])
        if row==0: ax.legend(frameon=False, fontsize=9, loc='lower right')
        ax.text(.04, .96, 'Positive = rightward', transform=ax.transAxes, fontsize=9, va='top', color='#555555')

        ax = fig.add_subplot(grid[row, 2])
        title(ax, 'CFIL'[row], 'MCA flux components')
        ax.plot(s['face_x'], 1000*s['advective'], color=COLORS['endpoint'], lw=2, label='Bound: advection (cortex)')
        ax.plot(s['face_x'], 1000*s['diffusive'], color=COLORS['diffusion'], lw=1.8, label='Bound: diffusion')
        ax.plot(s['face_x'], 1000*s['unbound_advective'], color=COLORS['unbound_advection'], lw=1.8, label='Unbound: advection (membrane)')
        ax.plot(s['face_x'], 1000*s['unbound_diffusive'], color=COLORS['unbound_diffusion'], lw=1.8, label='Unbound: diffusion')
        ax.plot(s['face_x'], 1000*s['total_flux'], color=COLORS['total'], lw=1.6, ls='--', label='Total MCA spatial flux')
        ax.axhline(0, color='#aaaaaa', lw=.8)
        ax.grid(alpha=.13)
        ax.set(xlim=(0,LENGTH), ylim=(flo-fpad, fhi+fpad), xlabel='Reference position (µm)', ylabel='MCA flux (10⁻³ conc. µm/s)', xticks=[0,10,20,30])
        ax.text(.04,.96,'Positive = rightward',transform=ax.transAxes,fontsize=9,va='top',color='#555555')
        if row==0: ax.legend(frameon=False, fontsize=8, loc='upper right', bbox_to_anchor=(1., .92))
        inset = ax.inset_axes([.50, .12, .45, .28])
        inset.plot(s['face_x'], s['diffusive']*1e6, color=COLORS['diffusion'], lw=1.4)
        inset.plot(s['face_x'], s['unbound_diffusive']*1e6, color=COLORS['unbound_diffusion'], lw=1.4)
        inset.axhline(0, color='#cccccc', lw=.6)
        inset.set(xlim=(0,LENGTH), ylim=(-diffusion_limit,diffusion_limit), xticks=[0,31], title='Diffusion (10⁻⁶ conc. µm/s)')
        inset.title.set_fontsize(8)
        inset.tick_params(labelsize=7, length=2)
        inset.patch.set_facecolor('#fafafa')
        inset.spines[['top','right']].set_visible(True)
        for spine in inset.spines.values(): spine.set_color('#bbbbbb')

    # Fix axes transformations before drawing screen-normal dashpot barrels.
    fig.canvas.draw()
    for ax, s in zip(mechanical_axes, snapshots):
        xm, xc, f = s['membrane_x'], s['cortex_x'], s['f']
        node_b = np.interp(nodes, x, f['mca_bound'])
        tensions = np.interp(.5*(nodes[:-1]+nodes[1:]), x, f['tension'])
        ax.plot(nodes, np.full(len(nodes),1.7),color='#cccccc',lw=.65,zorder=0)
        ax.scatter(nodes,np.full(len(nodes),1.7),color='#bbbbbb',s=9,zorder=1)
        for j in range(0,len(nodes),4): ax.plot([nodes[j],xm[j]],[1.66,1.27],color='#bbbbbb',ls=':',lw=.65,zorder=0)
        # Pale ribbons give the layers visible extent without adding mechanics.
        for positions, y, color in [(xm,1.15,'#D99C8A'),(xc,.05,COLORS['cortex'])]:
            ribbon=FancyBboxPatch((positions[0]-.17,y-.13),positions[-1]-positions[0]+.34,.26,boxstyle='round,pad=0.02,rounding_size=0.11',facecolor=color,edgecolor=color,alpha=.10,lw=.7,zorder=0)
            ax.add_patch(ribbon)
        for j in range(len(nodes)-1):
            draw_spring(ax,xm[j],xm[j+1],1.15,tmap(tnorm(tensions[j])))
            draw_dashpot(ax,[xc[j],.05],[xc[j+1],.05],COLORS['cortex'],width=1.15,halfwidth_pt=2.2)
        for j,value in enumerate(node_b):
            strength=p['linker_friction']*value
            draw_dashpot(ax,[xc[j],.14],[xm[j],1.06],bmap(bnorm(strength)),width=.7+2.2*value/.1,halfwidth_pt=2.5+3.0*value/.1)
        ax.scatter(xm,np.full(len(nodes),1.15),s=19,facecolor='white',edgecolor='#666666',lw=.7,zorder=5)
        ax.scatter(xc,np.full(len(nodes),.05),s=15,color=COLORS['cortex'],marker='s',zorder=5)
        v = np.interp(nodes,np.r_[0,x,LENGTH],np.r_[s['boundary'][0],f['velocity'],s['boundary'][1]])
        ax.quiver(xc,np.full(len(nodes),-.25),40*v,np.zeros(len(nodes)),angles='xy',scale_units='xy',scale=1,color=COLORS['cortex'],width=.0037,headwidth=3.5)
        ax.quiver(xc[[0,-1]],np.full(2,.05),40*s['boundary'],np.zeros(2),angles='xy',scale_units='xy',scale=1,color=COLORS['endpoint'],width=.0055,zorder=6)
        ax.annotate('',xy=(xm[0],1.48),xytext=(0,1.48),arrowprops=dict(arrowstyle='<->',color='#666666',lw=.8))
        ax.text(.015,.02,'Membrane displacement and 10-s relative slip ×2.5',transform=ax.transAxes,fontsize=8,color='#555555')

    fig.suptitle('Figure Y | Mechanics and MCA transport through cue withdrawal',x=.08,ha='left',y=.984,fontsize=21)
    fig.text(.08,.960,f'L = {LENGTH} µm · '+clock_description,fontsize=11,color='#555555')
    fig.text(.08,.941,'Left: deformed mechanics and native profiles. Middle: cortical-flow components. Right: signed bound + unbound MCA fluxes.',fontsize=11,color='#555555')
    handles=[Line2D([],[],color=COLORS['cortex'],lw=2,label='Cortex viscosity / flow'),Line2D([],[],color=COLORS['endpoint'],marker='>',ls='',label='Rac-controlled endpoint command'),Line2D([],[],color=COLORS['mca'],lw=3,label='MCA linker width and color ∝ local drag')]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.075,.927),ncol=3,frameon=False,fontsize=9,columnspacing=2)
    for rect,norm,cmap,label,ticks in [([.30,.896,.22,.007],tnorm,tmap,'Spring color: membrane tension (pN)',[-max_tension,0,max_tension]),([.68,.896,.22,.007],bnorm,bmap,'MCA linker drag (pN s / µm²)',[0,6,12])]:
        cax=fig.add_axes(rect)
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal',ticks=ticks)
        cb.ax.tick_params(labelsize=8,length=2)
        cb.set_label(label,fontsize=9,labelpad=1)
    fig.text(.08,.034,'All rows: geometry ×2.5; 10-s slip window; flow arrows: 4 µm = 0.1 µm/s. Profile positions and values are unexaggerated.',fontsize=10,color='#555555')
    fig.text(.08,.018,'Generated '+generated+' · source simulations complete · current report unchanged',fontsize=8,color='#777777')

    caption=('*Figure Y* Mechanical forcing sets cortical flow and bound/unbound MCA transport through cue withdrawal. '+clock_description+' All rows follow one shortest tested canonical 1%-strength-passing pulse from the same corrected independent-endpoint L = 31 µm checkpoint, with constant physical slope 0.24. '
        'A) Basal steady-state membrane/cortex illustration with Rac, Rho, bound MCA and tension profiles. B) Basal cortical-flow components. C) Basal bound/unbound MCA fluxes. '
        'D–F) The same three views 30 s after cue onset. G–I) The inherited state at cue withdrawal, before the first cue-free step; fields are continuous and none is reset. J–L) The same views 60 s after Rac at the incumbent left front drops below the activation-switch midpoint. '
        f'The downward threshold crossing is bracketed at {native_rt[k]:.1f}–{native_rt[k+1]:.1f} s after withdrawal on a native 0.5-s replay and linearly interpolated to {crossing:.6f} s; the final state is integrated to {target:.6f} s after withdrawal, exactly 60 s after that interpolated event. The threshold is {threshold:.12f}; it marks the midpoint of a smooth switch, not an abrupt binary loss of motion. '
        'In A, D, G and J, gray dots mark undeformed reference nodes. Membrane nodes are placed at reference position plus 2.5 times native membrane displacement; springs deform with these node positions and carry signed native tension on a common color scale. '
        'Cortex markers are registered to the membrane 10 s earlier and advanced by the integral of local cortical velocity over that window; their relative offset equals 2.5 times integrated cortical velocity minus the membrane displacement change. '
        'At basal steady state, the window is its stationary continuation: membrane velocity is zero and slip is 10 s times the native cortical velocity. The model has no accumulated cortex-displacement field; these are local kinematic markers, not permanent material particles or elastic MCA bond extensions. '
        'MCA dashpot barrel width, line weight and color encode the local bound-MCA drag density. Purple axial dashpots represent cortical viscosity. Purple arrows show native cortical flow; orange endpoint arrows show the prescribed command. '
        'The 17 drawing nodes resample the 201-cell solver; they are not a new coarse simulation. Underlying profiles use unexaggerated reference coordinates; concentration axes are labeled, with row-specific Rac limits and shared Rho/MCA/tension limits. The dotted Rac line is the original protrusion-switch midpoint. '
        'B, E, H and K separate endpoint forcing, Rho stress gradient and membrane-motion drag, retaining the native staggered forcing and actual final-step duration. Their sum reconstructs the saved cortical velocity. '
        'C, F, I and L show instantaneous face fluxes for both MCA pools: bound advection on cortical velocity, bound diffusion, unbound advection on membrane velocity, unbound diffusion, and the total of all four components. Positive velocity and flux mean rightward. '
        'Advective fluxes use the native upwind donor rule and their own endpoint velocities; membrane velocities are displacement increments divided by the actual preceding step duration. Diffusive flux is negative diffusivity times each pool’s spatial gradient, with homogeneous-Neumann boundary conditions. Bound diffusivity is 0.001 and unbound diffusivity 0.3 µm²/s. Insets show both diffusive contributions on a common labeled 10^-6 scale; main flux axes have a common 10^-3 scale. '
        'Local attachment exchange and uniform restoration are source terms, not spatial fluxes, and are not included in the plotted sum. Nonzero steady-state flux can balance these sources; it does not imply a changing steady profile. '
        'These are instantaneous diagnostic fluxes, not the prior split-step MCA polarity-rate budget or independently evolved pathway ablations. The 80-s dense local replay matches every retained state in its window, then a final short native step reaches the event-relative snapshot. No long integration or change to the main report. Scalar reversal strength and field-tail stability remain separate: this selected trial passes the canonical scalar gate and has field_tail_pass=False.')
    output(stem+'_caption.txt').write_text(caption+'\n')
    artifacts={}
    for ext in ['png','pdf','svg']:
        path=output(stem+'.'+ext)
        fig.savefig(path,dpi=200,facecolor='white')
        artifacts[layout.dev_artifact_key(OUT,path)]=sha(path)
    plt.close(fig)
    for name in [stem+'_caption.txt',stem+'_profiles.csv',stem+'_fluxes.csv',stem+'_drawing_nodes.csv',stem+'_event_trace.csv']:
        path=output(name);artifacts[layout.dev_artifact_key(OUT,path)]=sha(path)
    assert all(sha(REPORT/name)==digest for name,digest in report_hashes.items())
    audit=dict(passed=True,snapshots=[s['audit'] for s in snapshots],recovery_dense_replays=replay_audits,event=event,geometry_exaggeration=EXAGGERATION,slip_window_s=SLIP_WINDOW,checkpoint_hash_verified=True)
    record=dict(artifact_root='../..',generated_at=generated,producer='finished',legend=caption,sha256=artifacts,source_sha256=SOURCES,parameters=p,audit=audit,current_report_changed=False,new_long_integrations=0,
                selected_cue_duration_s=job['duration'],scalar_strength_pass=job['expected']['full_polarity_success'],field_tail_pass=job['expected']['field_tail_pass'],
                reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python -B '+str(Path(__file__))])
    output(stem+'.manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    print(OUT/(stem+'.png'))
    print(json.dumps(audit,indent=2))


if __name__=='__main__':
    main()
