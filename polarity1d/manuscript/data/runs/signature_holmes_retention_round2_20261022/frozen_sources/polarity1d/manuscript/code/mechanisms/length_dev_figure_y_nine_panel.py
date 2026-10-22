"""Nine-panel Figure Y: deformed mechanics, cortical flow and bound-MCA flux.

Default rows are basal steady state and two retained cue-on snapshots. The
withdrawal option uses short dense replays within retained recovery histories.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
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
              mca='#009E73', cortex='#762A83', rac='#C13C25', tension='#9D245D', diffusion='#009E73')
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


def snapshot(fields, times, i, x, parameters, label, when, basal=False):
    f = {key: value[i].copy() for key, value in fields.items()}
    dx = LENGTH/len(x)
    friction = parameters['linker_friction']*f['mca_bound']
    if basal:
        force_rho = f['rho']
        boundary = commands(f, parameters)
        membrane_v = np.zeros_like(x)
        slip = SLIP_WINDOW*f['velocity']
        endpoint_slip = SLIP_WINDOW*boundary
    else:
        assert abs(times[i]-times[i-1]-.5) < 1e-9
        j = int(np.argmin(abs(times-(times[i]-SLIP_WINDOW))))
        assert j > 0 and abs(times[i]-times[j]-SLIP_WINDOW) < 1e-9
        force_rho = fields['rho'][i-1]
        boundary = commands({key: value[i-1] for key, value in fields.items()}, parameters)
        membrane_v = (f['displacement']-fields['displacement'][i-1])/.5
        slip = np.trapezoid(fields['velocity'][j:i+1], times[j:i+1], axis=0)-(f['displacement']-fields['displacement'][j])
        boundary_history = commands({key: value[j-1:i] for key, value in fields.items()}, parameters)
        endpoint_slip = np.trapezoid(boundary_history, times[j:i+1], axis=0)-np.array([
            f['left_boundary_displacement'][0]-fields['left_boundary_displacement'][j, 0],
            f['right_boundary_displacement'][0]-fields['right_boundary_displacement'][j, 0]])
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
    divergence = -np.diff(advective+diffusive)/dx
    conservation_residual = float(abs(divergence.sum()*dx-((advective+diffusive)[0]-(advective+diffusive)[-1])))
    assert conservation_residual < 1e-12
    return dict(label=label, time=when, f=f, boundary=boundary, slip=slip, endpoint_slip=endpoint_slip,
        components=components, face_x=np.linspace(0, LENGTH, len(x)+1), advective=advective, diffusive=diffusive,
        audit=dict(label=label, velocity_reconstruction_residual_um_per_s=residual, native_tension_identity_residual_pN=tension_error, flux_conservation_residual=conservation_residual,
                   max_abs_advective_flux=float(abs(advective).max()), max_abs_diffusive_flux=float(abs(diffusive).max())))


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


def main(phase):
    generated = datetime.now(timezone.utc).isoformat()
    stem = 'figure_Y_cortical_motion_mca_9panel'+('_withdrawal' if phase == 'withdrawal' else '')
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
    replay_audits = []
    if phase == 'onset':
        source = register(ROOT/'data/runs/debelly_corrected_advection_20261008/L31_native.npz')
        summary = read(source.with_name('L31_summary.json'))
        assert sha(source)==summary['native_sha256'] and sha(checkpoint)==summary['checkpoint_sha256']
        with np.load(source) as z:
            t = z['time'].copy()
            fields = {key: z['field_'+key].copy() for key in FIELDS}
        for key in FIELDS:
            np.testing.assert_array_equal(fields[key][0], initial[key][0])
        for when, label in [(30., 'Early cue: 30 s after onset'), (200., 'Late cue: 200 s after onset')]:
            i = int(np.argmin(abs(t-when)))
            assert t[i] == when
            snapshots.append(snapshot(fields, t, i, x, p, label, when))
        clock_description = 'Basal steady state → 30 s and 200 s after cue onset; the cue remains on.'
    else:
        sys.path[:0] = [str(ROOT/'code'), str(ROOT/'code/core')]
        import debelly_corrected_reversal as corrected
        n = corrected.n
        for path in [Path(corrected.__file__), Path(n.__file__), Path(n.ds.__file__)]:
            register(path)
        assert not n.model_for('debelly').parameters.center_boundary_motion
        run = ROOT/'data/runs/length_narrative_memory_20261009'
        selection = read(run/'manifest.json')
        job = next(j for j in selection['jobs'] if j['case']['length']==LENGTH and j['case']['condition']=='fixed_slope')
        source = register(ROOT/job['source'] if job['source'] else run/f"workers/{selection['missing_jobs'].index(job):03d}.npz")
        assert sha(source)==read(source.with_suffix('.json'))['archive_sha256']
        with np.load(source) as z:
            t = z['recovery_time'].copy()
            fields = {key: z['recovery_field_'+key].copy() for key in FIELDS}
        for when, label in [(30., 'Early recovery: 30 s after withdrawal'), (500., 'Late recovery: 500 s after withdrawal')]:
            j = int(np.argmin(abs(t-(when-12.5))))
            assert abs(t[j]-(when-12.5)) < 1e-8
            state = {key: value[j].copy() for key, value in fields.items()}
            replay = n.run('debelly', LENGTH, state, n.ConstantInput(0), 12.5, interval=.5)
            local = {key: replay.fields[key].copy() for key in FIELDS}
            i = int(np.argmin(abs(t-when)))
            restart_error = max(float(np.max(abs(local[key][-1]-fields[key][i]))) for key in FIELDS)
            assert restart_error < 1e-9, restart_error
            snapshots.append(snapshot(local, replay.time, len(replay.time)-1, x, p, label, when))
            replay_audits.append(dict(time_after_withdrawal_s=when, replay_duration_s=12.5, maximum_native_endpoint_replay_error=restart_error))
        clock_description = f'Basal steady state → 30 s and 500 s after withdrawal of the {job["duration"]:.6f}-s shortest tested strength-passing cue.'

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
            data_rows.append(dict(snapshot=s['label'], time_s=s['time'], reference_position_um=float(xx), **{key: float(f[key][j]) for key in ['rac', 'rho', 'mca_bound', 'tension', 'displacement', 'velocity']},
                local_10s_slip_um=float(s['slip'][j]), **{key+'_velocity_um_per_s': float(value[j]) for key, value in s['components'].items()}))
        for j, xx in enumerate(s['face_x']):
            flux_rows.append(dict(snapshot=s['label'], time_s=s['time'], reference_face_position_um=float(xx), advective_flux=float(s['advective'][j]), diffusive_flux=float(s['diffusive'][j]), net_spatial_flux=float(s['advective'][j]+s['diffusive'][j])))
        for j, xx in enumerate(nodes):
            node_rows.append(dict(snapshot=s['label'], reference_position_um=float(xx), native_membrane_displacement_um=float(u[j]), local_10s_slip_um=float(slip[j]), plotted_membrane_position_um=float(s['membrane_x'][j]), plotted_cortex_marker_position_um=float(s['cortex_x'][j]), exaggeration=EXAGGERATION))
    csvwrite(stem+'_profiles.csv', data_rows)
    csvwrite(stem+'_fluxes.csv', flux_rows)
    csvwrite(stem+'_drawing_nodes.csv', node_rows)

    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False, 'pdf.fonttype':42, 'svg.fonttype':'none'})
    fig = plt.figure(figsize=(17.5, 16.5))
    grid = fig.add_gridspec(3, 3, left=.105, right=.975, bottom=.075, top=.815, width_ratios=[1.7, 1, 1], hspace=.33, wspace=.30)
    low = min(float(min(s['membrane_x'].min(), s['cortex_x'].min())) for s in snapshots)-2
    high = max(float(max(s['membrane_x'].max(), s['cortex_x'].max())) for s in snapshots)+2
    velocity_values = np.concatenate([np.r_[s['f']['velocity'], *s['components'].values()] for s in snapshots])
    vlo, vhi = float(velocity_values.min()), float(velocity_values.max())
    vpad = .12*max(vhi-vlo, .01)
    flux_values = np.concatenate([np.r_[s['advective'], s['diffusive'], s['advective']+s['diffusive']] for s in snapshots])*1000
    flo, fhi = float(flux_values.min()), float(flux_values.max())
    fpad = .12*max(fhi-flo, .01)
    diffusion_limit = max(float(abs(s['diffusive']).max()) for s in snapshots)*1.18e6
    max_tension = max(5., max(float(abs(s['f']['tension']).max()) for s in snapshots))
    tnorm = TwoSlopeNorm(vmin=-max_tension, vcenter=0, vmax=max_tension)
    bnorm = Normalize(0, 12)
    tmap, bmap = plt.get_cmap('coolwarm'), plt.get_cmap('YlGn')
    mechanical_axes = []
    for row, s in enumerate(snapshots):
        inner = grid[row, 0].subgridspec(5, 1, height_ratios=[3.1, 1, 1, 1, 1], hspace=.10)
        mech = fig.add_subplot(inner[0])
        title(mech, 'ADG'[row], s['label'])
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
        title(ax, 'BEH'[row], 'Cortical-flow components')
        for key, label in [('endpoint','Endpoint forcing'), ('rho','Rho stress gradient'), ('membrane','Membrane motion'), ('total','Total cortical flow')]:
            values = s['f']['velocity'] if key=='total' else s['components'][key]
            ax.plot(x, values, color=COLORS[key], ls='--' if key=='total' else '-', lw=2 if key=='total' else 1.7, label=label)
        ax.axhline(0, color='#aaaaaa', lw=.8)
        ax.grid(alpha=.13)
        ax.set(xlim=(0,LENGTH), ylim=(vlo-vpad, vhi+vpad), xlabel='Reference position (µm)', ylabel='Cortical velocity (µm/s)', xticks=[0,10,20,30])
        if row==0: ax.legend(frameon=False, fontsize=9, loc='lower right')
        ax.text(.04, .96, 'Positive = rightward', transform=ax.transAxes, fontsize=9, va='top', color='#555555')

        ax = fig.add_subplot(grid[row, 2])
        title(ax, 'CFI'[row], 'Bound-MCA flux components')
        ax.plot(s['face_x'], 1000*s['advective'], color=COLORS['endpoint'], lw=2, label='Net advective flux')
        ax.plot(s['face_x'], 1000*s['diffusive'], color=COLORS['diffusion'], lw=1.8, label='Diffusive flux')
        ax.plot(s['face_x'], 1000*(s['advective']+s['diffusive']), color=COLORS['total'], lw=1.6, ls='--', label='Sum of spatial fluxes')
        ax.axhline(0, color='#aaaaaa', lw=.8)
        ax.grid(alpha=.13)
        ax.set(xlim=(0,LENGTH), ylim=(flo-fpad, fhi+fpad), xlabel='Reference position (µm)', ylabel='MCA flux (10⁻³ conc. µm/s)', xticks=[0,10,20,30])
        ax.text(.04,.96,'Positive = rightward',transform=ax.transAxes,fontsize=9,va='top',color='#555555')
        if row==0: ax.legend(frameon=False, fontsize=9, loc='upper right', bbox_to_anchor=(1., .68))
        inset = ax.inset_axes([.50, .12, .45, .28])
        inset.plot(s['face_x'], s['diffusive']*1e6, color=COLORS['diffusion'], lw=1.4)
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

    fig.suptitle('Figure Y | Mechanics and MCA transport before and during reversal',x=.08,ha='left',y=.978,fontsize=21)
    fig.text(.08,.948,f'L = {LENGTH} µm · constant-slope cue 0.24 · '+clock_description,fontsize=11,color='#555555')
    fig.text(.08,.924,'Left: deformed spring/dashpot system and native profiles. Middle: cortical-flow components. Right: signed bound-MCA fluxes.',fontsize=11,color='#555555')
    handles=[Line2D([],[],color=COLORS['cortex'],lw=2,label='Cortex viscosity / flow'),Line2D([],[],color=COLORS['endpoint'],marker='>',ls='',label='Rac-controlled endpoint command'),Line2D([],[],color=COLORS['mca'],lw=3,label='MCA linker width and color ∝ local drag')]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.075,.910),ncol=3,frameon=False,fontsize=9,columnspacing=2)
    for rect,norm,cmap,label,ticks in [([.30,.871,.22,.009],tnorm,tmap,'Spring color: membrane tension (pN)',[-max_tension,0,max_tension]),([.68,.871,.22,.009],bnorm,bmap,'MCA linker drag (pN s / µm²)',[0,6,12])]:
        cax=fig.add_axes(rect)
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal',ticks=ticks)
        cb.ax.tick_params(labelsize=8,length=2)
        cb.set_label(label,fontsize=9,labelpad=1)
    fig.text(.08,.034,'All rows: geometry ×2.5; 10-s slip window; flow arrows: 4 µm = 0.1 µm/s. Profile positions and values are unexaggerated.',fontsize=10,color='#555555')
    fig.text(.08,.018,'Generated '+generated+' · source simulations complete · current report unchanged',fontsize=8,color='#777777')

    caption=('*Figure Y* Mechanical forcing sets cortical flow and bound-MCA transport before and during reversal. '+clock_description+' All rows use the same corrected independent-endpoint L = 31 µm checkpoint and constant-slope cue 0.24. '
        'A) Basal steady-state membrane/cortex illustration with Rac, Rho, bound MCA and tension profiles. B) Basal cortical-flow components. C) Basal bound-MCA fluxes. '
        'D–F) The same three views at the early time. G–I) The same views at the late time. '
        'In A, D and G, gray dots mark undeformed reference nodes. Membrane nodes are placed at reference position plus 2.5 times native membrane displacement; springs deform with these node positions and carry signed native tension on a common color scale. '
        'Cortex markers are registered to the membrane 10 s earlier and advanced by the integral of local cortical velocity over that window; their relative offset equals 2.5 times integrated cortical velocity minus the membrane displacement change. '
        'At basal steady state, the window is its stationary continuation: membrane velocity is zero and slip is 10 s times the native cortical velocity. The model has no accumulated cortex-displacement field; these are local kinematic markers, not permanent material particles or elastic MCA bond extensions. '
        'MCA dashpot barrel width, line weight and color encode the local bound-MCA drag density. Purple axial dashpots represent cortical viscosity. Purple arrows show native cortical flow; orange endpoint arrows show the prescribed command. '
        'The 17 drawing nodes resample the 201-cell solver; they are not a new coarse simulation. Underlying profiles use unexaggerated reference coordinates; concentration axes are labeled, with row-specific Rac limits and shared Rho/MCA/tension limits. The dotted Rac line is the original protrusion-switch midpoint. '
        'B, E and H separate endpoint forcing, Rho stress gradient and membrane-motion drag, retaining the native staggered forcing at nonzero cue times. Their sum reconstructs the saved cortical velocity. '
        'C, F and I show instantaneous face fluxes from the displayed bound-MCA and cortical-velocity fields: upwind advective flux, negative diffusivity times the bound-MCA spatial gradient, and their sum. Positive velocity and flux mean rightward. '
        'Diffusion uses homogeneous-Neumann boundary conditions; advective boundary flux uses the native commanded endpoint velocity and adjacent bound-MCA concentration. Diffusive-flux insets have common labeled 10^-6 units; main flux axes have common 10^-3 units. '
        'Local attachment exchange and uniform restoration are source terms, not spatial fluxes, and are not included in the plotted sum. Nonzero steady-state flux can balance these sources; it does not imply a changing steady profile. '
        'These are instantaneous diagnostic fluxes, not the prior split-step MCA polarity-rate budget or independently evolved pathway ablations. No change to the main report.')
    output(stem+'_caption.txt').write_text(caption+'\n')
    artifacts={}
    for ext in ['png','pdf','svg']:
        path=output(stem+'.'+ext)
        fig.savefig(path,dpi=200,facecolor='white')
        artifacts[layout.dev_artifact_key(OUT,path)]=sha(path)
    plt.close(fig)
    for name in [stem+'_caption.txt',stem+'_profiles.csv',stem+'_fluxes.csv',stem+'_drawing_nodes.csv']:
        path=output(name);artifacts[layout.dev_artifact_key(OUT,path)]=sha(path)
    assert all(sha(REPORT/name)==digest for name,digest in report_hashes.items())
    audit=dict(passed=True,phase=phase,snapshots=[s['audit'] for s in snapshots],recovery_dense_replays=replay_audits,geometry_exaggeration=EXAGGERATION,slip_window_s=SLIP_WINDOW,checkpoint_hash_verified=True)
    record=dict(artifact_root='../..',generated_at=generated,producer='finished',legend=caption,sha256=artifacts,source_sha256=SOURCES,parameters=p,audit=audit,current_report_changed=False,new_long_integrations=0,
                reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python -B '+str(Path(__file__))+' --phase '+phase])
    output(stem+'.manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    print(OUT/(stem+'.png'))
    print(json.dumps(audit,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--phase',choices=['onset','withdrawal'],default='onset')
    main(parser.parse_args().phase)
