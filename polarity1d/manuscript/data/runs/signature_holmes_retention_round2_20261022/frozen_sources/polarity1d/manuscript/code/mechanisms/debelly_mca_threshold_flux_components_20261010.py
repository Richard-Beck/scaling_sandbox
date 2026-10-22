"""Figure Y B/C-style cortical flow and MCA fluxes at audited native endpoints.

One 0.5-s native replay per endpoint resolves membrane velocities, which are
not stored in compact checkpoints. Original checkpoints are never changed.
"""
from pathlib import Path
from dataclasses import asdict
import copy
import shutil

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import debelly_mca_threshold_profiles_20261010 as profiles
import length_dev_output_layout as layout
from length_dev_basal_front_tension_decomposition import solve

ROOT, OUT, FIG = profiles.ROOT, profiles.OUT, profiles.FIG
RESULTS = OUT / 'results/flux_components_cortical_flow_4lengths'
STEM = 'mca_threshold_flux_components'
LENGTHS = [19., 21., 25., 31.]
FLOW_COMPONENTS = [
    ('endpoint', 'Endpoint forcing', '#D55E00', '-'),
    ('rho', 'Rho stress gradient', '#0072B2', '-'),
    ('membrane', 'Membrane motion', '#929292', '-'),
    ('total', 'Total cortical flow', '#222222', '--'),
]
COMPONENTS = [
    ('bound_advection', 'Bound: advection (cortex)', '#D55E00', '-'),
    ('bound_diffusion', 'Bound: diffusion', '#009E73', '-'),
    ('unbound_advection', 'Unbound: advection (membrane)', '#CC79A7', ':'),
    ('unbound_diffusion', 'Unbound: diffusion', '#0072B2', '-'),
    ('total_flux', 'Total MCA spatial flux', '#222222', '--'),
]


def main():
    n, e, manifest = profiles.configure()
    assert profiles.read(OUT / 'results/audit.json')['all_profiles_verified']
    parameters = n.model_for('debelly').parameters
    assert not parameters.center_boundary_motion
    inputs = copy.deepcopy(manifest['inputs'])
    for job in manifest['jobs']:
        p = OUT / 'workers' / f'{job["slot"]:03d}.npz'
        r = profiles.read(p.with_suffix('.json'))
        inputs.append(dict(length=job['length'], file=str(p.relative_to(OUT)),
                           sha256=r['trace_sha256'], convergence=r['convergence']))
    RESULTS.mkdir(parents=True, exist_ok=True)
    shutil.copy2(__file__, RESULTS / 'source.py')
    sources = {str(Path(__file__).relative_to(ROOT)): n.sha(__file__),
               str(Path(profiles.__file__).relative_to(ROOT)): n.sha(profiles.__file__),
               str((OUT / 'manifest.json').relative_to(ROOT)): n.sha(OUT / 'manifest.json'),
               str((OUT / 'results/audit.json').relative_to(ROOT)): n.sha(OUT / 'results/audit.json')}
    for path in [Path(layout.__file__), Path(__file__).with_name('length_dev_basal_front_tension_decomposition.py')]:
        sources[str(path.relative_to(ROOT))] = n.sha(path)
    snapshots, rows, flow_rows, summaries, audits = {}, [], [], [], []
    for length in LENGTHS:
        inp = next(record for record in inputs if record['length'] == length)
        p = OUT / inp['file']
        assert inp['convergence']['verified'] and n.sha(p) == inp['sha256']
        sources[str(p.relative_to(ROOT))] = n.sha(p)
        before, _ = n.read_checkpoint(p, 'debelly')
        calls = []
        original_transport = n.ds._transport

        def capture(values, velocity, diffusion, step, dx, **kwargs):
            calls.append(dict(velocity=velocity.copy(), diffusion=diffusion,
                              step=step, dx=dx, **kwargs))
            return original_transport(values, velocity, diffusion, step, dx, **kwargs)

        n.ds._transport = capture
        try:
            trial = n.run('debelly', length, before, None, .5, interval=.5)
        finally:
            n.ds._transport = original_transport
        assert len(calls) == 2 and trial.time.tolist() == [0., .5]
        after = trial.state
        dx = length / len(after['rac'])
        membrane_v = (after['displacement'] - before['displacement']) / .5
        membrane_end_v = np.array([
            (after[k][0] - before[k][0]) / .5
            for k in ['left_boundary_displacement', 'right_boundary_displacement']])
        cortex_end_v = np.asarray(calls[0]['boundary_velocity'])
        np.testing.assert_allclose(membrane_v, calls[1]['velocity'], atol=1e-14, rtol=0)
        np.testing.assert_allclose(membrane_end_v, calls[1]['boundary_velocity'], atol=1e-14, rtol=0)
        np.testing.assert_allclose(before['velocity'], calls[0]['velocity'], atol=1e-14, rtol=0)
        assert calls[0]['diffusion'] == parameters.mca_bound_diffusion
        assert calls[1]['diffusion'] == parameters.mca_unbound_diffusion

        def flux(pool, velocity, boundary_v, diffusion):
            face_v = .5 * (velocity[:-1] + velocity[1:])
            donor = np.where(face_v >= 0, pool[:-1], pool[1:])
            advective = np.r_[boundary_v[0] * pool[0], face_v * donor, boundary_v[1] * pool[-1]]
            diffusive = np.r_[0., -diffusion * np.diff(pool) / dx, 0.]
            return advective, diffusive

        b_adv, b_diff = flux(after['mca_bound'], after['velocity'], cortex_end_v,
                             parameters.mca_bound_diffusion)
        u_adv, u_diff = flux(after['mca_unbound'], membrane_v, membrane_end_v,
                             parameters.mca_unbound_diffusion)
        values = dict(bound_advection=b_adv, bound_diffusion=b_diff,
                      unbound_advection=u_adv, unbound_diffusion=u_diff,
                      total_flux=b_adv + b_diff + u_adv + u_diff)
        assert all(np.isfinite(a).all() for a in values.values())
        # Independently reconstruct the post-step cortical force balance and
        # native membrane strain. These are diagnostic consistency checks.
        friction = parameters.linker_friction * after['mca_bound']
        active_gradient = parameters.active_stress * parameters.rho_to_rac_coupling * n.ds._gradient(before['rho'], dx)
        reconstructed_v = n.ds._dirichlet_elliptic(
            friction, parameters.cortical_viscosity,
            friction * membrane_v + active_gradient, dx, *cortex_end_v)
        v_error = float(np.max(abs(reconstructed_v - after['velocity'])))
        # Independently assembled Figure Y B superposition, with the exact
        # staggered forcing state and endpoint commands used by this step.
        zero = np.zeros_like(friction)
        flow = dict(
            endpoint=solve(friction, parameters.cortical_viscosity, zero, dx, *cortex_end_v),
            rho=solve(friction, parameters.cortical_viscosity,
                      parameters.active_stress * parameters.rho_to_rac_coupling
                      * np.gradient(before['rho'], dx), dx),
            membrane=solve(friction, parameters.cortical_viscosity, friction * membrane_v, dx))
        flow_error = float(np.max(abs(sum(flow.values()) - after['velocity'])))
        assert flow_error < 1e-11
        flow['total'] = after['velocity'].copy()
        assert all(np.isfinite(value).all() for value in flow.values())
        tension = parameters.membrane_stiffness * n.ds._gradient(
            after['displacement'], dx, after['left_boundary_displacement'][0],
            after['right_boundary_displacement'][0])
        tension_error = float(np.max(abs(tension - after['tension'])))
        assert v_error < 1e-11 and tension_error < 1e-9
        max_changes = {key: float(np.max(abs(after[key] - before[key])))
                       for key in ['rac', 'rho', 'mca_bound', 'mca_unbound', 'velocity', 'tension']}
        replay_path = RESULTS / f'L{length:g}_replay.npz'
        n.save_checkpoint(replay_path, 'debelly', after, trial)
        audit = dict(length_um=length, checkpoint_sha256=n.sha(p),
                     replay_sha256=n.sha(replay_path), replay_step_s=.5,
                     cortical_force_balance_max_error=v_error,
                     independent_flow_superposition_max_error_um_per_s=flow_error,
                     cortical_boundary_velocity_um_per_s=cortex_end_v.tolist(),
                     tension_identity_max_error=tension_error,
                     native_membrane_velocity_max_error=float(abs(membrane_v - calls[1]['velocity']).max()),
                     replay_max_abs_field_changes=max_changes)
        audits.append(audit)
        x = np.linspace(0, length, len(after['rac']) + 1)
        weights = 2 * (np.arange(len(after['rac'])) + .5) / len(after['rac']) - 1
        summary = dict(length_um=length,
                       abs_rac_moment=float(abs(after['rac'] @ weights / after['rac'].sum())),
                       bound_mca_range=float(np.ptp(after['mca_bound'])),
                       max_abs_membrane_velocity_um_per_s=float(abs(membrane_v).max()))
        summary.update({f'max_abs_{key}': float(abs(value).max()) for key, value in values.items()})
        summaries.append(summary)
        cell_x = (np.arange(len(after['rac'])) + .5) * dx
        summary.update({f'max_abs_{key}_cortical_velocity_um_per_s': float(abs(value).max())
                        for key, value in flow.items()})
        snapshots[length] = dict(x=x, cell_x=cell_x, flow=flow, **values, summary=summary)
        for i, position in enumerate(cell_x):
            flow_rows.append(dict(length_um=length, cell_index=i, x_nominal_um=float(position),
                                  x_fraction=float(position / length),
                                  membrane_velocity_um_per_s=float(membrane_v[i]),
                                  **{key + '_cortical_velocity_um_per_s': float(value[i])
                                     for key, value in flow.items()}))
        for i, position in enumerate(x):
            rows.append(dict(length_um=length, face_index=i, x_nominal_um=float(position),
                             x_fraction=float(position / length),
                             **{key: float(value[i]) for key, value in values.items()}))
    n.csvwrite(RESULTS / 'spatial_flux_components.csv', rows)
    n.csvwrite(RESULTS / 'cortical_flow_components.csv', flow_rows)
    n.csvwrite(RESULTS / 'length_summary.csv', summaries)
    n.write(RESULTS / 'audit.json', dict(generated_at=n.now(), passed=True, cases=audits,
                                      sources=sources, frozen_model_sources=manifest['sources']))

    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42,
                         'svg.fonttype': 'none'})
    fig = plt.figure(figsize=(18, 10.8))
    flow_grid = fig.add_gridspec(1, 4, left=.065, right=.985, bottom=.59,
                                top=.825, wspace=.28)
    flux_grid = fig.add_gridspec(1, 4, left=.065, right=.985, bottom=.13,
                                top=.44, wspace=.28)
    all_flow = np.concatenate([s['flow'][key] for s in snapshots.values()
                               for key, *_ in FLOW_COMPONENTS])
    all_flux = np.concatenate([s[key] for s in snapshots.values()
                               for key, *_ in COMPONENTS]) * 1000
    all_bound_diffusion = np.concatenate([s['bound_diffusion'] for s in snapshots.values()]) * 1e6
    def limits(values):
        lo, hi = min(float(values.min()), 0), max(float(values.max()), 0)
        pad = .1 * max(hi - lo, 1e-9)
        return lo - pad, hi + pad
    for index, length in enumerate(LENGTHS):
        s = snapshots[length]
        cortex = fig.add_subplot(flow_grid[0, index])
        for key, _, color, ls in FLOW_COMPONENTS:
            cortex.plot(s['cell_x'], s['flow'][key], color=color, ls=ls,
                        lw=2 if key == 'total' else 1.7)
        cortex.set(xlim=(0, length), ylim=limits(all_flow), xlabel='Nominal position (µm)')
        cortex.set_title(f'{chr(65 + index)}  {length:g} µm', loc='left',
                         weight='bold', fontsize=12, pad=9)
        grid = flux_grid[0, index].subgridspec(2, 1, height_ratios=[2.25, 1], hspace=.13)
        ax = fig.add_subplot(grid[0])
        detail = fig.add_subplot(grid[1], sharex=ax)
        for key, _, color, ls in COMPONENTS:
            ax.plot(s['x'], 1000 * s[key], color=color, ls=ls,
                    lw=1.7 if key != 'total_flux' else 1.5,
                    zorder=4 if key == 'unbound_advection' else 3)
        detail.plot(s['x'], 1e6 * s['bound_diffusion'], color='#009E73', lw=1.7)
        ax.set_title(f'{chr(69 + index)}  {length:g} µm', loc='left',
                     weight='bold', fontsize=11, pad=7)
        ax.set(xlim=(0, length), ylim=limits(all_flux))
        detail.set(xlim=(0, length), ylim=limits(all_bound_diffusion), xlabel='Nominal position (µm)')
        ax.tick_params(axis='x', labelbottom=False)
        for a in [cortex, ax, detail]:
            a.axhline(0, color='#BBBBBB', lw=.65, zorder=0)
            a.grid(alpha=.14)
            a.tick_params(labelsize=9)
        if index % 4 == 0:
            cortex.set_ylabel('Cortical velocity (µm/s)')
            ax.set_ylabel('MCA flux (10⁻³ conc. µm/s)')
            detail.set_ylabel('Bound diffusion\n(10⁻⁶ conc. µm/s)')
        detail.text(.97, .12, 'Bound diffusion detail', transform=detail.transAxes,
                    va='bottom', ha='right', fontsize=8, color='#555555')
    handles = [Line2D([0], [0], color=color, ls=ls, lw=2, label=label)
               for _, label, color, ls in COMPONENTS]
    flow_handles = [Line2D([0], [0], color=color, ls=ls, lw=2, label=label)
                    for _, label, color, ls in FLOW_COMPONENTS]
    fig.suptitle('Steady cortical flow and MCA flux components across cell lengths', fontsize=18, weight='bold', y=.97)
    fig.text(.5, .923, 'Native polarized endpoints · positive velocity and flux are rightward · shared vertical scales across lengths',
             ha='center', fontsize=11, color='#444444')
    fig.text(.065, .873, 'Cortical-flow components', fontsize=12, weight='bold')
    fig.legend(handles=flow_handles, ncol=4, loc='upper center', bbox_to_anchor=(.61, .895),
               frameon=False, fontsize=10, columnspacing=2.5)
    fig.text(.065, .505, 'MCA flux components', fontsize=12, weight='bold')
    fig.legend(handles=handles, ncol=3, loc='upper center', bbox_to_anchor=(.61, .531),
               frameon=False, fontsize=10, columnspacing=2.5)
    fig.text(.065, .022,
             'Verified native endpoints, no cue; each sampled after one additional 0.5-s step to resolve membrane motion.\n'
             'Top: endpoint, Rho-gradient and membrane contributions sum to cortical flow (Figure Y B definitions). This is instantaneous source superposition.\n'
             'Below: MCA fluxes (Figure Y C definitions) and bound diffusion magnified 1,000×. Local exchange and replenishment are outside the spatial flux sum.',
             fontsize=10, linespacing=1.5)
    artifacts = []
    for ext in ['png', 'pdf', 'svg']:
        path = layout.dev_output_path(FIG, STEM + '.' + ext)
        fig.savefig(path, dpi=220)
        artifacts.append(path)
    plt.close(fig)
    for filename, source in [(STEM + '.csv', RESULTS / 'spatial_flux_components.csv'),
                             (STEM + '_cortical_flow.csv', RESULTS / 'cortical_flow_components.csv'),
                             (STEM + '_length_summary.csv', RESULTS / 'length_summary.csv'),
                             (STEM + '_audit.json', RESULTS / 'audit.json')]:
        dest = layout.dev_output_path(FIG, filename)
        shutil.copy2(source, dest)
        artifacts.append(dest)
    caption = (
        'Steady native deBelly cortical-flow and MCA spatial flux components at nominal lengths '
        + ', '.join(f'{v:g}' for v in LENGTHS) + ' microns. '
        'Top panels A-D follow Figure Y B: endpoint forcing (orange), direct Rho stress '
        'gradient (blue), membrane motion (gray), and native total cortical flow (black dashed). '
        'These are exact instantaneous linear source superpositions at fixed bound-MCA '
        'friction, not causal fractions of the coupled dynamics. Each component is solved '
        'with an independently assembled cell-centred Dirichlet operator. Endpoint forcing '
        'uses the actual commands and zero interior forcing; Rho and membrane contributions '
        'use zero endpoint velocities. The forcing Rho is from the beginning of the native '
        'step, while friction and total flow are from its end, preserving the staggered solve. '
        'Their sum is checked against the native post-step velocity below 1e-11 micron/s. '
        'Flow is plotted at native cell centres in microns/s, with a common vertical scale '
        'across all four lengths. Lower panels E-H '
        'use exactly the instantaneous face-flux definitions of Figure Y C: bound advection '
        'on cortical velocity (orange), bound diffusion (green), unbound advection on membrane '
        'velocity (magenta dotted), unbound diffusion (blue), and their sum (black dashed). '
        'Advection uses the native upwind donor and native independent endpoint velocities; '
        'diffusion uses minus diffusivity times the face gradient with zero diffusive endpoint '
        'flux. Bound diffusivity is 0.001 and unbound diffusivity is 0.3 square microns/s. '
        'Positive is rightward on the nominal reference domain. The four main flux axes '
        'share a vertical scale in 10^-3 concentration*micron/s. Bound-diffusion details '
        'share a common scale across all lengths in 10^-6 concentration*micron/s '
        '(1000-fold magnification); unbound diffusion remains visible on the main axes. '
        'All starting states passed the original convergence and holdout rules and retain '
        'asymmetry. Unlike a transport-step '
        'source budget, these instantaneous diagnostics use the displayed final MCA pools '
        'and cortical velocity. The membrane velocity and cortical boundary commands are '
        'taken from a single native 0.5-s no-cue replay from each exact endpoint; fields are '
        'displayed at the end of that step. This supplies velocities missing in compact '
        'checkpoints without modifying their bytes. The replay and maximum field changes '
        'are retained in the audit. Local attachment/detachment and uniform replenishment '
        'are sources, not spatial flux components, and are excluded. Spatial flux need not '
        'vanish at a steady concentration profile. No mass-balance panel is included.\n')
    path = layout.dev_output_path(FIG, STEM + '_caption.txt')
    path.write_text(caption)
    artifacts.append(path)
    n.write(layout.dev_output_path(FIG, STEM + '.manifest.json'), dict(
        generated_at=n.now(), artifact_root='../..', lengths_um=LENGTHS,
        diagnostic_replay_s=.5, parameters=asdict(parameters), sources=sources,
        native_evidence=str(RESULTS.relative_to(ROOT)), audit_passed=True,
        sha256={str(p.relative_to(FIG)): n.sha(p) for p in artifacts}))
    print('Audited and exported', len(snapshots), 'lengths:', ', '.join(f'{v:g}' for v in LENGTHS))
    for row in summaries:
        print(row)


if __name__ == '__main__':
    main()
