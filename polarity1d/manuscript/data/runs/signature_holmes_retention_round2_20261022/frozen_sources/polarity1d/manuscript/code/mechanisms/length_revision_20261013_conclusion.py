"""Explain the long-cell exception using retained native bound-MCA budgets.

No trajectories are integrated. The pre-replacement bound pool is recovered
from the stored update and its recorded uniform correction. A conservative
face product rule, matching revised Figure 9, splits advection exactly into
profile-motion and velocity-gradient terms. Direct observed differences check
every native-step and window budget independently.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.linalg import solve_banded

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'figures/length_mechanisms_20261008'
OUT = REPORT / 'revision_20261013/conclusion'
DEV = REPORT / 'dev/mca_speed_decomposition'
RUN = ROOT / 'data/runs/debelly_mca_speed_decomposition_20261010'
STEM = 'Figure_12_long_cell_reversal'
LENGTHS = [18, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39]
EXAMPLES = [31, 37, 39]
SOURCES = {}
LENGTH_COLOR = {31: '#762A83', 37: '#009E73', 39: '#D55E00'}
TERMS = ['velocity_gradient', 'profile_motion', 'binding_exchange',
         'uniform_replacement', 'bound_diffusion']
TERM_LABEL = {'velocity_gradient': 'Flow gradient', 'profile_motion': 'Profile advection',
              'binding_exchange': 'Net binding',
              'uniform_replacement': 'Uniform replacement', 'bound_diffusion': 'Bound diffusion'}
TERM_COLOR = {'velocity_gradient': '#D55E00', 'profile_motion': '#0072B2',
              'binding_exchange': '#009E73', 'uniform_replacement': '#CC79A7',
              'bound_diffusion': '#7F7F7F'}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def source(path):
    path = Path(path)
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    return path


def read(path):
    return json.loads(source(path).read_text())


def independent_bound_solve(bold, uold, vold, edge, p, dt, dx):
    """Direct banded FV check from original unbound-pool samples."""
    n = len(bold)
    velocity_faces = .5 * (vold[:-1] + vold[1:])
    positive = dt / dx * np.maximum(velocity_faces, 0.)
    negative = dt / dx * np.minimum(velocity_faces, 0.)
    diffusion = dt * p['mca_bound_diffusion'] / dx**2
    matrix = np.zeros((3, n))
    matrix[1] = 1. + dt * p['mca_off'] + 2. * diffusion
    matrix[1, [0, -1]] -= diffusion
    matrix[0, 1:] = -diffusion + negative
    matrix[2, :-1] = -diffusion - positive
    matrix[1, :-1] += positive
    matrix[1, 1:] -= negative
    matrix[1, 0] -= dt * edge[0] / dx
    matrix[1, -1] += dt * edge[1] / dx
    return solve_banded((1, 1), matrix, bold + dt * p['mca_on'] * uold)


def collect():
    parent = read(DEV / 'supporting/metadata/MCA_speed_decomposition_across_length.manifest.json')
    assert parent['passed'] and parent['producer'] == 'finished'
    # Preserve complete upstream provenance, including independent checkpoint,
    # selection, dense-replay, and physical exchange/transport validations.
    for name, digest in parent['source_sha256'].items():
        assert sha(source(ROOT / name)) == digest, name
    for name, digest in parent['sha256'].items():
        assert sha(source(DEV / name)) == digest, name
    parent_audit = read(DEV / 'supporting/metadata/MCA_speed_decomposition_across_length_audit.json')
    assert parent_audit['passed'] and parent_audit['intervals'] == 8064
    p = parent_audit['parameters']
    assert not p['center_boundary_motion']
    interval = pd.read_csv(source(DEV / 'supporting/data/MCA_speed_decomposition_across_length_native_interval_budgets.csv'))
    trajectories, windows, native_rows, persistence, checks, binding_checks = {}, [], [], [], [], []
    max_error = dict(advection_product_rule=0.,
                     independent_recorded_advection=0., physical_bound_diffusion=0.,
                     native_front_budget=0., integrated_window=0., outflux_reconstruction=0.,
                     independent_original_unbound_preclosure=0., independent_physical_binding=0.)
    for slot, L in enumerate(LENGTHS):
        rr = interval[interval.length_um == L].reset_index(drop=True)
        meta = read(RUN / 'workers' / f'{slot:03d}.json')
        path = source(RUN / 'workers' / f'{slot:03d}.npz')
        assert meta['passed'] and sha(path) == meta['npz_sha256']
        with np.load(path) as z:
            t, x, b, v, rac, rho, tension = [z[f'L{L}_{key}'].copy()
                for key in ['time', 'x', 'bound', 'velocity', 'rac', 'rho', 'tension']]
        assert len(t) == len(rr) + 1
        np.testing.assert_allclose(t[:-1], rr.time_start_s, atol=1e-12, rtol=0)
        np.testing.assert_allclose(t[1:], rr.time_end_s, atol=1e-12, rtol=0)
        dt = np.diff(t)
        dx = L / len(x)
        front = x >= .8 * L
        assert front.sum() == 40
        uniform = rr.bound_uniform_source_conc_s.to_numpy()
        bpre = b[1:] - dt[:, None] * uniform[:, None]
        vf = np.column_stack([rr.left_command_um_s,
            .5 * (v[:-1, :-1] + v[:-1, 1:]), rr.right_command_um_s])
        donor = np.column_stack([bpre[:, 0],
            np.where(vf[:, 1:-1] >= 0, bpre[:, :-1], bpre[:, 1:]), bpre[:, -1]])
        gradient = -.5 * (donor[:, 1:] + donor[:, :-1]) * np.diff(vf, axis=1) / dx
        profile = -.5 * (vf[:, 1:] + vf[:, :-1]) * np.diff(donor, axis=1) / dx
        adv = -np.diff(vf * donor, axis=1) / dx
        diff_flux = np.column_stack([np.zeros(len(dt)),
            -p['mca_bound_diffusion'] * np.diff(bpre, axis=1) / dx, np.zeros(len(dt))])
        diffusion = -np.diff(diff_flux, axis=1) / dx
        actual = np.diff(b, axis=0) / dt[:, None]
        local_rate = pd.DataFrame(dict(
            velocity_gradient=gradient[:, front].mean(axis=1),
            profile_motion=profile[:, front].mean(axis=1),
            binding_exchange=-rr.front_depletion_exchange_conc_s,
            uniform_replacement=uniform,
            bound_diffusion=diffusion[:, front].mean(axis=1),
            observed=actual[:, front].mean(axis=1)))
        independent_adv = -rr[[f'front_depletion_{k}_conc_s' for k in
            ['interior_advection', 'old_end_exit', 'new_end_exit']]].sum(axis=1).to_numpy()
        max_error['advection_product_rule'] = max(max_error['advection_product_rule'],
            float(np.max(abs(gradient + profile - adv))))
        max_error['independent_recorded_advection'] = max(max_error['independent_recorded_advection'],
            float(np.max(abs(adv[:, front].mean(axis=1) - independent_adv))))
        max_error['physical_bound_diffusion'] = max(max_error['physical_bound_diffusion'],
            float(np.max(abs(local_rate.bound_diffusion + rr.front_depletion_diffusion_conc_s))))
        max_error['native_front_budget'] = max(max_error['native_front_budget'],
            float(np.max(abs(local_rate[TERMS].sum(axis=1) - local_rate.observed))))
        left, right = -vf[:, 0] * bpre[:, 0], vf[:, -1] * bpre[:, -1]
        max_error['outflux_reconstruction'] = max(max_error['outflux_reconstruction'],
            float(np.max(abs(left - rr.left_bound_outflux))),
            float(np.max(abs(right - rr.right_bound_outflux))))
        assert np.all(left >= 0) and np.all(right >= 0)
        bfront = b[:, front].mean(axis=1)
        minimum = int(np.argmin(bfront))
        # Report all lengths rather than infer uniqueness from three examples.
        early_weights = np.maximum(0., np.minimum(t[1:], 180.) - np.maximum(t[:-1], 60.))
        loss_fraction = float(early_weights @ left / (early_weights @ (left + right)))
        persistent_low = np.flatnonzero(np.maximum.accumulate(left[::-1])[::-1] < 1e-6)
        low_time = float(t[persistent_low[0]]) if len(persistent_low) else None
        right_majority = np.mean(v > 0, axis=1)
        sustained_majority = np.flatnonzero(np.minimum.accumulate(right_majority[::-1])[::-1] >= .95)
        majority_time = float(t[sustained_majority[0]]) if len(sustained_majority) else None
        persistence.append(dict(length_um=L, cue_duration_s=float(t[-1]),
            old_end_loss_fraction_60_180=loss_fraction,
            first_sustained_left_bound_outflux_below_1e_6=low_time,
            first_sustained_95pct_rightward_cortex_s=majority_time,
            front_bound_MCA_initial=float(bfront[0]),
            front_bound_MCA_minimum=float(bfront[minimum]),
            front_bound_MCA_minimum_time_s=float(t[minimum]),
            front_bound_MCA_at_withdrawal=float(bfront[-1]),
            old_endpoint_Rac_initial=float(rac[0, 0]),
            old_endpoint_Rac_at_180_s=float(np.interp(180., t, rac[:, 0])),
            old_endpoint_Rac_at_withdrawal=float(rac[-1, 0]),
            fraction_rightward_cortex_at_60_s=float(np.interp(60., t, right_majority)),
            fraction_rightward_cortex_at_180_s=float(np.interp(180., t, right_majority)),
            native_field_tail_pass=bool(rr.native_field_tail_pass.iloc[0])))
        if L not in EXAMPLES:
            continue
        # Check the physical source definition against original U(t), rather
        # than infer it only from a closed concentration budget. These are
        # diagnostic half-second solves, not new trajectories.
        original_job = next(j for j in parent_audit['jobs'] if j['length_um'] == L)
        with np.load(source(ROOT / original_job['archive'])) as z:
            original_t = z['pulse_time'].copy()
            original_u = z['pulse_field_mca_unbound'].copy()
            original_b = z['pulse_field_mca_bound'].copy()
            original_v = z['pulse_field_velocity'].copy()
        for sample_s in [60., 120., 180.] + ([507., 600., 700.] if L == 39 else []):
            oi = int(np.argmin(abs(original_t - sample_s)))
            i = int(np.argmin(abs(t - sample_s)))
            assert abs(original_t[oi] - sample_s) < 1e-12
            assert abs(t[i] - sample_s) < 1e-12
            np.testing.assert_array_equal(original_b[oi], b[i])
            np.testing.assert_array_equal(original_v[oi], v[i])
            independent = independent_bound_solve(original_b[oi], original_u[oi],
                original_v[oi], [vf[i, 0], vf[i, -1]], p, dt[i], dx)
            pool_error = float(np.max(abs(independent - bpre[i])))
            physical_binding = p['mca_on'] * original_u[oi] - p['mca_off'] * independent
            binding_error = abs(float(physical_binding[front].mean()) - local_rate.binding_exchange.iloc[i])
            max_error['independent_original_unbound_preclosure'] = max(
                max_error['independent_original_unbound_preclosure'], pool_error)
            max_error['independent_physical_binding'] = max(
                max_error['independent_physical_binding'], binding_error)
            binding_checks.append(dict(length_um=L, time_start_s=sample_s,
                next_step_s=float(dt[i]), preclosure_bound_max_error=pool_error,
                original_front_unbound_mean=float(original_u[oi, front].mean()),
                front_binding_from_original_unbound=float(physical_binding[front].mean()),
                front_binding_retained_budget=float(local_rate.binding_exchange.iloc[i]),
                front_physical_binding_absolute_error=binding_error,
                original_bound_and_velocity_exact_match=True))
        for lo, hi, label in [(60., 180., 'common_60_180')] + (
                [(507., 700., 'L39_rebound_507_700')] if L == 39 else []):
            assert hi <= t[-1]
            weights = np.maximum(0., np.minimum(t[1:], hi) - np.maximum(t[:-1], lo))
            assert abs(weights.sum() - (hi - lo)) < 1e-10
            means = {k + '_conc_s': float(weights @ local_rate[k] / weights.sum())
                     for k in TERMS + ['observed']}
            direct = float((np.interp(hi, t, bfront) - np.interp(lo, t, bfront)) / (hi - lo))
            error = abs(direct - sum(means[k + '_conc_s'] for k in TERMS))
            max_error['integrated_window'] = max(max_error['integrated_window'], error)
            windows.append(dict(length_um=L, window=label, start_s=lo, end_s=hi,
                duration_s=hi-lo, front_MCA_start=float(np.interp(lo, t, bfront)),
                front_MCA_end=float(np.interp(hi, t, bfront)),
                direct_observed_rate_conc_s=direct, **means))
        trajectories[L] = dict(time=t, front=bfront, old_rac=rac[:, 0],
            rho=rho, tension=tension, bound=b, velocity=v, x=x,
            interval=rr, rate=local_rate, left_outflux=left, right_outflux=right)
        for i, row in rr.iterrows():
            native_rows.append(dict(length_um=L, time_start_s=float(t[i]),
                time_end_s=float(t[i+1]), step_s=float(dt[i]),
                left_bound_outflux=float(left[i]), right_bound_outflux=float(right[i]),
                **{k + '_conc_s': float(local_rate.iloc[i][k]) for k in TERMS + ['observed']}))
        checks.append(dict(length_um=L, intervals=len(dt), front_cells=int(front.sum()),
            cue_duration_s=float(t[-1]), native_field_tail_pass=bool(rr.native_field_tail_pass.iloc[0])))
    assert max(max_error.values()) < 1e-11, max_error
    assert next(p for p in persistence if p['length_um'] == 39)['front_bound_MCA_minimum_time_s'] == 507.
    audit = dict(passed=True, generated_at=datetime.now(timezone.utc).isoformat(), producer='finished',
        upstream_native_intervals_checked=8064, example_native_intervals_exported=len(native_rows),
        new_integrations=0, source_audits_passed=True, maximum_absolute_errors=max_error,
        selected_geometry='relative-20%', front_definition='outermost right 40 of 201 cells',
        velocity_gradient_split='negative arithmetic-mean face donor density times face-velocity difference / dx',
        profile_motion_split='negative arithmetic-mean face velocity times face-donor density difference / dx',
        conservative_product_rule_exact=True, upwind_discretization_included=True,
        causal_ablations_performed=False, examples=checks,
        independent_physical_binding_checks=binding_checks,
        physical_binding_definition='mca_on * original_old_unbound - mca_off * independently solved preclosure_bound',
        persistent_left_loss_cutoff_conc_um_s=1e-6,
        persistent_rightward_flow_fraction_cutoff=.95,
        rebound_window='507 s, the recorded minimum of outer-20% bound MCA, through 700 s; descriptive interval, not a response threshold',
        common_comparison_window='60–180 s at every eligible length')
    return trajectories, windows, native_rows, persistence, audit


def render(tr, windows):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5,
        'axes.spines.top': False, 'axes.spines.right': False,
        'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white'})
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 8.0))
    fig.subplots_adjust(left=.085, right=.985, top=.955, bottom=.14, hspace=.42, wspace=.28)
    for ax, letter in zip(axes.flat, 'ABCD'):
        ax.text(-.13, 1.025, letter, transform=ax.transAxes, fontsize=15, weight='bold')
        ax.grid(axis='y', color='#DFE3E5', lw=.65, alpha=.7)
        ax.tick_params(labelsize=10)
    for L in EXAMPLES:
        q, c = tr[L], LENGTH_COLOR[L]
        ax = axes[0, 0]
        for key, ls in [('right_outflux', '-'), ('left_outflux', '--')]:
            ax.plot(q['time'][1:], q[key] * 1e3, color=c, lw=1.9, ls=ls)
        for ax, key in [(axes[0, 1], 'old_rac'), (axes[1, 0], 'front')]:
            ax.plot(q['time'], q[key], color=c, lw=2.05, label=f'{L} µm')
            ax.plot(q['time'][-1], q[key][-1], 'o', ms=4.8, mfc='white', mec=c, mew=1.3)
    axes[0, 0].set(ylabel='Bound-MCA loss at each end\n(10⁻³ conc. µm/s)', ylim=(-.15, 13.5))
    handles = [Line2D([], [], color=LENGTH_COLOR[L], lw=2, label=f'{L} µm') for L in EXAMPLES]
    legend = axes[0, 0].legend(handles=handles, ncol=3, loc='upper left', frameon=False,
        fontsize=9.5, columnspacing=1.2, handlelength=1.5)
    axes[0, 0].add_artist(legend)
    axes[0, 0].legend(handles=[Line2D([], [], color='#555555', ls=ls, lw=1.8, label=label)
        for ls, label in [('--', 'Old left end'), ('-', 'New right end')]],
        loc='upper left', bbox_to_anchor=(0., .89), ncol=2, frameon=False, fontsize=9.5)
    axes[0, 1].set(ylabel='Rac at the old left endpoint', ylim=(.69, .925))
    axes[1, 0].set(ylabel='Mean bound MCA in rightmost 20%', ylim=(.038, .106))
    axes[1, 0].axvspan(60, 180, color='#CAD6DC', alpha=.32, lw=0)
    axes[1, 0].axvspan(507, 700, color='#E7D2DC', alpha=.30, lw=0)
    axes[1, 0].plot(507, tr[39]['front'][np.argmin(tr[39]['front'])], 'o', ms=5,
        color=LENGTH_COLOR[39])
    axes[1, 0].text(512, .0414, '507 s', color=LENGTH_COLOR[39], fontsize=9.5, va='top')
    for ax in [axes[0, 0], axes[0, 1], axes[1, 0]]:
        ax.set(xlim=(0, 835), xticks=[0, 200, 400, 600, 800], xlabel='Time under relative-20% cue (s)')
    ax = axes[1, 1]
    xpos = np.array([0., 1., 2., 3.45])
    pos, neg = np.zeros(4), np.zeros(4)
    for term in TERMS:
        y = np.array([r[term + '_conc_s'] for r in windows]) * 1e4
        bottom = np.where(y >= 0, pos, neg)
        ax.bar(xpos, y, width=.62, bottom=bottom, color=TERM_COLOR[term],
            edgecolor='white', linewidth=.35, label=TERM_LABEL[term], zorder=3)
        pos += np.maximum(y, 0)
        neg += np.minimum(y, 0)
    net = np.array([r['direct_observed_rate_conc_s'] for r in windows]) * 1e4
    ax.plot(xpos, net, 'D', color='#202528', ms=5.1, label='Net change', zorder=5)
    ax.axhline(0, color='#565E62', lw=.85)
    ax.set(xticks=xpos, xticklabels=['31 µm', '37 µm', '39 µm', '39 µm'],
        ylim=(-8.1, 6.8), ylabel='Bound-MCA change in rightmost 20%\n(10⁻⁴ conc./s; positive = accumulation)')
    ax.text(1, -.16, 'Common 60–180 s', transform=ax.get_xaxis_transform(),
        ha='center', va='top', fontsize=9.6)
    ax.text(3.45, -.16, 'Rebound\n507–700 s', transform=ax.get_xaxis_transform(),
        ha='center', va='top', fontsize=9.6)
    ax.legend(ncol=2, loc='upper center', bbox_to_anchor=(.5, 1.28), frameon=False,
        fontsize=8.8, columnspacing=1.0, handlelength=1.3, handletextpad=.4)
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(OUT / (STEM + '.' + ext), dpi=220)
    plt.close(fig)


def narrative(windows, persistence):
    w = {(r['length_um'], r['window']): r for r in windows}
    pp = {r['length_um']: r for r in persistence}
    long, rebound = pp[39], w[(39, 'L39_rebound_507_700')]
    early = w[(39, 'common_60_180')]
    paragraphs = [
        f"The 39-µm cell retains strong MCA loss at the old left end for hundreds of seconds after the right cue begins (Figure 12A). Both ends lose MCA initially at shorter lengths too, but old-end loss soon becomes negligible: it falls below 10^-6 concentration·µm/s by {pp[31]['first_sustained_left_bound_outflux_below_1e_6']:.0f} s at 31 µm and {pp[37]['first_sustained_left_bound_outflux_below_1e_6']:g} s at 37 µm, compared with {long['first_sustained_left_bound_outflux_below_1e_6']:.0f} s at 39 µm. During the same 60–180-s interval, the old end accounts for {100*long['old_end_loss_fraction_60_180']:.1f}% of bound-MCA boundary loss at 39 µm and only {100*pp[37]['old_end_loss_fraction_60_180']:.1f}% at 37 µm. The distinction is the persistence and dominance of old-end loss, rather than the mere presence of two outward fluxes.",
        f"The old Rac pole also persists at 39 µm (Figure 12B). At 180 s its endpoint Rac is {long['old_endpoint_Rac_at_180_s']:.3f}, essentially its starting {long['old_endpoint_Rac_initial']:.3f}; the corresponding values at 31 and 37 µm are {pp[31]['old_endpoint_Rac_at_180_s']:.3f} and {pp[37]['old_endpoint_Rac_at_180_s']:.3f}. The native endpoint law converts this retained Rac, together with local tension, into the continuing leftward command. At 180 s the cortex still moves leftward in {100*(1-long['fraction_rightward_cortex_at_180_s']):.1f}% of the 39-µm domain, while it moves rightward throughout the 31- and 37-µm cells. This is a prolonged mechanical configuration during the selected pulse; these trajectories do not establish an additional bifurcation in length.",
        f"During 60–180 s, the rightmost 20% of the 39-µm cell loses bound MCA at {abs(early['observed_conc_s'])*1e4:.2f}, with rates here and below expressed in 10^-4 native concentration/s (Figure 12C–D). Flow-gradient dilution contributes {early['velocity_gradient_conc_s']*1e4:.2f} and advection of the existing profile contributes {early['profile_motion_conc_s']*1e4:.2f}; both remove MCA from that region. At 31 and 37 µm, profile advection instead supplies the region and partly offsets flow-gradient dilution. Local MCA clearance at 39 µm therefore coexists with a retained old Rac pole.",
        f"The right-front MCA reaches its minimum, {long['front_bound_MCA_minimum']:.4f}, at 507 s, then rises to {long['front_bound_MCA_at_withdrawal']:.4f} at cue withdrawal (Figure 12C). During 507–700 s, profile advection changes sign and supplies {rebound['profile_motion_conc_s']*1e4:.2f}, offsetting {100*rebound['profile_motion_conc_s']/abs(rebound['velocity_gradient_conc_s']):.1f}% of continued flow-gradient dilution. Binding and uniform replacement supply a further {(rebound['binding_exchange_conc_s']+rebound['uniform_replacement_conc_s'])*1e4:.2f}, leaving observed accumulation of {rebound['observed_conc_s']*1e4:.2f} (Figure 12D). Bound diffusion contributes only {rebound['bound_diffusion_conc_s']*1e4:.3f}. This local balance explains the rebound: advection of the existing MCA profile and replenishment together outweigh the diluting flow gradient."
    ]
    return paragraphs


CAPTION = """*Figure 12* The 39-µm cell retains old-end MCA loss, then reaccumulates MCA near the challenged end. A) Bound-MCA loss at the old left boundary (dashed) and challenged right boundary (solid). B) Rac in the first native cell at the old left endpoint. C) Mean bound MCA in the rightmost 40/201 cells. Open circles mark cue withdrawal; the filled point marks the 39-µm minimum at 507 s. Shading identifies the common 60–180-s comparison and descriptive 39-µm rebound interval, 507–700 s. D) Signed contributions to the same local bound-MCA change; positive bars add MCA, negative bars remove it, and diamonds show the measured net rate. The conservative flow-gradient/profile-advection split matches Figure 9, including endpoint fluxes and native upwind effects. All examples use exact corrected basal checkpoints and shortest tested scalar-1%-strength-passing relative-20% pulses (285.28, 355.99 and 809.08 s). After 10000-s recovery, 31 and 37 µm pass the separate native-field tail check; 39 µm fails it. All twelve lengths are checked for old-end-loss persistence. These within-trajectory budgets describe coupled dynamics; they are not independently evolved ablations or evidence of an additional length bifurcation."""


METHODS = """The renderer reuses all 8064 audited native intervals from debelly_mca_speed_decomposition_20261010. The recorded bound pool is the implicit transport/exchange solution plus its uniform correction; subtracting the recorded uniform increment therefore exactly recovers the pre-correction bound pool. Let u_j be native face velocity and d_j its upwind bound-MCA donor. The exact finite-volume identity splits -difference(u*d)/dx into -mean(d)*difference(u)/dx (flow gradient) and -mean(u)*difference(d)/dx (profile motion). At endpoints u is the current prescribed command and d is the adjacent pre-correction bound density; interior u averages the previous recorded cortical velocity. This is the same staggered, conservative split as Figure 9. It includes native numerical upwind effects, so its terms are discrete counterparts of the continuous product-rule terms, rather than a separately evolved physical intervention.

The right-front average uses the same rightmost 40 of 201 cells at every length. Each interval's local flow-gradient and profile-motion rates, independently calculated bound diffusion, retained physical exchange and uniform replacement are checked against the direct recorded change of this average. Window means use exact overlap lengths of the native intervals and are checked against direct endpoint concentration differences. All upstream source hashes and dense replay audits are verified. Physical binding is additionally checked against original old-unbound MCA at 60, 120 and 180 s in all three examples and at 507, 600 and 700 s in the 39-µm example. A separately assembled banded transport solve recovers the pre-correction bound pool; physical net binding is then mca_on * old_unbound - mca_off * pre-correction_bound. This check does not infer exchange merely from closure. The 507-s boundary is the recorded minimum of this particular local average, not a response gate. The 700-s endpoint exposes most of the rebound before late slowing and is declared explicitly. The sustained old-end loss cutoff of 10^-6 concentration·µm/s is a descriptive, exported threshold checked at all twelve lengths. At 39 µm the corresponding time is 568 s; this descriptive time is not identified with the cue threshold or a bifurcation.
"""


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source(Path(__file__))
    tr, windows, native, persistence, audit = collect()
    render(tr, windows)
    for suffix, rows in [('windows', windows), ('native_interval_budgets', native),
                         ('all_length_persistence', persistence)]:
        pd.DataFrame(rows).to_csv(OUT / f'{STEM}_{suffix}.csv', index=False)
    traces = []
    for L, q in tr.items():
        for i, time in enumerate(q['time']):
            traces.append(dict(length_um=L, time_s=float(time),
                old_endpoint_Rac=float(q['old_rac'][i]), front_bound_MCA=float(q['front'][i]),
                old_endpoint_bound_MCA=float(q['bound'][i, 0]),
                old_endpoint_Rho=float(q['rho'][i, 0]), old_endpoint_tension=float(q['tension'][i, 0]),
                rightward_cortical_cell_fraction=float(np.mean(q['velocity'][i] > 0))))
    pd.DataFrame(traces).to_csv(OUT / f'{STEM}_traces.csv', index=False)
    paragraphs = narrative(windows, persistence)
    (OUT / f'{STEM}_prose.md').write_text('\n\n'.join(paragraphs) + '\n')
    (OUT / f'{STEM}_caption.txt').write_text(CAPTION + '\n')
    (OUT / 'METHODS.md').write_text(METHODS)
    chapter = dict(id='long_cell_reversal', figure_label='12',
        title='Persistent old-end loss and later MCA rebound at 39 µm',
        image=f'revision_20261013/conclusion/{STEM}.png', paragraphs=paragraphs,
        caption=CAPTION, manifest=f'revision_20261013/conclusion/{STEM}.manifest.json', sources=[
            (f'revision_20261013/conclusion/{STEM}_windows.csv', 'Local bound-MCA window budgets'),
            (f'revision_20261013/conclusion/{STEM}_all_length_persistence.csv', 'Old-end-loss persistence at all twelve lengths'),
            (f'revision_20261013/conclusion/{STEM}_native_interval_budgets.csv', 'Native-step local contributions'),
            ('revision_20261013/conclusion/METHODS.md', 'Exact split and window definitions'),
            ('revision_20261013/editorial_review/physical_binding_independent_audit.json',
             'Independent physical-binding check against original unbound MCA')])
    (OUT / 'chapter.json').write_text(json.dumps(chapter, indent=2, ensure_ascii=False) + '\n')
    (OUT / f'{STEM}_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    artifacts = {p.name: sha(p) for p in OUT.iterdir()
        if p.is_file() and not p.name.endswith('.manifest.json')}
    manifest = dict(generated_at=audit['generated_at'], producer='finished', passed=True,
        renderer=str(Path(__file__).relative_to(ROOT)), renderer_sha256=sha(Path(__file__)),
        artifacts=artifacts, sources=SOURCES, audit=audit,
        reproduction='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 '
            '/home/4473331/.conda/envs/scaling_sandbox/bin/python -B '
            'polarity1d/manuscript/code/mechanisms/length_revision_20261013_conclusion.py')
    (OUT / f'{STEM}.manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(OUT / f'{STEM}.png')
    print(json.dumps(audit['maximum_absolute_errors'], indent=2))


if __name__ == '__main__':
    main()
