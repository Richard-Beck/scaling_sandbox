"""Compact concluding native MCA-clearance figure for the revised report.

This renderer consumes audited retained budgets and trajectories.  It integrates
their exact native-step rates over a common window, checks the integrations
against direct concentration/centroid differences, and preserves source hashes.
No native integration, parameter change, or causal pathway removal is performed.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import ast
import __future__
from types import SimpleNamespace

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parents[1]
DEV = ROOT / 'figures/length_mechanisms_20261008/dev'
OUT = ROOT / 'figures/length_mechanisms_20261008/revision_20261010/conclusion'
RUN = ROOT / 'data/runs/debelly_mca_speed_decomposition_20261010'
STEM = 'Figure_13_MCA_clearance'
LENGTHS = [18, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39]
EXAMPLES = [23, 31, 39]
SOURCES = {}


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


def verify_parent(folder, stem):
    manifest = read(folder / 'supporting/metadata' / (stem + '.manifest.json'))
    assert manifest['passed'] and manifest['producer'] == 'finished'
    for name, digest in manifest['sha256'].items():
        path = source(folder / name)
        assert sha(path) == digest, path
    for name, digest in manifest['source_sha256'].items():
        path = ROOT / name
        assert sha(path) == digest, path
        SOURCES[name] = digest
    return read(folder / 'supporting/metadata' / (stem + '_audit.json'))


def csvwrite(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def collect():
    speed_audit = verify_parent(DEV / 'mca_speed_decomposition',
                                'MCA_speed_decomposition_across_length')
    flow_audit = verify_parent(DEV / 'mca_flow_budget',
                               'MCA_flow_redistribution_budget')
    assert speed_audit['passed'] and flow_audit['passed']
    assert not speed_audit['parameters']['center_boundary_motion']
    assert speed_audit['intervals'] == 8064
    parameters = speed_audit['parameters']
    native_path = source(ROOT / 'data/runs/length_narrative_memory_20261009/'
        'frozen_sources/code/core/models/debelly.py')
    native_ast = ast.parse(native_path.read_text())
    functions = [node for node in native_ast.body if isinstance(node, ast.FunctionDef)
        and node.name in ['smooth_threshold', '_feedback_terms']]
    assert len(functions) == 2
    namespace = {'np': np}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(native_path),
        'exec', flags=__future__.annotations.compiler_flag), namespace)
    native_feedback = namespace['_feedback_terms']
    native_parameters = SimpleNamespace(**parameters)
    intervals = pd.read_csv(source(DEV / 'mca_speed_decomposition/supporting/data/'
        'MCA_speed_decomposition_across_length_native_interval_budgets.csv'))
    late = pd.read_csv(source(DEV / 'mca_speed_decomposition/supporting/data/'
        'MCA_speed_decomposition_across_length_window_summary.csv'))
    flow = pd.read_csv(source(DEV / 'mca_flow_budget/supporting/data/'
        'MCA_flow_redistribution_budget_window_summary.csv'))
    trajectories, summaries, checks, traces, feedback_checks = {}, [], [], [], []
    direct_max = dict(centroid=0., local_front=0.,
                      native_step_centroid=0., native_step_local_front=0.,
                      mechanical_closure=0., late_centroid=0.)
    keys = ['centroid_speed_um_s', 'mean_cortical_velocity_um_s'] + [
        k for k in intervals.columns if k.startswith(('front_depletion_',
            'mechanical_speed_', 'speed_'))]
    for slot, L in enumerate(LENGTHS):
        rr = intervals[intervals.length_um == L]
        assert len(rr) and np.all(rr.step_s > 0)
        starts, ends = rr.time_start_s.to_numpy(), rr.time_end_s.to_numpy()
        weight = np.maximum(0., np.minimum(ends, 180.) - np.maximum(starts, 60.))
        assert abs(weight.sum() - 120.) < 1e-10
        means = {k: float(weight @ rr[k].to_numpy() / 120.) for k in keys}
        path = RUN / 'workers' / f'{slot:03d}.npz'
        meta = read(path.with_suffix('.json'))
        assert meta['passed'] and sha(source(path)) == meta['npz_sha256']
        with np.load(path) as z:
            t, c, x, b = [z[f'L{L}_{k}'].copy() for k in ['time', 'centroid', 'x', 'bound']]
            original_t = z[f'L{L}_original_time'].copy()
            original_c = z[f'L{L}_original_centroid'].copy()
            rac = z[f'L{L}_rac'].copy()
        front = x >= .8 * L
        oldfront = x <= .2 * L
        assert front.sum() == oldfront.sum() == 40
        assert np.max(abs(t[1:] - ends)) < 1e-10
        bfront, bold = b[:, front].mean(axis=1), b[:, oldfront].mean(axis=1)
        direct_centroid = (np.interp(60., t, c) - np.interp(180., t, c)) / 120.
        direct_front = (np.interp(60., t, bfront) - np.interp(180., t, bfront)) / 120.
        direct_old = (np.interp(180., t, bold) - np.interp(60., t, bold)) / 120.
        transport = sum(means['front_depletion_' + k + '_conc_s'] for k in
            ['interior_advection', 'new_end_exit', 'old_end_exit', 'diffusion'])
        refill = -sum(means['front_depletion_' + k + '_conc_s'] for k in
            ['exchange', 'uniform_replenishment'])
        assert transport > 0 and refill > 0 and transport > refill
        assert abs(direct_front - (transport - refill)) < 1e-12
        assert abs(direct_centroid - means['centroid_speed_um_s']) < 1e-12
        direct_max['centroid'] = max(direct_max['centroid'],
            abs(direct_centroid - means['centroid_speed_um_s']))
        direct_max['local_front'] = max(direct_max['local_front'],
            abs(direct_front - transport + refill))
        measured_step_centroid = -np.diff(c) / np.diff(t)
        measured_step_front = -np.diff(bfront) / np.diff(t)
        budget_step_front = rr[[k for k in rr if k.startswith('front_depletion_')]].sum(axis=1).to_numpy()
        direct_max['native_step_centroid'] = max(direct_max['native_step_centroid'],
            float(np.max(abs(measured_step_centroid - rr.centroid_speed_um_s))))
        direct_max['native_step_local_front'] = max(direct_max['native_step_local_front'],
            float(np.max(abs(measured_step_front - budget_step_front))))
        all_mechanical = sum(rr['mechanical_speed_' + k + '_um_s'].to_numpy() for k in
            ['left_drive', 'right_drive', 'rho', 'membrane', 'initial_remainder'])
        all_advection = sum(rr['speed_' + k + '_um_s'].to_numpy() for k in
            ['interior_advection', 'old_end_exit', 'new_end_exit'])
        direct_max['mechanical_closure'] = max(direct_max['mechanical_closure'],
            float(np.max(abs(all_mechanical - all_advection))))
        late_row = late[(late.length_um == L) & (late.window == 'late')].iloc[0]
        D = float(ends[-1])
        direct_late = (np.interp(D - 60., original_t, original_c) - original_c[-1]) / 60.
        direct_max['late_centroid'] = max(direct_max['late_centroid'],
            abs(direct_late - late_row.centroid_speed_um_s))
        row = dict(length_um=L, window_start_s=60., window_end_s=180.,
            cue_duration_s=D, native_field_tail_pass=bool(rr.native_field_tail_pass.iloc[0]),
            mean_cortical_velocity_um_s=means['mean_cortical_velocity_um_s'],
            centroid_speed_um_s=direct_centroid, net_front_transport_depletion_conc_s=transport,
            front_binding_refill_conc_s=-means['front_depletion_exchange_conc_s'],
            front_uniform_refill_conc_s=-means['front_depletion_uniform_replenishment_conc_s'],
            front_total_refill_conc_s=refill, observed_front_depletion_conc_s=direct_front,
            old_front_accumulation_conc_s=direct_old,
            front_at_start=float(bfront[0]), front_at_withdrawal=float(bfront[-1]),
            old_front_at_start=float(bold[0]), old_front_at_withdrawal=float(bold[-1]),
            late_rho_centroid_contribution_um_s=float(late_row.mechanical_speed_rho_um_s),
            late_endpoint_centroid_contribution_um_s=float(late_row.mechanical_speed_endpoint_um_s),
            late_observed_centroid_speed_um_s=float(late_row.centroid_speed_um_s))
        row['late_other_centroid_contribution_um_s'] = row['late_observed_centroid_speed_um_s'] - (
            row['late_rho_centroid_contribution_um_s'] + row['late_endpoint_centroid_contribution_um_s'])
        summaries.append(row)
        if 23 <= L <= 31:
            compare = flow[(flow.length_um == L) & (flow.window == 'later')].iloc[0]
            assert abs(compare.centroid_speed_um_s - direct_centroid) < 1e-12
            assert abs(compare.mean_cortical_velocity_um_s - means['mean_cortical_velocity_um_s']) < 1e-12
            assert abs(sum(compare[k] for k in compare.index if k.startswith('front_depletion_')) - direct_front) < 1e-12
        checks.append(dict(length_um=L, front_cells=int(front.sum()),
            interval_count=len(rr), common_window_duration_s=float(weight.sum()),
            direct_front_budget_error_conc_s=abs(direct_front - transport + refill),
            direct_centroid_budget_error_um_s=abs(direct_centroid - means['centroid_speed_um_s'])))
        if L in EXAMPLES:
            activation, _ = native_feedback(native_parameters, b, np.zeros_like(b))
            independent_activation = parameters['rac_mechanical_feedback'] * parameters['rho_to_rac_coupling'] * (
                1. - np.tanh((b - parameters['mca_rac_threshold']) / parameters['mca_switch_scale'])) / 2.
            assert np.max(abs(activation - independent_activation)) < 1e-15
            new_activation = activation[:, front].mean(axis=1)
            old_activation = activation[:, oldfront].mean(axis=1)
            assert new_activation[-1] > new_activation[0]
            assert old_activation[-1] < old_activation[0]
            native_job = next(job for job in speed_audit['jobs'] if job['length_um'] == L)
            archive = source(ROOT / native_job['archive'])
            with np.load(archive) as z:
                recovery_bound = z['recovery_field_mca_bound'][-1]
                recovery_rac = z['recovery_field_rac'][-1]
                assert z['recovery_time'][-1] == 10000.
            assert recovery_bound[front].mean() < recovery_bound[oldfront].mean()
            assert recovery_rac[-1] > recovery_rac[0]
            feedback_checks.append(dict(length_um=L,
                new_front_MCA_Rac_activation_initial=float(new_activation[0]),
                new_front_MCA_Rac_activation_withdrawal=float(new_activation[-1]),
                old_front_MCA_Rac_activation_initial=float(old_activation[0]),
                old_front_MCA_Rac_activation_withdrawal=float(old_activation[-1]),
                recovery_final_new_front_bound_MCA=float(recovery_bound[front].mean()),
                recovery_final_old_front_bound_MCA=float(recovery_bound[oldfront].mean()),
                recovery_final_new_endpoint_Rac=float(recovery_rac[-1]),
                recovery_final_old_endpoint_Rac=float(recovery_rac[0]),
                native_feedback_source_verified=True))
            trajectories[L] = dict(time=t, new_front=bfront, old_front=bold,
                new_rac=rac[:, -1], old_rac=rac[:, 0], duration=D)
            for j in range(len(t)):
                traces.append(dict(length_um=L, time_s=float(t[j]),
                    new_front_bound_MCA=float(bfront[j]), old_front_bound_MCA=float(bold[j]),
                    new_endpoint_Rac=float(rac[j, -1]), old_endpoint_Rac=float(rac[j, 0])))
    assert max(direct_max.values()) < 1e-10, direct_max
    anchor = next(r for r in summaries if r['length_um'] == 23)
    for row in summaries:
        for key in ['mean_cortical_velocity_um_s', 'centroid_speed_um_s',
                    'observed_front_depletion_conc_s']:
            row[key + '_relative_to_23um'] = row[key] / anchor[key]
    audit = dict(passed=True, native_intervals=8064, source_audits_passed=True,
        independent_maximum_errors=direct_max, common_window_checks=checks,
        reference_coordinates='nominal cell-centre coordinates',
        local_front_definition='rightmost 40/201 cells; width 40L/201',
        late_window='last 60 s of each length-specific selected cue; original recorded-interval interpolation',
        native_step_replays_in_source_campaign=[18, 37, 39],
        new_integrations=0, mechanical_terms_are_source_allocations=True,
        native_MCA_feedback_and_recovery_checks=feedback_checks,
        source_native_step_max_centroid_error_um_s=speed_audit['maximum_step_centroid_budget_error_um_s'])
    return summaries, trajectories, traces, audit


LENGTH_COLOR = {23: '#0072B2', 31: '#762a83', 39: '#D55E00'}
COLOR = dict(transport='#2B78A2', refill='#D59139', net='#24262A',
             flow='#168A80', centroid='#3B6FB6', clearance='#B56080',
             rho='#2869A4', endpoint='#D46D32', other='#92989B')


def render(rows, trajectories):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .8,
        'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white'})
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 8.7))
    fig.subplots_adjust(left=.09, right=.98, bottom=.095, top=.955,
                        hspace=.39, wspace=.29)
    for ax, letter in zip(axes.flat, 'ABCD'):
        ax.text(-.16, 1.025, letter, transform=ax.transAxes, fontsize=16, weight='bold')
        ax.grid(axis='y', color='#DEE2E3', lw=.65, alpha=.7, zorder=0)
        ax.tick_params(labelsize=10)
    ax = axes[0, 0]
    for L in EXAMPLES:
        tr = trajectories[L]
        for key, ls in [('new_front', '-'), ('old_front', '--')]:
            ax.plot(tr['time'], tr[key], color=LENGTH_COLOR[L], ls=ls, lw=2.2)
            ax.plot(tr['duration'], tr[key][-1], 'o', ms=5.5,
                    mec=LENGTH_COLOR[L], mfc='white', mew=1.5)
    ax.axvspan(60, 180, color='#CCD7DB', alpha=.3, zorder=0, lw=0)
    ax.set(xlim=(0, 850), ylim=(.04, .112), xticks=[0, 200, 400, 600, 800],
           yticks=[.04, .05, .06, .07, .08, .09, .10], xlabel='Time under cue (s)',
           ylabel='Mean bound MCA near each end')
    legend = ax.legend(handles=[Line2D([], [], color=LENGTH_COLOR[L], lw=2.2,
        label=f'{L} µm') for L in EXAMPLES], ncol=3, loc='upper right',
        frameon=False, fontsize=10, handlelength=1.5, columnspacing=1.3)
    ax.add_artist(legend)
    ax.legend(handles=[Line2D([], [], color='#555555', lw=2, ls=ls, label=label)
        for ls, label in [('-', 'New front'), ('--', 'Old front')]],
        frameon=False, loc='upper right', bbox_to_anchor=(1., .90),
        ncol=2, fontsize=10, handlelength=2., columnspacing=1.2)

    xs = np.array([r['length_um'] for r in rows])
    def vals(key):
        return np.array([r[key] for r in rows])
    ax = axes[0, 1]
    transport = vals('net_front_transport_depletion_conc_s') * 1e4
    refill = vals('front_total_refill_conc_s') * 1e4
    ax.fill_between(xs, refill, transport, color=COLOR['transport'], alpha=.14, lw=0)
    for yy, label, color, ls in [(transport, 'Transport depletion', COLOR['transport'], '-'),
            (refill, 'Binding + uniform refill', COLOR['refill'], '-'),
            (transport - refill, 'Observed depletion', COLOR['net'], '--')]:
        ax.plot(xs, yy, color=color, marker='o', ms=4, lw=2.2, ls=ls, label=label)
    ax.set(xlim=(17.5, 39.5), ylim=(0, 6.3), xticks=[18, 23, 27, 31, 35, 39],
        xlabel='Cell length (µm)', ylabel='New-front MCA rate (10⁻⁴ conc./s)')
    ax.legend(loc='upper left', frameon=False, fontsize=9.8, handlelength=2.1)

    ax = axes[1, 0]
    focused = [r for r in rows if 23 <= r['length_um'] <= 31]
    xf = [r['length_um'] for r in focused]
    for key, label, color, marker in [
            ('mean_cortical_velocity_um_s', 'Mean cortical flow', COLOR['flow'], 'o'),
            ('centroid_speed_um_s', 'MCA centroid speed', COLOR['centroid'], 's'),
            ('observed_front_depletion_conc_s', 'New-front MCA depletion', COLOR['clearance'], '^')]:
        ax.plot(xf, [r[key + '_relative_to_23um'] for r in focused], color=color,
                lw=2.3, marker=marker, ms=4.7, label=label)
    ax.plot(xf, np.array(xf) / 23., color='#9CA2A4', lw=1.4, ls=':', label='Cell length')
    ax.axhline(1, color='#8B9194', lw=.7, zorder=0)
    ax.set(xlim=(22.7, 32.15), ylim=(.75, 2.22), xticks=xf,
        xlabel='Cell length (µm)', ylabel='Change relative to 23 µm')
    ax.legend(loc='upper left', frameon=False, fontsize=9.8, handlelength=2.1)
    for key, color, dy in [('mean_cortical_velocity_um_s', COLOR['flow'], 0.),
            ('centroid_speed_um_s', COLOR['centroid'], .055),
            ('observed_front_depletion_conc_s', COLOR['clearance'], -.012)]:
        value = focused[-1][key + '_relative_to_23um']
        ax.text(31.2, value + dy, f'{value:.2f}×', color=color, fontsize=10.7,
                weight='bold', va='center')

    ax = axes[1, 1]
    for key, label, color, ls in [
            ('late_rho_centroid_contribution_um_s', 'Rho-driven advection', COLOR['rho'], '-'),
            ('late_endpoint_centroid_contribution_um_s', 'Endpoint-driven advection', COLOR['endpoint'], '-'),
            ('late_other_centroid_contribution_um_s', 'All remaining terms', COLOR['other'], ':'),
            ('late_observed_centroid_speed_um_s', 'Observed centroid speed', COLOR['net'], '--')]:
        ax.plot(xs, vals(key) * 1e3, color=color, lw=2.2 if 'observed' in key else 1.9,
                marker='o' if 'observed' in key else None, ms=4., ls=ls, label=label)
    ax.axhline(0, color='#8B9194', lw=.8, zorder=0)
    ax.set(xlim=(17.5, 39.5), ylim=(-23, 23), xticks=[18, 23, 27, 31, 35, 39],
        xlabel='Cell length (µm)', ylabel='Late centroid-speed contribution (nm/s)')
    ax.legend(loc='lower left', bbox_to_anchor=(.01, .025), frameon=False,
              fontsize=9.5, handlelength=2.1)
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(OUT / (STEM + '.' + ext), dpi=240)
    plt.close(fig)


def narrative(rows):
    byL = {r['length_um']: r for r in rows}
    a, b, long = byL[23], byL[31], byL[39]
    flowratio = b['mean_cortical_velocity_um_s_relative_to_23um']
    centroidratio = b['centroid_speed_um_s_relative_to_23um']
    front_ratio = b['observed_front_depletion_conc_s_relative_to_23um']
    text = f'''MCA clearance links mechanical redistribution to persistent Rac reversal. During successful pulses, bound MCA falls near the new right front and rises near the old left front (Figure 13A). The native feedback makes this consequential: lower bound MCA raises Rac activation, while higher MCA suppresses it. These changes favor Rac at the new front and weaken the old front. The attachment profiles remain partly reorganized at cue withdrawal; subsequent recovery establishes low MCA at the new front and high MCA at the old front.

The new front clears because transport removes more bound MCA than binding and uniform replacement return (Figure 13B). Most removal is cancelled by refill. Over the common 60–180-s window, transport depletion rises from {a['net_front_transport_depletion_conc_s']*1e4:.3f} to {b['net_front_transport_depletion_conc_s']*1e4:.3f} × 10^-4 concentration units/s between 23 and 31 µm. Refill rises further, from {a['front_total_refill_conc_s']*1e4:.3f} to {b['front_total_refill_conc_s']*1e4:.3f}, reducing the observed depletion from {a['observed_front_depletion_conc_s']*1e4:.3f} to {b['observed_front_depletion_conc_s']*1e4:.3f} on the same scale.

Faster cortex therefore does not guarantee faster MCA clearance (Figure 13C). From 23 to 31 µm, mean cortical flow increases {flowratio:.2f}-fold, centroid speed increases only {centroidratio:.2f}-fold, and actual new-front depletion becomes {(1-front_ratio)*100:.1f}% slower. Interior transport and boundary loss offset one another in the centroid budget. Uniform replacement can also move a right-biased centroid toward the midpoint while refilling the new front. Cortical speed, centroid speed and local clearance consequently measure different processes; the local budget exposes the small imbalance that clears MCA from the future front.

The longest cell shows a second limit (Figure 13D). Near withdrawal, Rho-associated advection supplies {long['late_rho_centroid_contribution_um_s']*1e3:.2f} nm/s toward the new rear, opposed by {abs(long['late_endpoint_centroid_contribution_um_s'])*1e3:.2f} nm/s from endpoint forcing. Net centroid movement falls to {long['late_observed_centroid_speed_um_s']*1e3:.3f} nm/s, and MCA returns to the new-front region (A). This cell clears that region rapidly early, then loses redistribution speed; its long cue does not reflect uniformly slow clearance. Larger cells thus face a larger redistribution burden without proportionate acceleration: stronger motion is offset by opposing mechanical drives, transport and continuous refill.
'''
    return text


CAPTION = '''*Figure 13* MCA clears the new front only through the small imbalance between transport and refill, so stronger cortical motion does not guarantee faster redistribution. A) Mean bound MCA in the new right-front region (solid) and old left-front region (dashed) during the shortest tested scalar-1%-strength-passing relative-20% cues at 23, 31 and 39 µm; open circles mark cue withdrawal. Each region contains the outermost 40 of 201 cells (width 40L/201). Gray shading marks the common 60–180-s window used in B–C. B) Local new-front transport depletion, binding plus uniform refill, and their observed difference over that common window at all twelve eligible lengths. Shading between transport and refill shows the small residual depletion. Transport includes interior and boundary advection and bound diffusion; refill combines physical binding exchange with the original uniform mass correction. C) Mean rightward cortical flow, leftward bound-MCA centroid speed, and actual local front depletion over 60–180 s, each normalized to its value at 23 µm; the dotted line gives relative cell length. The interval 23–31 µm is the matched range where rightward flow increases strongly while local clearance slows. D) Signed Rho- and endpoint-associated contributions to bound-MCA centroid speed during the final 60 s of each selected cue. Positive contributions move the centroid left, toward the new rear. All remaining terms comprise membrane-associated advection, diffusion, binding, uniform replacement and the negligible initial reconstruction remainder; the three contributions sum to the observed centroid speed. Mechanical curves are exact source allocations on the native coupled trajectory with common upwind donors, rather than independently evolved ablations. Late windows use length-specific cue durations. All lengths use the corrected independent-endpoint model and exact original basal checkpoints. Scalar success is assessed after 10000 s of recovery; native-field tail flags are separate and exported. Common-window and native-step identities are independently checked against recorded MCA profiles. No new trajectories or parameter changes.
'''


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source(Path(__file__))
    rows, trajectories, traces, audit = collect()
    render(rows, trajectories)
    csvwrite(OUT / (STEM + '_summary.csv'), rows)
    csvwrite(OUT / (STEM + '_traces.csv'), traces)
    (OUT / (STEM + '_caption.txt')).write_text(CAPTION)
    prose = narrative(rows)
    (OUT / (STEM + '_prose.md')).write_text(prose)
    paragraphs = [paragraph.strip() for paragraph in prose.strip().split('\n\n')]
    assert len(paragraphs) == 4
    chapter = dict(title='Continuous refill limits MCA clearance',
        image='revision_20261010/conclusion/' + STEM + '.png',
        paragraphs=paragraphs, caption=CAPTION.strip(), sources=[])
    (OUT / 'chapter.json').write_text(json.dumps(chapter, indent=2, ensure_ascii=False) + '\n')
    generated = datetime.now(timezone.utc).isoformat()
    audit.update(generated_at=generated, producer='finished')
    (OUT / (STEM + '_audit.json')).write_text(json.dumps(audit, indent=2) + '\n')
    artifacts = {p.name: sha(p) for p in OUT.glob(STEM + '*')
        if p.is_file() and not p.name.endswith('.manifest.json')}
    artifacts['chapter.json'] = sha(OUT / 'chapter.json')
    manifest = dict(generated_at=generated, producer='finished', passed=True,
        renderer=str(Path(__file__).relative_to(REPO)), renderer_sha256=sha(Path(__file__)),
        artifacts=artifacts, sources=SOURCES, audit=audit,
        selected_geometry='relative-20%', lengths_um=LENGTHS, new_integrations=0,
        reproduction='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 '
            '/home/4473331/.conda/envs/scaling_sandbox/bin/python -B '
            'polarity1d/manuscript/code/mechanisms/length_revision_20261010_conclusion.py')
    (OUT / (STEM + '.manifest.json')).write_text(json.dumps(manifest, indent=2) + '\n')
    print(OUT / (STEM + '.png'))
    print(json.dumps(audit['independent_maximum_errors'], indent=2))


if __name__ == '__main__':
    main()
