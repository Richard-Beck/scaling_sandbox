"""Development Figure Y: native mechanics, MCA advection, distance and speed.

Reuses corrected, audited archives. Three small elliptic solves reconstruct
one saved cortical velocity; no new trajectory or report rebuild is performed.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

import length_dev_output_layout as layout
from length_dev_basal_front_tension_decomposition import solve

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'figures/length_mechanisms_20261008'
OUT = REPORT / 'dev/resting_state'
BASE = ROOT / 'data/runs/debelly_corrected_three_test_20261008'
FLOW = ROOT / 'data/runs/debelly_corrected_advection_20261008'
QUALITY = ROOT / 'data/runs/length_quality_debelly_traces_20261009'
MEMORY = ROOT / 'data/runs/length_mechanisms_20261008/narrative_20261009/memory'
STEM = 'figure_Y_cortical_motion_mca'
COLORS = dict(endpoint='#D55E00', rho='#0072B2', membrane='#888888',
              nonadvection='#AA4499', total='#222222', mca='#009E73', cortex='#762a83')
SOURCES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def register(path):
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    return path


def read(path):
    return json.loads(register(path).read_text())


def rows(path):
    with register(path).open(newline='') as handle:
        return list(csv.DictReader(handle))


def output(name):
    return layout.dev_output_path(OUT, name)


def writecsv(name, records):
    with output(name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def dashpot(ax, start, end, color, width=1.1, size=.065):
    start, end = np.asarray(start), np.asarray(end)
    delta = end - start
    normal = np.array([-delta[1], delta[0]]) / np.linalg.norm(delta)
    a, b, piston = start + .32 * delta, start + .70 * delta, start + .50 * delta
    for line in [np.array([start, a]), np.array([a+size*normal, b+size*normal, b-size*normal, a-size*normal]), np.array([piston-size*normal, piston+size*normal]), np.array([piston, end])]:
        ax.plot(*line.T, color=color, lw=width, solid_capstyle='round')


def spring(ax, left, right, y, color):
    xx = np.r_[left, np.linspace(left+.12*(right-left), right-.12*(right-left), 13), right]
    yy = np.full_like(xx, y)
    yy[2:-2] += .045 * np.resize([1, -1], len(xx)-4)
    ax.plot(xx, yy, color=color, lw=1.7)


def panel(ax, letter, title):
    ax.set_title(title, loc='left', fontsize=12, pad=11)
    ax.text(-.055 if letter == 'A' else -.13, 1.045, letter, transform=ax.transAxes, fontsize=15, weight='bold')


def main():
    generated = datetime.now(timezone.utc).isoformat()
    report_hashes = {name: sha(REPORT / name) for name in ['index.html', 'length_mechanisms.pdf']}
    for path in [Path(__file__), Path(layout.__file__), Path(__file__).with_name('length_dev_basal_front_tension_decomposition.py')]:
        register(path)
    parameters = read(OUT / 'supporting/metadata/basal_endpoint_metrics_vs_length.manifest.json')['parameters']
    assert not parameters['center_boundary_motion']
    base = read(BASE / 'manifest.json')
    assert not base['center_boundary_motion']
    assert read(BASE / 'results/independent_audit.json')['passed']
    assert read(FLOW / 'independent_audit.json')['passed']
    quality = read(QUALITY / 'manifest.json')
    assert read(QUALITY / 'status.json')['independent_audit_passed']
    selection = read(ROOT / 'data/runs/length_narrative_memory_20261009/manifest.json')
    memory = [r for r in rows(MEMORY / 'memory_measurements.csv') if r['condition'] == 'fixed_slope']
    budgets = [r for r in rows(MEMORY / 'cue_timing_budget.csv') if r['panel'] == 'A']
    budget_meta = read(REPORT / 'narrative_20261009/memory/debelly_budget.manifest.json')
    assert budget_meta['passed']
    lengths = sorted({int(c['length']) for c in base['cases'] if c['model'] == 'debelly' and c['eligible']})
    memory.sort(key=lambda r: float(r['length']))
    budgets.sort(key=lambda r: float(r['length']))
    assert lengths == [int(r['length']) for r in memory] == [int(r['length']) for r in budgets]
    measured, checks = [], []
    for r, budget in zip(memory, budgets):
        length = int(r['length'])
        case = next(c for c in base['cases'] if c['model'] == 'debelly' and c['length'] == length and c['condition'] == 'fixed_slope')
        checkpoint = register(ROOT / case['checkpoint'])
        assert sha(checkpoint) == case['checkpoint_sha256'] == r['checkpoint_sha256']
        with np.load(checkpoint) as initial:
            b = initial['state_mca_bound'].copy()
            x = initial['x'].copy()
        initial_distance = float(b @ (x-length/2) / b.sum())
        slot = next(i for i, j in enumerate(quality['jobs']) if j['case']['length'] == length and j['case']['condition'] == 'fixed_slope')
        archive = register(QUALITY / f'workers/{slot:03d}.npz')
        assert sha(archive) == read(archive.with_suffix('.json'))['archive_sha256']
        with np.load(archive) as z:
            times = z['continuous_time'].copy()
            bound = z['continuous_field_mca_bound'].copy()
        np.testing.assert_allclose(bound[0], b, atol=1e-12)
        polarity = bound @ (2*x/length-1) / bound.sum(axis=1)
        i = int(np.flatnonzero(polarity <= 0)[0])
        assert i > 0
        zero_time = float(times[i-1] + (times[i]-times[i-1])*polarity[i-1]/(polarity[i-1]-polarity[i]))
        speed = initial_distance / zero_time
        np.testing.assert_allclose([initial_distance, zero_time, speed], [float(r['centroid_distance_to_zero']), float(r['mca_zero']), float(r['mean_centroid_speed_to_zero'])], rtol=1e-11, atol=1e-11)
        selected = next(j for j in selection['jobs'] if j['case']['length'] == length and j['case']['condition'] == 'fixed_slope')
        assert selected['duration'] > 60 and float(budget['actual_time']) == 60
        rates = {k: float(budget['rate_'+k]) for k in ['boundary', 'rho', 'membrane', 'diffusion', 'exchange', 'replenishment']}
        net = float(budget['reorientation_rate'])
        residual = abs(sum(rates.values())-net)
        assert residual < 1e-11 and float(budget['mechanical_residual']) < 1e-11
        measured.append(dict(length_um=length, initial_centroid_distance_um=initial_distance, uninterrupted_mca_zero_s=zero_time, mean_distribution_speed_um_per_s=speed,
            shortest_strength_passing_cue_s=selected['duration'], instantaneous_budget_time_s=60,
            endpoint_advection_rate_per_s=rates['boundary'], rho_advection_rate_per_s=rates['rho'], membrane_advection_rate_per_s=rates['membrane'],
            diffusion_rate_per_s=rates['diffusion'], exchange_rate_per_s=rates['exchange'], restoration_rate_per_s=rates['replenishment'],
            nonadvection_rate_per_s=sum(rates[k] for k in ['diffusion', 'exchange', 'replenishment']), net_rate_per_s=net,
            source_archive=str(archive.relative_to(ROOT)), checkpoint_sha256=sha(checkpoint)))
        checks.append(dict(length_um=length, checkpoint_hash_verified=True, continuous_history_hash_verified=True, initial_distance_recomputed=True, continuous_zero_clock_recomputed=True, budget_sum_residual_per_s=residual))

    # Reconstruct the staggered native velocity at 60 s in the illustrated cell.
    length = 31
    native_path = register(FLOW / 'L31_native.npz')
    summary = read(FLOW / 'L31_summary.json')
    assert sha(native_path) == summary['native_sha256']
    with np.load(native_path) as z:
        t = z['time'].copy()
        i = int(np.argmin(abs(t-60)))
        assert t[i] == 60 and t[i]-t[i-1] == .5
        x = z['x'].copy()
        f = {k: z['field_'+k][i-1:i+1].copy() for k in ['rac', 'rho', 'mca_bound', 'displacement', 'velocity', 'tension', 'left_boundary_displacement', 'right_boundary_displacement']}
    dx = length/len(x)
    b = f['mca_bound'][1]
    friction = parameters['linker_friction'] * b
    threshold = .899061648345644
    switch = .5*(1+np.tanh((f['rac'][0, [0, -1]]-threshold)/parameters['protrusion_switch_scale']))
    commands = parameters['rac_to_rho_coupling']*parameters['protrusion_speed']*switch/(1+f['tension'][0, [0, -1]]**2/parameters['protrusion_tension_scale'])
    commands[1] *= -1
    membrane_v = (f['displacement'][1]-f['displacement'][0])/.5
    velocities = dict(
        endpoint=solve(friction, parameters['cortical_viscosity'], np.zeros_like(x), dx, *commands),
        rho=solve(friction, parameters['cortical_viscosity'], parameters['active_stress']*parameters['rho_to_rac_coupling']*np.gradient(f['rho'][0], dx), dx),
        membrane=solve(friction, parameters['cortical_viscosity'], friction*membrane_v, dx))
    total = f['velocity'][1]
    mechanical_residual = float(np.max(abs(sum(velocities.values())-total)))
    assert mechanical_residual < 1e-11
    profile_rows = [dict(length_um=31, time_s=60, reference_position_um=float(xx), membrane_displacement_um=float(f['displacement'][1, j]), tension_pN=float(f['tension'][1, j]),
        rho_forcing_concentration=float(f['rho'][0, j]), bound_mca_concentration=float(b[j]), membrane_velocity_um_per_s=float(membrane_v[j]),
        endpoint_cortical_velocity_um_per_s=float(velocities['endpoint'][j]), rho_cortical_velocity_um_per_s=float(velocities['rho'][j]), membrane_cortical_velocity_um_per_s=float(velocities['membrane'][j]), total_cortical_velocity_um_per_s=float(total[j])) for j, xx in enumerate(x)]
    writecsv(STEM + '_length_measurements.csv', measured)
    writecsv(STEM + '_snapshot_profiles.csv', profile_rows)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig = plt.figure(figsize=(13.5, 8.7))
    grid = fig.add_gridspec(2, 3, left=.075, right=.975, bottom=.17, top=.88, wspace=.43, hspace=.62, height_ratios=[1.17, 1])
    a = fig.add_subplot(grid[0, :2])
    panel(a, 'A', 'MCA couples the moving cortex to the membrane')
    a.text(.02, .98, 'Native reversal snapshot: L = 31 µm, 60 s after cue onset', transform=a.transAxes, fontsize=9, va='top', color='#555555')
    nodes = np.linspace(0, length, 17)
    ul, ur = f['left_boundary_displacement'][1, 0], f['right_boundary_displacement'][1, 0]
    u = np.interp(nodes, np.r_[0, x, length], np.r_[ul, f['displacement'][1], ur])
    positions = nodes+u
    tension = np.interp(.5*(nodes[:-1]+nodes[1:]), x, f['tension'][1])
    bound = np.interp(nodes, x, b)
    vv = np.interp(nodes, np.r_[0, x, length], np.r_[commands[0], total, commands[1]])
    a.scatter(nodes, np.full(17, 1.48), s=9, color='#bcbcbc')
    for left, right, value in zip(positions[:-1], positions[1:], tension):
        spring(a, left, right, 1.02, plt.get_cmap('Reds')(Normalize(0, 10)(value)))
        dashpot(a, [left, .05], [right, .05], COLORS['cortex'], size=.055)
    for position, value in zip(positions, bound):
        dashpot(a, [position, .1], [position, .96], COLORS['mca'], width=.8+9*value, size=.09)
    a.scatter(positions, np.full(17, 1.02), s=15, facecolor='white', edgecolor='#666666', lw=.7, zorder=5)
    a.scatter(positions, np.full(17, .05), s=12, color=COLORS['cortex'], marker='s', zorder=5)
    a.quiver(positions, np.full(17, -.24), 30*vv, np.zeros(17), angles='xy', scale_units='xy', scale=1, color=COLORS['cortex'], width=.0035)
    a.quiver(positions[[0, -1]], np.full(2, .05), 30*commands, np.zeros(2), angles='xy', scale_units='xy', scale=1, color=COLORS['endpoint'], width=.006, zorder=6)
    physical_x = x+f['displacement'][1]
    for field, base_y, scale, color, unit in [(f['rho'][0], -.98, .7, COLORS['rho'], '1'), (b, -1.93, 7, COLORS['mca'], '0.1')]:
        a.axhline(base_y, color='#dddddd', lw=.7)
        a.fill_between(physical_x, base_y, base_y+scale*field, color=color, alpha=.13)
        a.plot(physical_x, base_y+scale*field, color=color, lw=1.7)
        a.plot([34.4, 34.4], [base_y, base_y+.7], color=color, lw=1.2)
        a.text(34.8, base_y+.35, unit, color=color, va='center', fontsize=8)
    a.annotate('MCA drag', xy=(positions[8], .60), xytext=(positions[8]+4, .71), color=COLORS['mca'], fontsize=9, arrowprops=dict(arrowstyle='->', color=COLORS['mca'], lw=.9))
    a.text(9, -.42, 'Cortical flow: 3 µm arrow = 0.1 µm/s', fontsize=8, color=COLORS['cortex'])
    a.set(xlim=(-5.2, 36), ylim=(-2.08, 2.65), xlabel='Physical position (µm)', xticks=[0, 10, 20, 30], yticks=[1.48, 1.02, .05, -.69, -1.59], yticklabels=['Reference', 'Membrane', 'Cortex', 'Rho', 'Bound MCA'])
    a.tick_params(axis='y', length=0, labelsize=9)
    a.spines['left'].set_visible(False)
    a.legend(handles=[Line2D([], [], color='#b51b1b', lw=2, label='Membrane springs (tension)'), Line2D([], [], color=COLORS['cortex'], lw=2, label='Cortical dashpots'), Line2D([], [], color=COLORS['endpoint'], marker='>', ls='', label='Rac-controlled endpoint motion')], loc='upper center', bbox_to_anchor=(.52, .91), ncol=3, frameon=False, fontsize=8, columnspacing=1.6)

    ax = fig.add_subplot(grid[0, 2])
    panel(ax, 'B', 'Components of cortical flow')
    for key, label, ls in [('endpoint', 'Endpoint forcing', '-'), ('rho', 'Rho stress gradient', '-'), ('membrane', 'Membrane motion', '-'), ('total', 'Total flow', '--')]:
        y = total if key == 'total' else velocities[key]
        ax.plot(x, y, color=COLORS[key], ls=ls, lw=2 if key == 'total' else 1.6, label=label)
    ax.axhline(0, color='#aaaaaa', lw=.7)
    ax.grid(alpha=.13)
    ax.set(xlabel='Reference position (µm)', ylabel='Cortical velocity (µm/s)', xlim=(0, 31))
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    ax.text(.04, .96, 'Same native state as A\nPositive velocity is rightward', transform=ax.transAxes, va='top', fontsize=8, color='#555555')

    ax = fig.add_subplot(grid[1, 0])
    panel(ax, 'C', 'MCA reorientation at 60 s')
    curves = [('endpoint_advection_rate_per_s', 'endpoint', 'Endpoint advection', '-'), ('rho_advection_rate_per_s', 'rho', 'Rho advection', '-'), ('membrane_advection_rate_per_s', 'membrane', 'Membrane-motion advection', '-'), ('nonadvection_rate_per_s', 'nonadvection', 'Diffusion + exchange + restoration', '-'), ('net_rate_per_s', 'total', 'Net MCA reorientation', '--')]
    for field, color, label, ls in curves:
        ax.plot(lengths, [r[field]*1e4 for r in measured], color=COLORS[color], ls=ls, lw=1.8, marker='o', ms=3, label=label)
    ax.axhline(0, color='#aaaaaa', lw=.7)
    ax.grid(alpha=.13)
    ax.set(xlabel='Cell length (µm)', ylabel='MCA reorientation rate (10⁻⁴ s⁻¹)', xlim=(17, 40), xticks=[18, 23, 31, 39])
    ax.text(.04, .96, 'Positive rates erase the original MCA bias', transform=ax.transAxes, va='top', fontsize=8, color='#555555')
    ax.set_ylim(-2.8, 8)
    budget_handles = ax.get_legend_handles_labels()[0]

    ax = fig.add_subplot(grid[1, 1])
    panel(ax, 'D', 'MCA must redistribute farther')
    ax.plot(lengths, [r['initial_centroid_distance_um'] for r in measured], color=COLORS['mca'], marker='o', ms=4, lw=1.8)
    ax.set(xlabel='Cell length (µm)', ylabel='Initial MCA centroid distance\nfrom the midpoint (µm)', xlim=(17, 40), xticks=[18, 23, 31, 39], ylim=(0, 3))
    ax.grid(alpha=.13)
    ax.text(.05, .96, 'Uninterrupted cue: MCA zero occurs\nwhen the centroid reaches the midpoint', transform=ax.transAxes, fontsize=8, va='top', color='#555555')

    ax = fig.add_subplot(grid[1, 2])
    panel(ax, 'E', 'Speed does not keep pace with length')
    speed = np.array([r['mean_distribution_speed_um_per_s'] for r in measured])
    reference = speed[lengths.index(23)] * np.array(lengths)/23
    ax.plot(lengths, 1000*speed, color=COLORS['mca'], marker='o', ms=4, lw=1.8, label='Measured distribution speed')
    ax.plot(lengths, 1000*reference, color='#777777', ls=':', lw=1.6, label='Speed proportional to length\n(anchored at L = 23 µm)')
    ax.set(xlabel='Cell length (µm)', ylabel='Centroid distance / MCA-zero time\n(10⁻³ µm/s)', xlim=(17, 40), xticks=[18, 23, 31, 39], ylim=(0, 9))
    ax.grid(alpha=.13)
    ax.legend(frameon=False, fontsize=8, loc='upper left')

    fig.suptitle('Figure Y | Mechanics sets cortical flow and MCA redistribution', fontsize=18, x=.075, ha='left', y=.975)
    fig.text(.075, .932, 'Corrected native model · uninterrupted constant-slope cue 0.24 · all 12 eligible lengths', fontsize=11, color='#555555')
    fig.legend(budget_handles, [item[2] for item in curves], loc='lower center', bbox_to_anchor=(.51, .036), ncol=3, frameon=False, fontsize=9, title='C: contributions to the MCA reorientation rate', title_fontsize=9, columnspacing=2)
    fig.text(.075, .018, 'Generated '+generated+' · retained simulations complete · distribution speed is not molecular or cortical velocity', fontsize=8, color='#666666')

    caption = ('*Figure Y* Cortical mechanics drives MCA redistribution, while greater redistribution distance in longer cells is not matched by a proportional increase in distribution speed. '
        'A) One native spring/dashpot illustration at L = 31 µm, 60 s after onset of the uninterrupted constant-slope cue (0.24). Membrane springs are positioned at x + displacement and colored by tension (pale to dark red, 0–10 pN). '
        'Green MCA dashpots couple membrane and cortex; their line width increases with local bound-MCA density. Purple axial dashpots depict cortical viscosity, purple arrows show cortical velocity, and orange arrows show prescribed Rac- and tension-dependent endpoint commands. '
        'Rho and bound-MCA profiles share the physical horizontal coordinate; vertical bars represent 1 and 0.1 native concentration units respectively. The cortex drawing is registered locally below the membrane; connector inclination or length does not represent permanent bond extension. '
        'B) Endpoint-driven, Rho-gradient-driven and membrane-motion-driven components sum to the same native cortical flow shown in A. Positive velocity is rightward. The solver is staggered: native Rho and endpoint forcing at 59.5 s act with MCA friction at 60 s; membrane velocity is the displacement difference over that half-second. '
        'C) Exact native bound-MCA polarity-rate budgets at 60 s for all twelve eligible lengths, using the previously audited half-second source identity. Positive rates erase the original attachment bias. Endpoint, Rho and membrane-motion advection are separate; diffusion, bound/unbound exchange and uniform restoration are grouped. '
        'These are instantaneous source contributions, not independently evolved pathway removals; the chemical and mechanical state remains coupled. '
        'D) Initial bound-MCA centroid distance to the nominal midpoint, recomputed from each exact checkpoint. E) That distance divided by the first uninterrupted-cue MCA-polarity-zero time, independently recomputed with linear interpolation between retained 1-s samples. '
        'The dotted reference is speed proportional to cell length, anchored at L = 23 µm. D–E use nominal reference coordinates and all twelve eligible lengths. This speed measures distribution movement, not tracked molecular speed or cortical velocity. '
        'The continuous-cue crossing time is not a finite-pulse success threshold. Cues remain on for these measurements; no finite withdrawal or recovery is mixed into the clocks. No new long integrations; current report unchanged.')
    output(STEM+'_caption.txt').write_text(caption+'\n')
    artifacts = {}
    for extension in ['png', 'pdf', 'svg']:
        path = output(STEM+'.'+extension)
        fig.savefig(path, dpi=220, facecolor='white')
        artifacts[layout.dev_artifact_key(OUT, path)] = sha(path)
    plt.close(fig)
    for name in [STEM+'_caption.txt', STEM+'_length_measurements.csv', STEM+'_snapshot_profiles.csv']:
        path = output(name)
        artifacts[layout.dev_artifact_key(OUT, path)] = sha(path)
    assert all(sha(REPORT/name) == digest for name, digest in report_hashes.items())
    audit = dict(passed=True, lengths=lengths, native_snapshot_length_um=31, snapshot_time_s=60,
        maximum_velocity_reconstruction_residual_um_per_s=mechanical_residual, checks=checks)
    manifest = dict(artifact_root='../..', generated_at=generated, producer='finished', legend=caption, sha256=artifacts,
        source_sha256=SOURCES, parameters=parameters, audit=audit, new_long_integrations=0, current_report_changed=False,
        reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python -B '+str(Path(__file__))])
    output(STEM+'.manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(OUT/(STEM+'.png'))
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
