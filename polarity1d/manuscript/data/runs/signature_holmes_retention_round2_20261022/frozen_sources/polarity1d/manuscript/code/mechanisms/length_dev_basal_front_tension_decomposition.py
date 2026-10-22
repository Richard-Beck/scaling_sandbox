"""Plot the frozen-state stationary front-tension balance; no integrations."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json

import numpy as np
from scipy.linalg import solve_banded
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import length_dev_output_layout as layout

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'figures/length_mechanisms_20261008'
OUT = REPORT / 'dev/resting_state'
BASE = ROOT / 'data/runs/debelly_corrected_three_test_20261008'
STEM = 'basal_endpoint_front_tension_decomposition'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def solve(diagonal, coefficient, rhs, dx, left=0., right=0.):
    """Independently assemble the native cell-centred Dirichlet operator."""
    n = len(diagonal)
    ratio = coefficient / dx**2
    matrix = np.zeros((3, n))
    matrix[1] = diagonal + 2 * ratio
    matrix[1, [0, -1]] += ratio
    matrix[0, 1:] = -ratio
    matrix[2, :-1] = -ratio
    adjusted = rhs.copy()
    adjusted[0] += 2 * ratio * left
    adjusted[-1] += 2 * ratio * right
    return solve_banded((1, 1), matrix, adjusted)


def gradient(values, dx, left, right):
    result = np.empty_like(values)
    result[1:-1] = (values[2:] - values[:-2]) / (2 * dx)
    result[0] = (values[1] + values[0] - 2 * left) / (2 * dx)
    result[-1] = (2 * right - values[-1] - values[-2]) / (2 * dx)
    return result


def main():
    sources = {}

    def register(path):
        sources[str(path.relative_to(ROOT))] = sha(path)
        return path

    def read(path):
        return json.loads(register(path).read_text())

    report_hashes = {name: sha(REPORT / name) for name in ['index.html', 'length_mechanisms.pdf']}
    meta = read(OUT / 'supporting/metadata/basal_endpoint_metrics_vs_length.manifest.json')
    metrics_path = register(OUT / 'supporting/data/basal_endpoint_metrics_vs_length.csv')
    assert sha(metrics_path) == meta['sha256']['supporting/data/basal_endpoint_metrics_vs_length.csv']
    parameters = meta['parameters']
    assert not parameters['center_boundary_motion']
    campaign = read(BASE / 'manifest.json')
    assert not campaign['center_boundary_motion']
    assert read(BASE / 'results/independent_audit.json')['passed']
    with metrics_path.open(newline='') as handle:
        metrics = sorted(csv.DictReader(handle), key=lambda row: float(row['length_um']))
    eligible = sorted({int(c['length']) for c in campaign['cases'] if c['model'] == 'debelly' and c['eligible']})
    assert eligible == [int(row['length_um']) for row in metrics]
    k = parameters['membrane_stiffness']
    eta = parameters['cortical_viscosity']
    zeta = parameters['active_stress'] * parameters['rho_to_rac_coupling']
    rows, audits = [], []
    for row in metrics:
        length = int(row['length_um'])
        checkpoint = register(ROOT / row['checkpoint'])
        assert sha(checkpoint) == row['checkpoint_sha256']
        siblings = [c for c in campaign['cases'] if c['model'] == 'debelly' and c['length'] == length]
        assert len(siblings) == 3 and all(c['eligible'] and c['checkpoint_sha256'] == sha(checkpoint) for c in siblings)
        with np.load(checkpoint) as archive:
            fields = {key[6:]: archive[key].copy() for key in archive.files if key.startswith('state_')}
        dx = length / len(fields['rac'])
        zero = np.zeros_like(fields['rac'])
        friction = parameters['linker_friction'] * np.maximum(fields['mca_bound'], 1e-10)
        left_u = float(fields['left_boundary_displacement'][0])
        right_u = float(fields['right_boundary_displacement'][0])
        switch = .5 * (1 + np.tanh((fields['rac'][[0, -1]] - float(row['rac_switch_threshold'])) / parameters['protrusion_switch_scale']))
        commands = parameters['rac_to_rho_coupling'] * parameters['protrusion_speed'] * switch / (1 + fields['tension'][[0, -1]]**2 / parameters['protrusion_tension_scale'])
        commands[1] *= -1
        np.testing.assert_allclose(commands, [float(row['derived_left_command_signed_um_per_s']), float(row['derived_right_command_signed_um_per_s'])], atol=1e-15)
        endpoint_v = solve(friction, eta, zero, dx, *commands)
        rho_v = solve(friction, eta, zeta * np.gradient(fields['rho'], dx), dx)
        extension_u = solve(zero, k, zero, dx, left_u, right_u)
        endpoint_drag_u = solve(zero, k, friction * endpoint_v, dx)
        rho_drag_u = solve(zero, k, friction * rho_v, dx)
        extension_t = k * gradient(extension_u, dx, left_u, right_u)
        endpoint_drag_t = k * gradient(endpoint_drag_u, dx, 0., 0.)
        rho_drag_t = k * gradient(rho_drag_u, dx, 0., 0.)
        total_t = extension_t + endpoint_drag_t + rho_drag_t
        np.testing.assert_allclose(extension_t, k * (right_u - left_u) / length, atol=1e-10, rtol=1e-10)
        errors = {
            'displacement_max_abs_error_um': float(np.max(abs(extension_u + endpoint_drag_u + rho_drag_u - fields['displacement']))),
            'velocity_max_abs_error_um_per_s': float(np.max(abs(endpoint_v + rho_v - fields['velocity']))),
            'tension_max_abs_error_pN': float(np.max(abs(total_t - fields['tension']))),
        }
        assert max(errors.values()) < 1e-5, (length, errors)
        rows.append(dict(length_um=length, endpoint_stretching_pN=float(extension_t[0]), endpoint_flow_drag_pN=float(endpoint_drag_t[0]), rho_driven_flow_pN=float(rho_drag_t[0]), total_front_tension_pN=float(fields['tension'][0]), reconstructed_front_tension_pN=float(total_t[0]), front_tension_residual_pN=float(total_t[0] - fields['tension'][0]), checkpoint=row['checkpoint'], checkpoint_sha256=sha(checkpoint)))
        audits.append(dict(length_um=length, checkpoint_hash_verified=True, shared_start_verified=True, **errors))

    csvpath = layout.dev_output_path(OUT, STEM + '.csv')
    with csvpath.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig, ax = plt.subplots(figsize=(9, 5.6), layout='constrained')
    for field, label, color, width in [
        ('endpoint_stretching_pN', 'Endpoint stretching', '#0072B2', 2),
        ('endpoint_flow_drag_pN', 'Drag from endpoint-driven flow', '#D55E00', 2),
        ('rho_driven_flow_pN', 'Rho-driven flow', '#009E73', 2),
        ('total_front_tension_pN', 'Total front tension', '#252525', 2.8),
    ]:
        ax.plot(eligible, [row[field] for row in rows], color=color, lw=width, marker='o', ms=4.5, label=label)
    ax.axhline(0, color='#888888', lw=.85, zorder=0)
    ax.grid(alpha=.14, lw=.6)
    ax.set(xlabel='Cell length (µm)', ylabel='Old-front tension contribution (pN)', title='Resting old-front tension: mechanical contributions', xlim=(17.5, 39.5), xticks=eligible)
    ax.tick_params(axis='x', labelsize=9)
    ax.legend(frameon=False, loc='lower left', fontsize=10)
    artifacts = {layout.dev_artifact_key(OUT, csvpath): sha(csvpath)}
    for extension in ['png', 'pdf', 'svg']:
        path = layout.dev_output_path(OUT, STEM + '.' + extension)
        fig.savefig(path, dpi=240, facecolor='white')
        artifacts[layout.dev_artifact_key(OUT, path)] = sha(path)
    plt.close(fig)
    caption = ('Figure D1 Increasing endpoint-driven cortical drag produces a downturn in resting old-front tension despite increasing endpoint stretching. '
        'One facet shows all 12 eligible corrected basal checkpoints (L = 18–39 µm). Blue: endpoint stretching, the tension of the affine displacement profile connecting the stored endpoints without distributed drag. '
        'Orange: membrane tension generated by MCA drag from cortical flow driven solely by the native endpoint commands, with zero Rho stress gradient. '
        'Green: membrane tension generated by MCA drag from Rho-gradient-driven cortical flow with zero endpoint commands. '
        'The orange and green membrane responses each have zero endpoint displacement. All responses retain each length’s native bound-MCA friction and Rho profile. '
        'Black: stored native tension in the first finite-volume cell at the old left front. The three contributions sum to its stationary reconstruction. '
        'This is a linear mechanical decomposition at frozen native chemistry, not independently evolved removal of coupled pathways. Only eligible starts are included; no new simulation or change to the report.')
    caption_path = layout.dev_output_path(OUT, STEM + '_caption.txt')
    caption_path.write_text(caption + '\n')
    artifacts[layout.dev_artifact_key(OUT, caption_path)] = sha(caption_path)
    register(Path(__file__))
    register(Path(layout.__file__))
    assert all(sha(REPORT / name) == digest for name, digest in report_hashes.items())
    record = dict(artifact_root='../..', generated_at=datetime.now(timezone.utc).isoformat(), producer='finished', legend=caption, sha256=artifacts, source_sha256=sources, parameters=parameters, new_integrations=0, current_report_changed=False, audit=dict(passed=True, lengths=eligible, checks=audits, maximum_front_tension_residual_pN=max(abs(row['front_tension_residual_pN']) for row in rows)), reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python -B ' + str(Path(__file__))])
    layout.dev_output_path(OUT, STEM + '.manifest.json').write_text(json.dumps(record, indent=2) + '\n')
    print(OUT / (STEM + '.png'))
    print(json.dumps(record['audit'], indent=2))


if __name__ == '__main__':
    main()
