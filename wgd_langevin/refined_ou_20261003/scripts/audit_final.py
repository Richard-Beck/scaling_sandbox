"""Independent final integrity and objective audit; run after all producers stop."""
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
ROOT = RUN.parents[1]
sys.path.insert(0, str(HERE))
import fit
import model


def main():
    checks = {}
    prov = json.loads((RUN / 'data/source_provenance.json').read_text())
    for name, expected in prov['source_sha256'].items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        assert actual == expected, name
    checks['reviewed_input_source_hashes'] = len(prov['source_sha256'])
    curve = pd.read_csv(RUN / 'data/xuan_msd_curve.csv')
    assert len(curve) == 510
    for group, part in curve.groupby('group'):
        assert set(part.lag_bin) == {0, 1, 2}
        assert np.isclose(part.objective_weight.sum(), 1)
        assert np.allclose(part.groupby('lag_bin').objective_weight.sum(), 1 / 3)
        assert part.lag_h.min() < 1 and part.lag_h.max() > 10
    checks['full_curves_and_balanced_phase_weights'] = True
    expected_prov = model.provenance()
    records = []
    for index in range(40):
        path = RUN / 'results' / f'fit_{index:02}.json'
        obj = json.loads(path.read_text())
        assert obj['provenance'] == expected_prov, path.name
        z = np.asarray(obj['z'])
        assert np.all(np.isfinite(z))
        assert np.all(z >= model.LOW) and np.all(z <= model.HIGH)
        residual = fit.residual(z, obj['scenario'])
        score = float(residual @ residual)
        assert np.isclose(score, obj['total_score'], rtol=1e-10, atol=1e-8), path.name
        assert obj['scenario'] == list(model.SCENARIOS)[index // 8]
        assert obj['start'] == index % 8
        records.append(dict(file=path.name, scenario=obj['scenario'],
                            success=obj['success'], stored_score=obj['total_score'],
                            recomputed_score=score))
    checks['all_40_bounded_start_objectives_recomputed'] = True
    checks['converged_start_count'] = sum(r['success'] for r in records)
    selected = json.loads((RUN / 'results/selected_fits.json').read_text())
    sys.path.insert(0, str(HERE / 'hosny_sensitivity'))
    from fast_residuals import compile_residuals
    compiled = compile_residuals(fit)
    rng = np.random.default_rng(20261003)
    base_z = np.asarray(selected['primary']['z'])
    for index in range(12):
        z = base_z.copy()
        if index:
            z[:20] += rng.normal(0, 0.5, 20)
            z[20:] += rng.normal(0, 2, 2)
        prediction = fit.analytic(z)
        original, rows = fit.data_residuals(prediction)
        fast, identities = compiled(prediction)
        assert np.array_equal(original, fast)
        assert [(r['study'], r['group'], r['metric'], r.get('quantile')) for r in rows] == [
            (r['study'], r['group'], r['metric'], r.get('quantile')) for r in identities]
    checks['compiled_sensitivity_residuals_exactly_equal_at_twelve_parameter_checks'] = True
    for scenario in model.SCENARIOS:
        candidates = [r for r in records if r['scenario'] == scenario]
        assert np.isclose(selected[scenario]['total_score'],
                          min(r['stored_score'] for r in candidates))
        validation = json.loads((RUN / 'results' / f'validation_{scenario}.json').read_text())
        assert validation['n'] == 16384
        assert validation['seeds'] == [74101, 92317]
        assert np.array_equal(validation['z'], selected[scenario]['z'])
        rows = validation['comparison']
        assert sum(r['metric'] == 'msd_curve' for r in rows) == 510
        early = [r for r in rows if r['study'] == 'Xuan2018' and r['metric'] == 'msd']
        assert len(early) == 2
        assert all(r['role'] == 'correlated curve check' for r in early)
    checks['selected_minima_and_two_frozen_ensembles_per_scenario'] = True
    movie = pd.read_csv(RUN / 'diagnostics/movie_track_checks.csv')
    assert len(movie) == 48 and movie.cell_id.nunique() == 8
    assert movie.timing_scenario.nunique() == 3
    assert set(movie.centroid_allowance_px) == {0, 3}
    checks['movie_checks_cover_eight_tracks_three_clocks_two_allowances'] = True
    for name, key, base in [
        ('figures/main_03_migration_OU_provenance.json', 'source_sha256', ROOT),
        ('figures/xuan_population_MSD_diagnostics_provenance.json', 'source_sha256', RUN),
        ('figures/supp_04_Hosny_sensitivity_provenance.json', 'inputs_sha256', ROOT),
        ('diagnostics/movie_provenance.json', 'provenance', ROOT),
    ]:
        provenance = json.loads((RUN / name).read_text())
        for path, expected in provenance[key].items():
            assert hashlib.sha256((base / path).read_bytes()).hexdigest() == expected, path
    checks['all_figure_and_movie_producer_source_hashes'] = True
    sensitivity = RUN / 'results/hosny_sensitivity'
    hosny = json.loads((sensitivity / 'hosny_analysis.json').read_text())
    assert not hosny['smoke_run']
    assert hosny['surface']['all_drift_optimizers_successful']
    assert hosny['surface']['all_fitted_parameters_within_model_bounds']
    assert hosny['surface']['max_optimized_minus_fixed_total_delta'] <= 1e-7
    assert hosny['profile_optimization_audit']['converged_profile_start_count'] == 288
    assert hosny['leaveout_optimization_audit']['selected_fit_successful']
    assert hosny['leaveout_optimization_audit']['all_parameters_within_model_bounds']
    assert hashlib.sha256((RUN / 'results' / hosny['fit']).read_bytes()).hexdigest() == hosny['fit_sha256']
    profiles = pd.read_csv(sensitivity / 'hosny_drift_profiles.csv')
    for row in profiles.itertuples():
        z = np.asarray(json.loads(row.z))
        residual = fit.residual(z, 'primary')
        assert np.isclose(residual @ residual, row.total_score, rtol=1e-10, atol=1e-7)
    checks['full_hosny_surface_profiles_leaveout_and_288_recomputed_objectives'] = True
    focused = json.loads((sensitivity / 'hosny_surface_zoom.json').read_text())
    assert focused['validation']['all_points_converged']
    assert focused['validation']['all_values_finite']
    assert focused['validation']['all_parameters_within_model_bounds']
    assert len(pd.read_csv(sensitivity / 'hosny_surface_zoom.csv')) == 3721
    assert json.loads((sensitivity / 'status.json').read_text())['state'] == 'complete'
    checks['focused_3721_point_surface_and_stopped_producer'] = True
    for stem in ['main_03_migration_OU', 'supp_04_Hosny_sensitivity']:
        for suffix in ['png', 'pdf', 'svg']:
            path = RUN / 'figures' / f'{stem}.{suffix}'
            assert path.stat().st_size > 1000, str(path)
    checks['six_manuscript_exports_exist'] = True
    transient = [str(p.relative_to(RUN)) for p in RUN.rglob('*')
                 if p.is_file() and (p.suffix in {'.tmp', '.pyc'} or p.name.endswith('~'))]
    assert not transient, transient
    checks['no_temporary_or_bytecode_files'] = True
    result = dict(generated_utc=datetime.now(timezone.utc).isoformat(),
                  passed=True, checks=checks, objective_records=records)
    (RUN / 'results/final_audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
