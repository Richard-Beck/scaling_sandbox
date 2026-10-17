"""Independently check plotted OLS fits against the authoritative phase/case data."""
from pathlib import Path
import csv, gzip, json
import numpy as np
import runner as r


def main():
    out = r.CODE / 'results'
    records = json.load(gzip.open(out / 'records.json.gz', 'rt'))
    groups = {}
    for prep in records['preparations']:
        for i, phase in enumerate(prep['phases']):
            if not phase['settled']:
                continue
            if i == 0 and prep['assay'] == 'gradient_reversal':
                metric = 'settle_conditioning'
            elif i == 1 and prep['assay'] == 'finite_countercue':
                metric = 'settle_finite_zero'
            else:
                continue
            groups.setdefault((prep['model'], metric), []).append((prep['length'], phase['time_to_settle']))
    for case in records['cases']:
        if case['status'] == 'success':
            groups.setdefault((case['model'], case['assay']), []).append((case['length'], case['endpoint']))
    fits = list(csv.DictReader((out / 'length_regressions.csv').open()))
    eligible = 0
    for fit in fits:
        points = groups.get((fit['model'], fit['metric']), [])
        assert int(fit['n_valid']) == len(points), fit
        assert (fit['fit_eligible'] == 'True') == (len(points) > 3), fit
        if len(points) <= 3:
            assert fit['slope'] == '' and fit['normalized_slope'] == '', fit
            continue
        eligible += 1
        x, y = np.asarray(points).T
        matrix = np.column_stack([x, np.ones(len(x))])
        slope, intercept = np.linalg.lstsq(matrix, y, rcond=None)[0]
        np.testing.assert_allclose([float(fit['slope']), float(fit['intercept'])], [slope, intercept], rtol=1e-10, atol=1e-10)
        length = (x.min() + x.max()) / 2
        response = slope * length + intercept
        np.testing.assert_allclose([float(fit['characteristic_length']), float(fit['characteristic_response'])], [length, response], rtol=1e-10, atol=1e-10)
        if abs(response) > 1e-12:
            np.testing.assert_allclose(float(fit['normalized_slope']), slope * length / response, rtol=1e-10, atol=1e-10)
        else:
            assert fit['normalized_slope'] == '', fit
        spread = np.sum((y - y.mean()) ** 2)
        if spread > 1e-20:
            expected_r2 = 1 - np.sum((matrix @ np.array([slope, intercept]) - y) ** 2) / spread
            np.testing.assert_allclose(float(fit['r_squared']), expected_r2, rtol=1e-10, atol=1e-10)
        else:
            assert fit['r_squared'] == '', fit
    assert len(fits) == 28
    result = dict(generated_at=r.now(), passed=True, series_checked=len(fits), eligible_fits_checked=eligible,
                  estimator='Independent numpy.linalg.lstsq against complete records',
                  records_sha256=r.sha(out / 'records.json.gz'), fits_sha256=r.sha(out / 'length_regressions.csv'))
    r.write(out / 'regression_validation.json', result)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
