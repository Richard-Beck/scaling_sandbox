"""Descriptive parameter comparison with the selected historical manuscript fit."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

RUN = Path(__file__).resolve().parent.parent
ROOT = RUN.parents[1]
OLD = ROOT / 'wgd_langevin/gaussian_hierarchy/results/fit_05.json'


def main():
    old = json.loads(OLD.read_text())
    selected_path = RUN / 'results/selected_fits.json'
    new = json.loads(selected_path.read_text())['primary']
    parameters = [np.exp(np.array(o['z'][:20]).reshape(5, 2, 2)) for o in (old, new)]
    rows = []
    studies = ['Xuan2018', 'Mallin2023', 'Jemaa2017', 'Jemaa2023', 'Hosny_baseline']
    for i, study in enumerate(studies):
        for j, name in enumerate(['tau_h', 'sigma_um_h_minus_3_over_2']):
            a, b = parameters
            rows.append(dict(study=study, parameter=name,
                previous_reference=a[i, 0, j], refined_reference=b[i, 0, j],
                previous_wgd=a[i, 1, j], refined_wgd=b[i, 1, j],
                previous_fold=a[i, 1, j] / a[i, 0, j],
                refined_fold=b[i, 1, j] / b[i, 0, j]))
    pd.DataFrame(rows).to_csv(RUN / 'results/comparison_to_previous.csv', index=False)
    metadata = dict(generated_utc=datetime.now(timezone.utc).isoformat(),
        source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in [OLD, selected_path, Path(__file__).resolve()]},
        previous_pooled_folds=np.exp(np.log(parameters[0][:, 1] / parameters[0][:, 0]).mean(axis=0)).tolist(),
        refined_pooled_folds=np.exp(np.log(parameters[1][:, 1] / parameters[1][:, 0]).mean(axis=0)).tolist(),
        note='Parameter comparison only. Objective scores are not compared because the data targets and weights changed. No historical fit supplies a start or constraint.')
    (RUN / 'results/comparison_to_previous.json').write_text(json.dumps(metadata, indent=2) + '\n')


if __name__ == '__main__':
    main()
