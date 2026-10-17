"""Full-native-field convergence diagnostics shared by corrected basal assays."""
import numpy as np

def monitored(tr, key, length):
    out = {k: np.asarray(v) for k, v in tr.fields.items()}
    if key == 'debelly':
        out = {k: out[k] for k in ['rac', 'rho', 'mca_bound', 'mca_unbound',
                                  'displacement', 'velocity', 'tension', 'boundary_separation']}
        out['displacement'] = out['displacement'] - out['displacement'].mean(axis=1)[:, None]
        out['physical_length'] = length + out.pop('boundary_separation')
    return out

def deviations(fields, target, scales):
    count = len(next(iter(fields.values())))
    distance = np.zeros(count)
    for name, y in fields.items():
        error = np.abs(y - target[name]).reshape(count, -1).max(axis=1) / scales[name]
        distance = np.maximum(distance, error)
    return distance

def diagnostics(tr, key, length, scales, rule):
    keep = tr.time >= tr.time[-1] - rule['chunk']
    t = tr.time[keep]
    metrics = {}
    for name, y in monitored(tr, key, length).items():
        y = y[keep]
        variation = float(np.ptp(y, axis=0).max()) / scales[name]
        drift = float(abs(y[-1] - y[-2]).max()) / (t[-1] - t[-2]) * rule['chunk'] / scales[name]
        metrics[name] = dict(relative_range=variation, projected_relative_drift=drift)
    return dict(steady=all(max(v.values()) <= rule['tolerance'] for v in metrics.values()), metrics=metrics)
