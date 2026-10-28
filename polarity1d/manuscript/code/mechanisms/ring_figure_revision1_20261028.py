"""Shared display conventions for the R1/R2 revision; no data transformation."""
import numpy as np

SIZE_COLORS = {20: '#E69F00', 40: '#009E73', 60: '#0072B2',
               1: '#E69F00', 2: '#009E73', 3: '#0072B2'}


def pair_rotation(matrices):
    """One angular shift for a size pair, centring its initial dominant domain.

    Initial domain centres receive equal weight across sizes. The same angular
    shift, rounded to each display grid, is applied at every time and size.
    """
    peaks = []
    for v in matrices:
        activity = v[0] - v[0].min()
        mask = activity > .5 * activity.max()
        # Find the dominant connected half-height region on a periodic grid.
        start = int(np.argmin(mask))
        runs = np.split(np.roll(np.arange(len(mask)), -start),
                        np.flatnonzero(np.diff(np.roll(mask, -start))) + 1)
        regions = [r for r in runs if mask[r[0]]]
        region = max(regions, key=lambda r: activity[r].sum())
        centre = np.sum(activity[region] * np.exp(2j*np.pi*(region+.5)/len(mask)))
        peaks.append(np.angle(centre)/(2*np.pi))
    peaks = np.array(peaks)
    centre = np.mean(np.exp(2j * np.pi * peaks))
    phase = peaks[0] if abs(centre) < 1e-12 else np.angle(centre) / (2*np.pi)
    return float((.5 - phase + .5) % 1 - .5)


def rotate_display(values, turns):
    return np.roll(values, int(np.rint(turns * values.shape[1])), axis=1)
