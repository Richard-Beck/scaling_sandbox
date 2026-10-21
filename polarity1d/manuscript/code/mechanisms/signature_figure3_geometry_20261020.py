"""Prospective Figure 3 geometry; continuous centres, unmodified native fields.

The 10% edge is an interpolated contrast contour of the exact initial basal.
Nothing here reclassifies overlapping placements as independent challengers.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import find_peaks

GAPS = {
    "wave_pinning": [0, .5, 1, 2, 3, 4, 5, 6, 8, 10, 12],
    "goryachev": [0, 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20],
    "holmes_model3": [0, .5, 1, 2, 3, 4, 5, 6, 8, 10, 12],
    "spring": [0, .05, .1, .2, .3, .4, .5, .6, .7],
}


def arc(x, center, circumference):
    return (np.asarray(x) - center + circumference / 2) % circumference - circumference / 2


def interval_weights(x, circumference, center, left, right):
    """Finite-volume intersection with one unwrapped physical interval."""
    a = arc(x, center, circumference)
    dx = circumference / len(x)
    weights = np.zeros(len(x))
    for shift in (-circumference, 0, circumference):
        weights += np.maximum(0, np.minimum(a + dx/2, right + shift)
                              - np.maximum(a - dx/2, left + shift))
    return weights


def components(mask):
    mask = np.asarray(mask, bool)
    if not mask.any():
        return []
    if mask.all():
        return [np.arange(len(mask))]
    result = []
    for start in np.flatnonzero(mask & ~np.roll(mask, 1)):
        group, j = [], int(start)
        while mask[j]:
            group.append(j)
            j = (j + 1) % len(mask)
        result.append(np.asarray(group))
    return result


def contour(x, circumference, primary, fraction=.1):
    """Linear contour crossings around the single initial active region.

    Coordinates are cell centres; finite-volume weights are derived only after
    the physical crossing has been interpolated. A seam-straddling site is
    unwrapped around its actual maximum, never shifted in the saved state.
    """
    p = np.asarray(primary)
    contrast = float(np.ptp(p))
    if contrast <= 0:
        raise ValueError("No active-region contrast")
    peak = int(np.argmax(p))
    origin = float(x[peak])
    s = arc(x, origin, circumference)
    order = np.argsort(s)
    s, u = s[order], p[order]
    level = float(p.min() + fraction * contrast)
    inside = np.flatnonzero(u >= level)
    if len(inside) == 0 or np.any(np.diff(inside) != 1):
        raise ValueError("Initial contrast contour is not one connected region")
    a, b = int(inside[0]), int(inside[-1])
    if a == 0 or b == len(s) - 1:
        raise ValueError("Initial contour reaches the opposite seam")
    left = float(np.interp(level, u[a-1:a+1], s[a-1:a+1]))
    right = float(np.interp(level, u[b:b+2][::-1], s[b:b+2][::-1]))
    weights = interval_weights(x, circumference, origin, left, right)
    np.testing.assert_allclose(weights.sum(), right-left, atol=1e-10)
    return dict(origin=origin, peak_cell=peak, left=left, right=right,
                width=right-left, background=float(p.min()), contrast=contrast,
                threshold_fraction=fraction, threshold=level, weights=weights)


def serial_contour(value):
    return {k: v for k, v in value.items() if k != "weights"}


def architecture(primary):
    p = np.asarray(primary)
    contrast = float(np.ptp(p))
    peaks, properties = find_peaks(np.tile(p, 3), prominence=.05)
    selected = (peaks >= len(p)) & (peaks < 2*len(p))
    return dict(contrast=contrast, half_contrast_regions=len(components(p >= p.min()+.5*contrast))
                if contrast > 1e-12 else 0,
                resolved_peaks=int(selected.sum()), peak_cells=(peaks[selected]-len(p)).tolist(),
                peak_prominences=properties["prominences"][selected].tolist())


def classification(left, right, center_distance, support):
    near, far = center_distance-support/2, center_distance+support/2
    overlap = max(0., min(right, far)-max(left, near))
    if overlap <= 1e-10:
        return "external"
    if near >= left-1e-10 and far <= right+1e-10:
        return "internal"
    return "borderline"


def compile_placements(model, x, circumference, primary, support):
    edges = {f: contour(x, circumference, primary, f) for f in (.1, .5)}
    edge = edges[.1]
    maximum = circumference/2-edge["right"]-support/2
    if maximum < -1e-10:
        raise ValueError("No external cue placement at the antipode")
    prescribed = [(float(g), "absolute_gap_grid") for g in GAPS[model] if g < maximum-1e-10]
    prescribed.append((float(maximum), "antipode"))
    prescribed += [(g, "requested_borderline") for g in ((-.07, -.13) if model == "spring" else (-1., -2.))]
    if model != "goryachev":
        # Exactly centred internal reference, if the support actually fits.
        prescribed.append((float(-support/2-edge["right"]), "internal_reference"))
    rows = []
    for gap, role in prescribed:
        distance = gap+edge["right"]+support/2
        if distance < -1e-10 or distance > circumference/2+1e-10:
            rows.append(dict(placement_index=len(rows), requested_role=role,
                placement_status="geometrically_infeasible_requested_gap",
                requested_gap=float(gap), edge_gap_10=None, edge_gap_50=None,
                classification_10="unavailable", classification_50="unavailable",
                claim_eligible=False, basal_circumference=circumference, support=support,
                requested_center_distance=distance,
                attainable_gap_range=[float(-edge["right"]-support/2), float(maximum)],
                closest_attainable_gap=float(np.clip(gap, -edge["right"]-support/2, maximum)),
                reason="Requested gap requires negative centre distance or passes the antipode; no surrogate placement is run",
                initial_10_percent_contour=serial_contour(edge),
                initial_50_percent_contour=serial_contour(edges[.5])))
            continue
        if abs(distance) < 1e-10:
            distance = 0.
        near = distance-support/2
        far = distance+support/2
        center = float((edge["origin"]+distance) % circumference)
        cls = classification(edge["left"], edge["right"], distance, support)
        cls50 = classification(edges[.5]["left"], edges[.5]["right"], distance, support)
        rows.append(dict(placement_index=len(rows), requested_role=role, placement_status="feasible",
            edge_gap_10=float(near-edge["right"]), edge_gap_50=float(near-edges[.5]["right"]),
            classification_10=cls, classification_50=cls50, claim_eligible=cls == "external",
            cue_center_physical=center, cue_center_fraction=center/circumference,
            center_distance=distance, target_angle_radians=2*np.pi*distance/circumference,
            short_route=distance, long_route=circumference-distance,
            short_edge_gap=float(near-edge["right"]),
            long_edge_gap=float(circumference-distance+edge["left"]-support/2),
            cue_near_edge_relative=near, cue_far_edge_relative=far,
            support=support, basal_circumference=circumference,
            initial_10_percent_contour=serial_contour(edge),
            initial_50_percent_contour=serial_contour(edges[.5]),
            center_evaluation="continuous physical/material centre; no mesh snapping",
            requested_class_matches_actual=(role != "requested_borderline" or cls == "borderline")
                and (role != "internal_reference" or cls == "internal"),
            internal_reference_support_fits=role != "internal_reference" or cls == "internal"))
    return rows


def region_observables(x, circumference, primary, basal_primary, placement, fraction=.5):
    """Frozen incumbent/cue-support excess amounts, no lineage interpretation.

    The baseline is each initial basal's minimum primary concentration. The
    fixed incumbent region uses the specified initial contrast contour. Spring
    physical weights scale with instantaneous circumference through x and C.
    """
    initial_C = placement["basal_circumference"]
    initial_x = np.asarray(x) * initial_C/circumference
    core = contour(initial_x, initial_C, basal_primary, fraction)
    weights = core["weights"] * circumference/initial_C
    center = placement["cue_center_fraction"] * circumference
    cue_weights = interval_weights(x, circumference, center, -placement["support"]/2, placement["support"]/2)
    background = core["background"]
    old = float((primary-background) @ weights)
    challenger = float((primary-background) @ cue_weights)
    initial_old = float((basal_primary-background) @ core["weights"])
    return dict(incumbent_excess=old, challenger_excess=challenger,
                basal_incumbent_excess=initial_old,
                incumbent_fraction=old/initial_old, challenger_fraction=challenger/initial_old,
                excess_ratio=challenger/old if old > 0 else None,
                dominance=(challenger-old)/(challenger+old) if challenger+old > 0 else None)


def measured_regions(x, circumference, primary, incumbent_origin, cue_origin, fraction=.5):
    """Actual two-region gap at withdrawal, with ambiguous geometries explicit."""
    p = np.asarray(primary)
    level = float(p.min()+fraction*np.ptp(p))
    groups = components(p >= level)
    if len(groups) != 2:
        return dict(available=False, reason="not_exactly_two_regions", regions=len(groups), fraction=fraction)
    descriptors = []
    dx = circumference/len(x)
    for group in groups:
        # Crossing positions follow the same linear convention as initial edges.
        start, stop = int(group[0]), int(group[-1])
        prev, nxt = (start-1) % len(x), (stop+1) % len(x)
        left = float(x[prev]+(level-p[prev])/(p[start]-p[prev])*dx)
        length_to_stop = ((stop-start) % len(x))*dx
        right = float(x[start]+length_to_stop+(level-p[stop])/(p[nxt]-p[stop])*dx)
        left_unwrapped = float(x[start]-(1-(level-p[prev])/(p[start]-p[prev]))*dx)
        width = right-left_unwrapped
        midpoint = (left_unwrapped+width/2) % circumference
        descriptors.append(dict(left=left % circumference, width=width, midpoint=midpoint,
                                peak_cell=int(group[np.argmax(p[group])]), cells=group.tolist()))
    d = abs(float(arc(descriptors[1]["midpoint"], descriptors[0]["midpoint"], circumference)))
    gap = d-(descriptors[0]["width"]+descriptors[1]["width"])/2
    old = int(np.argmin([abs(float(arc(r["midpoint"], incumbent_origin, circumference))) for r in descriptors]))
    challenger = int(np.argmin([abs(float(arc(r["midpoint"], cue_origin, circumference))) for r in descriptors]))
    available = old != challenger
    return dict(available=available, reason="two_spatially_associated_regions" if available else "association_ambiguous",
                regions=2, fraction=fraction, measured_gap=gap, center_distance=d,
                long_center_distance=circumference-d, incumbent_index=old, challenger_index=challenger,
                region_descriptors=descriptors)


def activity_vector(primary, fraction_positions, direction_floor=1e-10):
    p = np.asarray(primary)
    denominator = float(p.sum())
    z = complex(p @ np.exp(2j*np.pi*np.asarray(fraction_positions))/denominator) if denominator else 0j
    return dict(real=z.real, imag=z.imag, magnitude=abs(z),
                heading_radians=float(np.angle(z)) if abs(z) > direction_floor else None,
                direction_available=abs(z) > direction_floor,
                magnitude_eligibility_005=abs(z) >= .05,
                magnitude_eligibility_010=abs(z) >= .1,
                magnitude_eligibility_020=abs(z) >= .2)
