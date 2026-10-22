"""Proposed complete-active-region readouts for review before dependent runs.

Two resolved native peaks define two periodic watershed basins bounded by the
minimum on each connecting route. Integrate the complete piecewise-linear
excess profile in each basin, above the common instantaneous field minimum.
Cue support is never a challenger-amount ROI. One moving or merged region is
reported separately and is never invented as a two-region competition event.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import csv
import json

import numpy as np
from scipy.signal import find_peaks
import signature_figure3_geometry_20261020 as geometry
import signature_figure3_native_20261020 as native


def interval_nodes(primary, circumference, left, width):
    """Linear periodic reconstruction at exact fractional interval boundaries."""
    p = np.asarray(primary, float)
    n, C = len(p), float(circumference)
    dx = C/n
    if width < -1e-12 or width > C+1e-10:
        raise ValueError("Interval width must lie within one circumference")
    left = float(left % C)
    right = left+float(width)
    x = (np.arange(-1, 2*n+2)+.5)*dx
    y = p[np.arange(-1, 2*n+2) % n]
    inside = (x > left) & (x < right)
    xx = np.r_[left, x[inside], right]
    yy = np.r_[np.interp(left, x, y), y[inside], np.interp(right, x, y)]
    return xx, yy


def excess_integral(primary, circumference, left, width, floor, threshold=None):
    """Exact integral of the linear reconstructed profile in a fractional ROI.

    For a thresholded sensitivity, retain the original excess p-floor only
    where p exceeds the declared contour. Exact crossing points are inserted
    in each linear segment; there is no cell-centre inclusion approximation.
    """
    x, y = interval_nodes(primary, circumference, left, width)
    if threshold is None:
        return float(np.sum(np.diff(x)*(y[:-1]+y[1:]-2*floor)/2))
    answer = 0.
    for a, b, u, v in zip(x[:-1], x[1:], y[:-1], y[1:]):
        if u >= threshold and v >= threshold:
            answer += (b-a)*(u+v-2*floor)/2
        elif u < threshold and v < threshold:
            continue
        else:
            crossing = a+(threshold-u)/(v-u)*(b-a)
            answer += ((crossing-a)*(u+threshold-2*floor)/2 if u >= threshold else
                       (b-crossing)*(threshold+v-2*floor)/2)
    return float(answer)


def linear_interval_weights(circumference, cells, left, width):
    """Weights for the same periodic linear reconstruction as excess_integral.

    Fractional basin boundaries stay fixed in material coordinates when these
    physical weights are scaled by a later Spring circumference ratio.
    """
    C, n = float(circumference), int(cells)
    x, _ = interval_nodes(np.zeros(n), C, left, width)
    trap = np.zeros(len(x))
    trap[:-1] += np.diff(x)/2
    trap[1:] += np.diff(x)/2
    coordinate = x/(C/n)-.5
    lower = np.floor(coordinate).astype(int)
    fraction = coordinate-lower
    weights = np.zeros(n)
    np.add.at(weights, lower % n, trap*(1-fraction))
    np.add.at(weights, (lower+1) % n, trap*fraction)
    np.testing.assert_allclose(weights.sum(), width, rtol=1e-13, atol=1e-12)
    assert weights.min() >= -1e-13
    return weights


def frozen_basin_exchange_budget(adapter, state, spatial_readout):
    """Exact cue-free RHS in complete withdrawal basins, with floor frozen.

    This instantaneous exchange observable freezes the measured watershed
    boundaries and common concentration floor at withdrawal. It does not
    differentiate a moving minimum, peak finder or time-dependent watershed.
    Spring material boundaries are fixed while physical arc weights expand.
    """
    if not spatial_readout["two_region_available"]:
        return dict(available=False, reason=spatial_readout["state_class"],
                    exact_cue_free_rhs_required=True)
    C = adapter.circumference(state)
    p = adapter.primary(state)
    rhs = adapter.rhs(state)
    floor = spatial_readout["common_current_floor"]
    if adapter.spring:
        u, v, _ = adapter.model.fields(state)
        from spring_ring_closure_20261017 import native as spring_native
        reaction = spring_native.auth["f"](u, v, 4., 5, 6, 2, 3.+80.*(C-adapter.size))
        diffusion = .01*(np.roll(u, 1)-2*u+np.roll(u, -1))*(adapter.cells/C)**2
        dilution = -u*rhs[-1]/C
    else:
        fields = np.asarray(state).reshape(-1, adapter.cells)
        reaction = adapter.model.reaction(0., adapter.grid.x, fields, np.zeros(adapter.cells))[0]
        diffusion = adapter.model.diffusivities[0]*native.chemical.periodic_laplacian(fields[0], adapter.grid.dx)
        dilution = np.zeros(adapter.cells)
    records = []
    for role in ("incumbent", "challenger"):
        basin = spatial_readout[role+"_basin"]
        weights = linear_interval_weights(C, adapter.cells, basin["basin_left_physical"], basin["basin_width_physical"])
        excess = float((p-floor) @ weights)
        np.testing.assert_allclose(excess, spatial_readout[role+"_excess"], rtol=1e-12, atol=1e-12)
        field_rate = float(rhs[:adapter.cells] @ weights)
        dilation_rate = float(rhs[-1]/C*excess) if adapter.spring else 0.
        exact = field_rate+dilation_rate
        epsilon = 1e-5/max(1., float(np.max(abs(rhs)/adapter.scales)))
        def frozen_amount(y):
            scale = adapter.circumference(y)/C if adapter.spring else 1.
            return float((adapter.primary(y)-floor) @ (weights*scale))
        difference = (frozen_amount(state+epsilon*rhs)-frozen_amount(state-epsilon*rhs))/(2*epsilon)
        error = abs(difference-exact)
        assert error <= 1e-8+1e-5*abs(exact)
        chemistry_rate, diffusion_rate, dilution_rate = [float(v @ weights) for v in (reaction, diffusion, dilution)]
        assert abs(chemistry_rate+diffusion_rate+dilution_rate+dilation_rate-exact) <= 1e-9*(1+abs(exact))
        records.append(dict(role=role, frozen_withdrawal_floor=floor,
            frozen_basin_left_fraction=basin["basin_left_physical"]/C,
            frozen_basin_width_fraction=basin["basin_width_physical"]/C,
            frozen_material_roi_weights=(weights/C).tolist(), withdrawal_excess=excess,
            exact_rate=exact, chemical_rate=chemistry_rate, diffusion_rate=diffusion_rate,
            concentration_dilution_rate=dilution_rate, physical_weight_dilation_rate=dilation_rate,
            directional_finite_difference=difference, finite_difference_error=error,
            finite_difference_epsilon=epsilon))
    return dict(available=True, exact_cue_free_rhs=True, regions=records,
        convention="Complete watershed basins and common floor frozen at withdrawal; exact piecewise-linear weights, Spring physical weights follow circumference; no moving-boundary or floor derivative")


def resolved_peaks(primary, prominence=.05):
    p = np.asarray(primary)
    peaks, properties = find_peaks(np.tile(p, 3), prominence=prominence)
    use = (peaks >= len(p)) & (peaks < 2*len(p))
    return peaks[use]-len(p), properties["prominences"][use]


def valley(primary, first, second):
    """Continuous subcell saddle position, with nonunique active flats flagged.

    A boundary at the discrete minimum cell jumps when neighboring values
    exchange order. The three-point parabolic vertex moves continuously
    through their midpoint at equality. Integrals still use the exact linear
    native-profile reconstruction, not the fitted parabola.
    """
    p = np.asarray(primary)
    n = len(p)
    count = (int(second)-int(first)) % n
    route = int(first)+np.arange(count+1)
    values = p[route % n]
    minimum = float(values.min())
    tied = route[np.isclose(values, minimum, rtol=0., atol=1e-14)]
    # A midpoint among tied minimum cells is a fractional boundary. If distant
    # disjoint ties occur, choose the first tied minimum and retain ambiguity.
    contiguous = len(tied) <= 1 or np.all(np.diff(tied) == 1)
    position = float((tied[0]+tied[-1])/2 if contiguous else tied[0])
    curvature, offset = None, 0.
    if len(tied) == 1:
        k = int(tied[0])
        left, center, right = float(p[(k-1) % n]), float(p[k % n]), float(p[(k+1) % n])
        curvature = left-2*center+right
        if curvature > 0:
            offset = float(np.clip(.5*(left-right)/curvature, -.5, .5))
            position += offset
    active = minimum-float(p.min()) > 1e-8*max(1., float(np.ptp(p)))
    ambiguous = active and (not contiguous or len(tied) > 2 or
        (curvature is not None and curvature <= 64*np.finfo(float).eps*max(1., abs(minimum))))
    return dict(cell_coordinate=position, floor=minimum, tied_minimum_cells=len(tied),
                tied_minima_contiguous=bool(contiguous), subcell_vertex_offset=offset,
                local_curvature=curvature, boundary_convention="Three-point parabolic subcell minimum; midpoint of contiguous flat ties",
                ambiguous_nonunique_active_saddle=bool(ambiguous))


def whole_field_check(primary, circumference, incumbent_fraction, cue_fraction):
    p = np.asarray(primary)
    floor = float(p.min())
    C = float(circumference)
    signed = float(geometry.arc(cue_fraction, incumbent_fraction, 1.))
    if abs(signed) < 1e-12:
        return dict(available=False, reason="coincident_incumbent_and_cue_axes")
    # The cue Voronoi half-ring begins at the short-route angular bisector.
    # Reverse the chosen boundary for a negative directed placement.
    midpoint = (incumbent_fraction+signed/2) % 1
    left = midpoint*C if signed > 0 else (midpoint*C-C/2) % C
    challenger = excess_integral(p, C, left, C/2, floor)
    total = C*float(np.mean(p)-floor)
    incumbent = total-challenger
    return dict(available=total > 1e-12, common_current_floor=floor,
        incumbent_excess=incumbent, challenger_excess=challenger, total_excess=total,
        dominance=(challenger-incumbent)/total if total > 1e-12 else None,
        convention="Complete-field circular Voronoi halves around initial incumbent and cue axes, no beam-footprint ROI; diagnostic only")


def readout(primary, circumference, incumbent_fraction, cue_fraction, prominence=.05, detailed=True):
    p = np.asarray(primary, float)
    C, n = float(circumference), len(p)
    q = (np.arange(n)+.5)/n
    floor, contrast = float(p.min()), float(np.ptp(p))
    peaks, prominences = resolved_peaks(p, prominence)
    result = dict(convention="Complete excess in two periodic watersheds above common current minimum; continuous subcell-parabolic saddle boundaries, piecewise-linear fractional integration",
        common_current_floor=floor, full_active_amount=float(C*np.mean(p)),
        full_excess_amount=float(C*(np.mean(p)-floor)), contrast=contrast,
        resolved_peak_prominence=prominence, resolved_peaks=len(peaks),
        peak_cells=peaks.tolist(), peak_values=p[peaks].tolist(),
        peak_prominences=prominences.tolist(),
        global_vector=geometry.activity_vector(p, q) if detailed else None,
        whole_field_diagnostic=whole_field_check(p, C, incumbent_fraction, cue_fraction) if detailed else None,
        threshold_topology={str(f): len(geometry.components(p >= floor+f*contrast))
                            if contrast > 1e-12 else 0 for f in (.25, .5, .75)} if detailed else None,
        two_region_available=False, regional_dominance=None)
    if contrast < .05 or len(peaks) == 0:
        result.update(state_class="homogeneous_or_unresolved_low_contrast")
        return result
    old_distance = np.abs(geometry.arc(q[peaks], incumbent_fraction, 1.))
    cue_distance = np.abs(geometry.arc(q[peaks], cue_fraction, 1.))
    if len(peaks) == 1:
        if abs(old_distance[0]-cue_distance[0]) <= 1e-10:
            identity = "single_region_between_sites"
        else:
            identity = "single_region_cue_associated" if cue_distance[0] < old_distance[0] else "single_region_incumbent_associated"
        result.update(state_class=identity, single_peak_fraction=float(q[peaks[0]]),
            distance_to_incumbent=float(C*old_distance[0]), distance_to_cue=float(C*cue_distance[0]),
            note="One-region relocation/merger is not a two-region excess handover")
        if identity in ("single_region_cue_associated", "single_region_incumbent_associated"):
            incumbent = result["full_excess_amount"] if identity == "single_region_incumbent_associated" else 0.
            challenger = result["full_excess_amount"]-incumbent
            result.update(incumbent_excess=incumbent, challenger_excess=challenger,
                allocated_region_dominance=-1. if incumbent > 0 else 1.,
                allocation_note="One complete region assigned by spatial proximity only; a change of this label is a relocation/merger event, not a continuous two-region handover")
        return result
    if len(peaks) != 2:
        result.update(state_class="multiple_resolved_regions")
        return result
    old, challenger = int(np.argmin(old_distance)), int(np.argmin(cue_distance))
    if old == challenger or abs(old_distance[0]-old_distance[1]) <= 1e-10 or abs(cue_distance[0]-cue_distance[1]) <= 1e-10:
        result.update(state_class="two_peaks_spatial_association_ambiguous")
        return result
    a, b = int(peaks[0]), int(peaks[1])
    first, second = valley(p, a, b), valley(p, b, a)
    if first["ambiguous_nonunique_active_saddle"] or second["ambiguous_nonunique_active_saddle"]:
        result.update(state_class="two_peaks_with_ambiguous_active_saddle", valley_descriptors=[first, second])
        return result
    dx = C/n
    # Peak a belongs to the basin from the outer-route valley to inner valley.
    v1 = (first["cell_coordinate"]+.5)*dx
    v2 = (second["cell_coordinate"]+.5)*dx
    basins = [(v2 % C, (v1-v2) % C), (v1 % C, (v2-v1) % C)]
    descriptors = []
    for index, (left, width) in enumerate(basins):
        amount = excess_integral(p, C, left, width, floor)
        sensitivities = {str(f): excess_integral(p, C, left, width, floor,
            threshold=floor+f*(float(p[peaks[index]])-floor)) for f in (.25, .5, .75)} if detailed else {}
        descriptors.append(dict(peak_cell=int(peaks[index]), peak_fraction=float(q[peaks[index]]),
            basin_left_physical=left, basin_width_physical=width,
            full_basin_excess=amount, local_peak_contour_excess=sensitivities))
    E_I, E_C = descriptors[old]["full_basin_excess"], descriptors[challenger]["full_basin_excess"]
    total = E_I+E_C
    np.testing.assert_allclose(total, result["full_excess_amount"], rtol=1e-12, atol=1e-12)
    if E_I <= 1e-12 or E_C <= 1e-12:
        result.update(state_class="two_peaks_with_unresolved_excess_amount")
        return result
    contour_sensitivity = {}
    for f in ((.25, .5, .75) if detailed else ()):
        old_amount = descriptors[old]["local_peak_contour_excess"][str(f)]
        cue_amount = descriptors[challenger]["local_peak_contour_excess"][str(f)]
        contour_sensitivity[str(f)] = dict(incumbent_excess=old_amount, challenger_excess=cue_amount,
            excess_ratio=cue_amount/old_amount if old_amount > 1e-12 else None,
            dominance=(cue_amount-old_amount)/(cue_amount+old_amount) if cue_amount+old_amount > 1e-12 else None)
    result.update(state_class="two_resolved_spatially_associated_regions", two_region_available=True,
        incumbent_excess=E_I, challenger_excess=E_C, challenger_to_incumbent_ratio=E_C/E_I,
        regional_dominance=(E_C-E_I)/total, incumbent_basin=descriptors[old],
        allocated_region_dominance=(E_C-E_I)/total,
        challenger_basin=descriptors[challenger], valley_descriptors=[first, second],
        amounts_partition_complete_excess=True, contour_sensitivity=contour_sensitivity,
        peak_center_distance=float(C*abs(geometry.arc(q[peaks[old]], q[peaks[challenger]], 1.))),
        contour_measured_gaps={str(f): geometry.measured_regions(q*C, C, p,
            incumbent_fraction*C, cue_fraction*C, f) for f in (.25, .5, .75)} if detailed else {})
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--circumference", type=float, required=True)
    parser.add_argument("--incumbent-fraction", type=float, required=True)
    parser.add_argument("--cue-fraction", type=float, required=True)
    parser.add_argument("--cells", type=int, required=True)
    args = parser.parse_args()
    with np.load(args.archive) as z:
        state = z["states"][-1]
    print(json.dumps(readout(state[:args.cells], args.circumference,
                            args.incumbent_fraction, args.cue_fraction), indent=2))


if __name__ == "__main__":
    main()
