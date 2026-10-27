"""Imaging-style periodic-region tracking for the final Figure 4 brief.

Observed regions are connected components above a current contrast contour.
Their amounts cover complete native-peak basins, including activity below
that contour. Identity follows maximum overlap; splitting and merging remain
explicit observed topology events.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import linear_sum_assignment
from scipy.signal import find_peaks


def components(mask):
    mask = np.asarray(mask, bool)
    if mask.ndim != 1 or not len(mask):
        raise ValueError("A nonempty periodic one-dimensional mask is required")
    if not mask.any():
        return []
    if mask.all():
        return [np.arange(len(mask))]
    starts = np.flatnonzero(mask & ~np.roll(mask, 1))
    result = []
    for start in starts:
        cells = [int(start)]
        while mask[(cells[-1]+1) % len(mask)]:
            cells.append((cells[-1]+1) % len(mask))
        result.append(np.asarray(cells))
    return result


def _complete_peak_basins(profile, circumference):
    """Partition complete excess at continuous subcell native-peak saddles.

    The prominence tolerance removes floating-point ripples only; it is not
    the imaging contour. Thus an incumbent below the imaging threshold keeps
    its measured amount while its native peak remains resolved. Active flat
    or disjoint equal saddles make the associated amount explicitly ambiguous.
    """
    p, C = np.asarray(profile), float(circumference)
    n, floor, contrast = len(p), float(p.min()), float(np.ptp(p))
    epsilon = max(1e-12, 1e-8*contrast)
    peaks, properties = find_peaks(np.tile(p, 3), prominence=epsilon)
    peaks = peaks[(peaks >= n) & (peaks < 2*n)]-n
    if not len(peaks):
        return []
    total = C*float(np.mean(p)-floor)
    if len(peaks) == 1:
        return [dict(peak_cell=int(peaks[0]), mask=np.ones(n, bool), excess=total,
                     amount_available=True, boundary_ambiguous=False)]
    boundaries, ambiguous = [], []
    for a, b in zip(peaks, np.roll(peaks, -1)):
        route = int(a)+np.arange((int(b)-int(a)) % n+1)
        values = p[route % n]
        minimum = float(values.min())
        tied = route[np.isclose(values, minimum, rtol=0., atol=1e-14)]
        contiguous = len(tied) == 1 or np.all(np.diff(tied) == 1)
        active_flat = len(tied) > 1 and minimum > floor+epsilon
        uncertainty = bool((not contiguous or active_flat) and minimum > floor+epsilon)
        if len(tied) > 1:
            location = float((tied[0]+tied[-1])/2) if contiguous else float(tied[0])
        else:
            middle = int(tied[0])
            left, center, right = p[(middle-1) % n], p[middle % n], p[(middle+1) % n]
            curvature = left-2*center+right
            offset = .5*(left-right)/curvature if curvature > 0 else 0.
            location = middle+float(np.clip(offset, -.5, .5))
        boundaries.append((location+.5)*C/n % C)
        ambiguous.append(uncertainty)
    dx = C/n
    prefix = np.r_[0., np.cumsum(.5*dx*(p+np.roll(p, -1)))]

    def primitive(x):
        coordinate = (x-.5*dx)/dx
        index = math.floor(coordinate)
        fraction = coordinate-index
        cycles, k = divmod(index, n)
        return (cycles*prefix[-1]+prefix[k]
                +dx*(p[k]*fraction+.5*(p[(k+1) % n]-p[k])*fraction*fraction))

    q = (np.arange(n)+.5)*dx
    result = []
    for j, peak in enumerate(peaks):
        left, right = boundaries[j-1], boundaries[j]
        width = (right-left) % C
        excess = primitive(left+width)-primitive(left)-floor*width
        uncertain = bool(ambiguous[j-1] or ambiguous[j])
        result.append(dict(peak_cell=int(peak), mask=(q-left) % C < width,
            excess=float(excess), amount_available=not uncertain,
            boundary_ambiguous=uncertain))
    if not np.isclose(sum(r["excess"] for r in result), total, rtol=1e-10, atol=1e-10):
        raise ValueError("Complete native peak basins do not partition the activity excess")
    return result


def frame_regions(profile, circumference, basal_contrast, threshold=.5, smoothing=0.):
    """Reproduce native contour components with the brief's relative qualifier.

    Width uses native contour-cell quadrature. Amount includes the complete
    valley-bounded basins of every native peak inside that component, rather
    than discarding the lower activity outside its imaging contour.
    Smoothing is a periodic Gaussian with a declared physical
    standard deviation, converted using this frame's actual circumference.
    """
    profile = np.asarray(profile, float)
    C = float(circumference)
    if (profile.ndim != 1 or len(profile) < 3 or not np.isfinite(profile).all()
            or not np.isfinite(C) or C <= 0 or not np.isfinite(basal_contrast) or basal_contrast <= 0
            or not 0 < threshold < 1 or not np.isfinite(smoothing) or smoothing < 0):
        raise ValueError("Finite native profile, circumference, basal contrast and contour required")
    p = gaussian_filter1d(profile, smoothing*len(profile)/C, mode="wrap") if smoothing else profile
    floor, contrast = float(p.min()), float(np.ptp(p))
    qualified = contrast > 1e-12 and contrast >= .1*basal_contrast
    regions, atoms = [], _complete_peak_basins(p, C) if qualified else []
    assigned_atoms = set()
    if qualified:
        for cells in components(p >= floor+threshold*contrast):
            mask = np.zeros(len(p), bool)
            mask[cells] = True
            excess = p[cells]-floor
            angles = 2*np.pi*(cells+.5)/len(p)
            moment = np.sum(excess*np.exp(1j*angles))
            selected = [i for i, atom in enumerate(atoms) if mask[atom["peak_cell"]]]
            if not selected:
                raise ValueError("A qualified region has no resolved native peak")
            assigned_atoms.update(selected)
            full_mask = np.logical_or.reduce([atoms[i]["mask"] for i in selected])
            regions.append(dict(mask=mask, cells=cells, full_mask=full_mask,
                excess=sum(atoms[i]["excess"] for i in selected),
                contour_excess=float(excess.sum()*C/len(p)),
                amount_available=all(atoms[i]["amount_available"] for i in selected),
                width=C*len(cells)/len(p), width_fraction=len(cells)/len(p),
                heading=float(np.angle(moment) % (2*np.pi)) if abs(moment) > 1e-12 else None,
                peak=float(p[cells].max())))
    return dict(regions=regions, qualified=qualified, contrast=contrast,
                baseline=floor, smoothed_mean=float(p.mean()),
                latent_atoms=[atom for i, atom in enumerate(atoms) if i not in assigned_atoms],
                native_peak_basins=len(atoms),
                low_contrast=contrast < .1*basal_contrast)


def wrapped_angle(value):
    return float((value+np.pi) % (2*np.pi)-np.pi)


def _dominant(regions):
    if not regions:
        return None
    # An ambiguous amount can change the ordering. It cannot silently lose
    # the dominant identity merely because a numerical saddle was selected.
    if any(not region["amount_available"] for region in regions):
        return None
    order = sorted(regions, key=lambda r:r["excess"], reverse=True)
    if len(order) > 1 and math.isclose(order[0]["excess"], order[1]["excess"],
                                      rel_tol=1e-12, abs_tol=1e-12):
        return None
    return order[0]


def survival_curve(tracks, *, event_clock="observed"):
    """Incident-track Kaplan–Meier with explicit observation-interval sensitivity.

    Tracks already present at analysis start have unknown birth time and are
    excluded from the incident cohort. End-of-window, unresolved low contrast
    and merger identity loss are right censored. Native disappearance times
    have lower/upper observation brackets; report both alongside the nominal
    first-missing minus first-seen clock.
    """
    if event_clock not in ("observed", "lower", "upper"):
        raise ValueError("Unknown survival observation clock")
    sample = []
    for track in tracks:
        if track["left_censored"]:
            continue
        event = not track["right_censored"]
        key = "duration" if event_clock == "observed" or not event else "duration_"+event_clock
        sample.append((float(track[key]), event))
    survival = 1.
    answer = [dict(time=0., survival=1., at_risk=len(sample), events=0, censored=0)]
    for time in sorted({t for t, event in sample}):
        at_risk = sum(t >= time for t, event in sample)
        deaths = sum(t == time and event for t, event in sample)
        censored = sum(t == time and not event for t, event in sample)
        survival *= 1-deaths/at_risk
        answer.append(dict(time=time, survival=survival, at_risk=at_risk,
                           events=deaths, censored=censored))
    return dict(clock=event_clock, incident_tracks=len(sample),
                excluded_unknown_birth=sum(t["left_censored"] for t in tracks), curve=answer)


def track_profiles(time, profiles, circumference, *, basal_contrast, threshold=.5,
                   smoothing=0., turning_intervals=(1, 2, 4, 8)):
    """Track an already-selected analysis window; preserve every censoring reason."""
    time, profiles = np.asarray(time, float), np.asarray(profiles, float)
    if (time.ndim != 1 or len(time) < 2 or profiles.ndim != 2
            or len(profiles) != len(time) or not np.isfinite(time).all()
            or np.any(np.diff(time) <= 0)):
        raise ValueError("Aligned native profiles and strictly increasing time required")
    C = np.broadcast_to(np.asarray(circumference, float), time.shape)
    if not np.isfinite(C).all() or np.any(C <= 0):
        raise ValueError("Actual circumference must be positive in every frame")
    tracks, previous, previous_dominant = {}, [], None
    frames, births, takeover_events, pairwise_events, overlap_events = [], [], [], [], []
    next_id = 0
    for frame, (t, p, length) in enumerate(zip(time, profiles, C)):
        descriptor = frame_regions(p, length, basal_contrast, threshold, smoothing)
        current = descriptor["regions"]
        overlap = np.zeros((len(previous), len(current)), int)
        for i, old in enumerate(previous):
            for j, new in enumerate(current):
                overlap[i, j] = np.count_nonzero(old["mask"] & new["mask"])
        assigned, retained = {}, set()
        if overlap.size:
            ii, jj = linear_sum_assignment(-overlap)
            for i, j in zip(ii, jj):
                if overlap[i, j] > 0:
                    assigned[int(j)] = previous[int(i)]["track_id"]
                    retained.add(previous[int(i)]["track_id"])
            for i, old in enumerate(previous):
                children = np.flatnonzero(overlap[i] > 0)
                if len(children) > 1:
                    overlap_events.append(dict(time=float(t), kind="split",
                        parent_id=old["track_id"], current_component_indices=children.tolist(),
                        note="Observed threshold-component split; biological ancestry is not inferred"))
        for j, region in enumerate(current):
            parents = [previous[i]["track_id"] for i in range(len(previous)) if overlap[i, j] > 0]
            if j in assigned:
                identity = assigned[j]
            else:
                identity = next_id
                next_id += 1
                dominant_id = previous_dominant["track_id"] if previous_dominant is not None else None
                angle = (abs(wrapped_angle(region["heading"]-previous_dominant["heading"]))
                         if region["heading"] is not None and previous_dominant is not None
                         and previous_dominant["heading"] is not None else None)
                tracks[identity] = dict(track_id=identity, birth_observed=float(t),
                    birth_lower=float(time[frame-1]) if frame else None,
                    last_seen=float(t), left_censored=frame == 0, parent_ids=parents,
                    previous_dominant_id=dominant_id, takeover_streak=0, takeover_recorded=False,
                    pairwise_streak=0, pairwise_recorded=False,
                    takeover_unscorable_reason=None)
                if frame:
                    birth = dict(track_id=identity, time=float(t), birth_lower=float(time[frame-1]),
                        parent_ids=parents, previous_dominant_id=dominant_id,
                        distance_from_dominant_radians=angle,
                        distance_from_dominant_physical=angle*length/(2*np.pi) if angle is not None else None,
                        observed_split_from_dominant=dominant_id is not None and dominant_id in parents,
                        no_overlap_with_previous_regions=not parents,
                        opposite_hemisphere_birth=bool(not parents and angle is not None and angle >= np.pi/2))
                    births.append(birth)
            region["track_id"] = identity
            tracks[identity]["last_seen"] = float(t)
            tracks[identity]["last_full_mask"] = region["full_mask"].copy()
            if len(parents) > 1:
                overlap_events.append(dict(time=float(t), kind="merge", parents=parents,
                    retained_track_id=identity,
                    note="Maximum-overlap ancestry; no invented newborn identity at merger"))
        for i, old in enumerate(previous):
            identity = old["track_id"]
            if identity in retained:
                continue
            merged = bool(overlap.shape[1] and np.any(overlap[i] > 0))
            reason = "merger_identity_loss" if merged else "low_contrast_unresolved" if not descriptor["qualified"] else "disappearance"
            record = tracks[identity]
            record.update(end_observed=float(t), end_reason=reason,
                right_censored=reason != "disappearance",
                duration=float((record["last_seen"] if reason != "disappearance" else t)-record["birth_observed"]),
                duration_lower=float(record["last_seen"]-record["birth_observed"]),
                duration_upper=float(t-(record["birth_lower"] if record["birth_lower"] is not None else record["birth_observed"])))
        amounts = {region["track_id"]:region["excess"] for region in current if region["amount_available"]}
        amount_source = {identity:"qualified_complete_region" for identity in amounts}
        references = {tracks[region["track_id"]]["previous_dominant_id"] for region in current
            if not tracks[region["track_id"]]["left_censored"]
            and not tracks[region["track_id"]]["takeover_recorded"]}
        hidden = [identity for identity in references if identity is not None and identity not in amounts
            and tracks[identity].get("end_reason") == "disappearance"]
        assigned_hidden = {identity:[] for identity in hidden}
        ambiguous_hidden = set()
        for atom in descriptor["latent_atoms"]:
            scores = [(int(np.count_nonzero(atom["mask"] & tracks[identity]["last_full_mask"])), identity)
                      for identity in hidden if tracks[identity]["last_full_mask"][atom["peak_cell"]]]
            if not scores:
                continue
            maximum = max(score for score, identity in scores)
            selected = [identity for score, identity in scores if score == maximum and score > 0]
            if len(selected) != 1:
                ambiguous_hidden.update(selected)
            else:
                assigned_hidden[selected[0]].append(atom)
        for identity, atoms in assigned_hidden.items():
            if atoms and identity not in ambiguous_hidden and all(a["amount_available"] for a in atoms):
                amounts[identity] = sum(a["excess"] for a in atoms)
                amount_source[identity] = "resolved_native_peak_basins_below_imaging_contour"
                tracks[identity]["last_full_mask"] = np.logical_or.reduce([a["mask"] for a in atoms])
        dominant = _dominant(current)
        for region in current:
            identity = region["track_id"]
            record = tracks[identity]
            old_id = record["previous_dominant_id"]
            if record["left_censored"] or old_id is None or record["takeover_recorded"]:
                continue
            old = tracks[old_id]
            if old_id not in amounts or identity not in amounts:
                record.update(takeover_streak=0, pairwise_streak=0,
                    takeover_unscorable_reason="previous_dominant_or_newborn_complete_amount_unresolved")
                continue
            old_amount, new_amount = amounts[old_id], amounts[identity]
            dominance = (new_amount-old_amount)/(new_amount+old_amount) if new_amount+old_amount > 1e-12 else None
            pairwise_success = dominance is not None and dominance > .5
            global_success = (pairwise_success and dominant is not None
                              and dominant["track_id"] == identity)
            record["pairwise_streak"] = record["pairwise_streak"]+1 if pairwise_success else 0
            record["takeover_streak"] = record["takeover_streak"]+1 if global_success else 0
            if dominant is None:
                record["takeover_unscorable_reason"] = "dominant_ordering_among_qualified_regions_unresolved"
            event = dict(time=float(t), first_qualifying_frame_time=float(time[frame-1]),
                    newborn_track_id=identity, previous_dominant_id=old_id,
                    dominance=dominance, required_consecutive_frames=2,
                    previous_dominant_amount_source=amount_source[old_id],
                    previous_dominant_excess=old_amount, newborn_excess=new_amount)
            if record["pairwise_streak"] == 2 and not record["pairwise_recorded"]:
                record["pairwise_recorded"] = True
                pairwise_events.append(dict(event, globally_dominant_required=False))
            if record["takeover_streak"] == 2:
                record["takeover_recorded"] = True
                takeover_events.append(dict(event, globally_dominant_required=True))
        frames.append(dict(time=float(t), circumference=float(length), contrast=descriptor["contrast"],
            low_contrast=descriptor["low_contrast"], qualified=descriptor["qualified"], region_count=len(current),
            dominant_track_id=dominant["track_id"] if dominant else None,
            dominant_heading=dominant["heading"] if dominant else None,
            native_peak_basins=descriptor["native_peak_basins"],
            latent_peak_basins=len(descriptor["latent_atoms"]),
            regions=[{k:v for k, v in region.items() if k not in ("mask", "full_mask", "cells")} for region in current]))
        previous, previous_dominant = current, dominant
    for region in previous:
        record = tracks[region["track_id"]]
        duration = float(time[-1]-record["birth_observed"])
        record.update(end_observed=float(time[-1]), end_reason="analysis_window_end",
            right_censored=True, duration=duration, duration_lower=duration,
            duration_upper=float(time[-1]-(record["birth_lower"] if record["birth_lower"] is not None else record["birth_observed"])))
    dt = np.diff(time)
    exposure = float(np.sum(.5*(C[:-1]+C[1:])*dt))
    count = np.asarray([f["region_count"] for f in frames])
    histogram = {str(int(n)):float(np.sum(dt[count[:-1] == n])/dt.sum()) for n in np.unique(count)}
    width_numerator = sum(sum(r["width"] for r in f["regions"])*duration for f, duration in zip(frames[:-1], dt))
    width_fraction_numerator = sum(sum(r["width_fraction"] for r in f["regions"])*duration for f, duration in zip(frames[:-1], dt))
    region_time = float(np.sum(count[:-1]*dt))
    turning = []
    for interval in turning_intervals:
        interval = int(interval)
        if interval < 1:
            raise ValueError("Turning intervals must be positive frame counts")
        angles = []
        if interval < len(time):
            for first, last in zip(frames[:-interval], frames[interval:]):
                if first["dominant_heading"] is not None and last["dominant_heading"] is not None:
                    angles.append(wrapped_angle(last["dominant_heading"]-first["dominant_heading"]))
        turning.append(dict(frame_interval=interval, usable_pairs=len(angles),
            possible_pairs=max(0, len(time)-interval), angles_radians=angles,
            mean_absolute_angle=float(np.mean(np.abs(angles))) if angles else None))
    all_tracks = [{k:v for k, v in track.items() if k != "last_full_mask"} for track in tracks.values()]
    return dict(protocol=dict(region_contour_fraction=threshold,
        qualification="Current contrast >=0.1 times exact basal contrast; all complete contour components",
        amount="Complete native-peak watershed basins above current minimum; every resolved peak inside a qualified contour belongs to its complete region. Subcell parabolic saddle, piecewise-linear integral.",
        low_peak_amount="An incumbent below the imaging contour retains its resolved native-peak basin amount by spatial overlap; absent/unresolved/merged identity is never assigned zero.",
        overlap="Maximum native material-cell overlap, one-to-one; splitting/merging explicit",
        smoothing_physical_gaussian_sigma=smoothing,
        takeover="Newborn is dominant among qualified observed regions AND D=(E_new−E_previous_dominant)/(E_new+E_previous_dominant)>0.5 in each of two consecutive observed frames. Pairwise-only events are separate.",
        opposite_birth="No previous overlap and at least pi/2 from previous dominant; retain continuous distance distribution",
        lifetime="Incident tracks only; end-of-window/merger/low-contrast identity loss right censored; observed death interval sensitivity retained"),
        frames=frames, tracks=all_tracks, births=births, takeovers=takeover_events,
        pairwise_takeovers=pairwise_events,
        overlap_events=overlap_events, turning=turning,
        survival={clock:survival_curve(all_tracks, event_clock=clock) for clock in ("observed", "lower", "upper")},
        summary=dict(analysis_duration=float(time[-1]-time[0]), perimeter_time_exposure=exposure,
            observed_births=len(births), births_per_time_per_perimeter=len(births)/exposure,
            takeovers=len(takeover_events), pairwise_takeovers=len(pairwise_events),
            region_count_mean=float(np.sum(count[:-1]*dt)/dt.sum()),
            region_count_time_distribution=histogram,
            low_contrast_time_fraction=float(np.sum(dt[[f["low_contrast"] for f in frames[:-1]]])/dt.sum()),
            mean_region_width=width_numerator/region_time if region_time else None,
            mean_region_width_fraction=width_fraction_numerator/region_time if region_time else None,
            incident_tracks=sum(not t["left_censored"] for t in all_tracks),
            initial_unknown_birth_tracks=sum(t["left_censored"] for t in all_tracks)))


def imaging_variants(time, profiles, circumference, *, basal_contrast, kernel_length,
                     basal_profile, basal_circumference,
                     thresholds=(.3, .5, .7), downsampling=(1, 2, 4, 8)):
    """Raw and fixed-physical-kernel Gaussian imaging views at every cadence."""
    time, profiles = np.asarray(time), np.asarray(profiles)
    C = np.broadcast_to(np.asarray(circumference), time.shape)
    basal = np.asarray(basal_profile, float)
    if (basal.ndim != 1 or len(basal) != profiles.shape[1]
            or not math.isclose(float(np.ptp(basal)), basal_contrast, rel_tol=1e-9, abs_tol=1e-12)
            or not np.isfinite(basal_circumference) or basal_circumference <= 0):
        raise ValueError("Imaging sensitivities require the exact native basal profile and physical circumference")
    optical_basal_contrast = float(np.ptp(gaussian_filter1d(
        basal, kernel_length*len(basal)/basal_circumference, mode="wrap")))
    for threshold in thresholds:
        for factor in downsampling:
            if int(factor) != factor or factor < 1:
                raise ValueError("Cadence downsampling must be a positive integer")
            indices = np.unique(np.r_[np.arange(0, len(time), int(factor)), len(time)-1])
            for smooth in (False, True):
                yield dict(threshold=threshold, downsample_factor=int(factor), smoothed=smooth,
                    observed_basal_contrast=optical_basal_contrast if smooth else basal_contrast,
                    result=track_profiles(time[indices], profiles[indices], C[indices],
                        basal_contrast=optical_basal_contrast if smooth else basal_contrast, threshold=threshold,
                        smoothing=kernel_length if smooth else 0.))
