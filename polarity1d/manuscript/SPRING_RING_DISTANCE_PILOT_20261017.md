# Challenger separation changes the transient coexistence of active sites in the Spring ring analogue

The first Spring distance pilot supports a message about **near-site bridging
versus separated competitor growth and old-site extinction**. The nearer cue
forms an active bridge to the incumbent; the two farther cues grow separate
active regions while the incumbent peak decays near its original position.
The target-region maximum exceeds activity one at the first saved sample,
t = 0.1, in all three cases. Under this single long cue, every challenged
position ultimately retains one basal cap centred at the target after withdrawal.

![Spring circular-distance pilot: full active-Rac fields, raw original-core and target-region activity, and evolving circumference for three challenger positions.](figures/length_mechanisms_20261008/dev/spring_ring_distance_pilot_20261017/spring_ring_distance_pilot.png)

[Vector PDF](figures/length_mechanisms_20261008/dev/spring_ring_distance_pilot_20261017/supporting/vectors/spring_ring_distance_pilot.pdf)
· [Caption](figures/length_mechanisms_20261008/dev/spring_ring_distance_pilot_20261017/supporting/metadata/caption.txt)
· [Raw region time courses](figures/length_mechanisms_20261008/dev/spring_ring_distance_pilot_20261017/supporting/data/fixed_region_timecourses.csv)
· [Numerical summary](data/runs/spring_ring_distance_pilot_20261017/results/summary.json)
· [Saved-state audit](data/runs/spring_ring_distance_pilot_20261017/results/audit.json)
· [Fresh-context scientific review](../../current%20manuscript%20figs/reviews/spring_distance_r1.md)

The independent review accepts this bounded competition result. This revision
distinguishes near-site bridging from old-site extinction at the farther
targets, separates threshold connectivity from resolved maxima, and reports
activation crossings as sampling brackets. The review recommends using early
spatial histories and selected raw profiles for a main-text candidate, with
the complete continuation and mechanical checks retained here as support.

## Message and displayed evidence

The same exact incumbent is challenged at initial material offsets 0.3, 0.4 and
0.5 of the circumference. Its actual initial circumference is 2.04418 in native
Spring units, giving shortest initial arc distances 0.6133, 0.8177 and 1.0221.
Every cue has peak 4, total instantaneous physical support 0.2 and duration 1000,
followed by 13000 units without extra cue. All initial cue supports lie outside
the original half-contrast core; no state is moved or modified between trials.

| Initial offset, d/C | Initial target mean | Target mean at t = 0.1 | Target mean at t = 0.2 | Last sample with two half-contrast components | Last sample with two resolved maxima |
|---|---:|---:|---:|---:|---:|
| 0.3 | 0.7799 | 0.9732 | 1.1399 | 0.1 | 0.8 |
| 0.4 | 0.7119 | 0.8974 | 1.0597 | 6.6 | 7.7 |
| 0.5 | 0.7082 | 0.8932 | 1.0548 | 6.9 | 8.2 |

The mean activity-one crossing is bracketed within (0.1, 0.2] in every case;
these samples do not establish identical activation kinetics. The target-region
maxima are already 1.0946, 1.0266 and 1.0233 at t = 0.1. Activity one is a
descriptive raw concentration readout, not a takeover or stable-front gate.

Spatial connectivity and peak counts answer different questions. At offset
0.3, the regions above the frozen original half-contrast level are connected
by sampled t = 0.2, while two local maxima remain resolved through t = 0.8.
At offsets 0.4 and 0.5, two above-half-contrast components remain through
t = 6.6 and 6.9. The second component disappears as the incumbent decays;
the incumbent maximum remains near its original position rather than moving
into the target peak. Peak counts use circular peak prominence at least 0.05
native concentration, the same resolved-peak convention as the preparation
pilot. All quoted endpoints are sampled observations, not exact event times.

At t = 5, the more separated cases retain incumbent local maxima
of 1.0690 and 1.0812, with prominence 0.3235 and 0.3711. These are distinct
spatial maxima alongside the target peak, not a first-moment interpretation of
a single moving cap. In the nearer case only the target peak remains resolved
at this time, and the above-half-contrast activity is already connected.

All three final cue-free states have the initial basal profile translated to
their respective target positions: the maximum full chemical-field residual is
below 9.5 × 10⁻¹¹, and final circumference differs from its initial value by
less than 7.1 × 10⁻¹⁴. These final comparisons use complete active and inactive
profiles. The final 1000-unit full-field and circumference stationarity checks
also pass on the frozen preparation scales. These trajectories therefore
distinguish different transient routes to the same translated basal state;
they do not identify different final fates at the tested cue duration.

## Why this figure conveys that message

The full spatial histories show whether the activity remains connected,
contains simultaneous distinct peaks, or moves to the target. A direction score
alone cannot make those distinctions. All heatmaps share the same raw Rac
scale; their time axes are linear through 10 native units and logarithmic
thereafter, making both the early two-site transient and prolonged cue-free
retention visible. The white dotted line marks withdrawal at 1000. The fields
are cyclically shifted for display only, placing the original incumbent near
material position 0.25; the two displayed edges remain connected.

The middle row exposes the raw activity sampled at the original core and the
target, without amplitude normalization. These are **fixed material regions**.
Their late values must be interpreted alongside the spatial fields: a broad
cap translated to a nearby target overlaps more of the original core than one
translated to the opposite side. Higher final original-core activity therefore
does not establish survival of the incumbent as a second active site.

The circumference trace makes the retained mechanical feedback explicit. All
trials start with identical circumference, but their subsequent extension can
differ. A test that silently clamped circumference would remove part of this
Spring analogue’s mechanism. Its observed circumference range is approximately
2.04157–2.04432 over the present trajectories.

This is a development evidence figure. A main-text candidate could retain the
early spatial histories and one activity comparison, with the full continuation
and circumference checks linked as support; the final manuscript layout is not
selected here.

## Scope and limitations

This ring uses the [reviewed circular global-extension analogue](SPRING_RING_CLOSURE_20261017.md):
the native local reaction, diffusion, force switch, extension-feedback and
dilution coefficients are unchanged, while force is averaged around a uniformly
deforming periodic domain. This spatial averaging rule is an explicit new
closure. It is not the original 2D continuum mechanics, a local cortical-flow
model or a protrusive-gap model.

Moving the cue changes its adjacency to the incumbent’s existing flank as well
as the two available arc lengths. Although the initial cue and half-contrast
core do not overlap, the near target begins with higher Rac (table above).
Consequently these three trajectories do not isolate a pure travel-time law
from an identical local chemical background. The same global starting state is
preserved, as intended by the circular-location experiment.

The 1000-unit cue is long compared with the observed early competition. It was
specified before running these trials and was not optimized. Its shared final
outcome does not establish separation-independent pulse thresholds, responses
to brief cues or spontaneous migration persistence. There is one resting size,
one grid and one cue dose. No noise-based formation probability, force measured
in cells or universal relation between separation and coexistence is inferred.

## Raw source and verification

Fresh frozen sources and all full-field trajectories are under
[spring_ring_distance_pilot_20261017](data/runs/spring_ring_distance_pilot_20261017/),
SLURM array **21083218**, four independent one-CPU workers with QoS `small` and
3 GB each, without a concurrency throttle. No parameter search, repeated
preparation or additional integrations were performed.

The input is an exact copy of the [qualified Spring checkpoint](data/runs/spring_ring_preparation_pilot_20261017/workers/qualified_basal.npz),
SHA-256 `ae07e97a6038c2ddf02cc76892889f0ae3d57000b239c396222072773716327f`.
Each worker stores `maintained_cue.npz` and `recovery.npz`: `states` has active
Rac in columns 0–199, inactive Rac in columns 200–399, and circumference in
column 400. Physical arc coordinates are material coordinates multiplied by
each sample’s circumference. The cue centre stays at a fixed material offset;
its support follows the instantaneous physical circumference.

The BDF solver uses a full dense numerical Jacobian, preserving the global
dependencies introduced by mean force and dilution. Output spacing is 0.1 for
the first 10 units, 1 to 100 and 5 thereafter, restarted at cue withdrawal.
The audit checks immutable source and input hashes, exact initial-state and
withdrawal continuity, full sample grids, finite nonnegative fields, physical
amount and final stationarity. All four workers pass, with empty stderr and
maximum sampled physical-amount error below 2.5 × 10⁻¹²; concentrations were
not clipped or renormalized.

A reflected −0.4 offset checks numerical symmetry. After reflection its complete
chemical fields agree with the +0.4 trajectory within 1.9 × 10⁻¹⁰ and its
circumference within 5.0 × 10⁻¹⁴. This check does not establish spatial
convergence. The new trajectories are scientifically distinct from evidence
already contained in the original length-mechanism PDF.
