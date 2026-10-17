# Phase 1 simulation campaign: experimental design

This document defines what each of the 688 tasks in
`phase1_production_v1` actually simulated. The campaign is not one experiment:
it contains three assays applied to the four Jilkine-reviewed models and four
assays applied to the De Belly mechanochemical model.

![Four stimulation protocols](results/campaigns/phase1_production_v1/campaign_protocols.png)

## Scientific scope

The wave-pinning, Otsuji, Goryachev, and LEGI implementations are frozen as
equation-level reconstructions with the assumptions in `MODEL_ASSUMPTIONS.md`.
The results below describe those implementations under the declared campaign
protocols. They are not claimed to reproduce length thresholds or numerical
branches in the Jilkine publication. Exact reproduction was abandoned because
the source does not document a reproducible symmetry-breaking mechanism and
omits other consequential numerical details.

The optimized finite-volume De Belly port is the accepted Phase 1 screening
implementation. A native-Gridap trajectory cross-check was deliberately left
out of scope.

## Campaign inventory

There are 680 production tasks and eight numerical-refinement tasks:

| Model group | Assay | Production tasks |
|---|---|---:|
| Four Jilkine-reviewed models | two-site challenge | 320 |
| Four Jilkine-reviewed models | gradient reversal | 80 |
| Four Jilkine-reviewed models | alternating-gradient tracking | 64 |
| De Belly, six coupling variants | single-site propagation | 48 |
| De Belly, six coupling variants | two-site challenge | 48 |
| De Belly, six coupling variants | gradient reversal | 48 |
| De Belly, six coupling variants | alternating-gradient tracking | 72 |

All tasks are listed explicitly in
`campaigns/phase1_production_v1.json`. The manifest fixes the model parameters,
protocol, length, mesh, integration time, replicate/seed, and phenotype
settings for every task.

## Jilkine-model experiments

The four models are wave-pinning, Otsuji, Goryachev, and LEGI. Their reported
observables are active GTPase, active GTPase, active Cdc42, and LEGI response,
respectively.

### J1 — two transient end stimuli

Question: after both ends are stimulated, does the model sustain two response
domains, pass through two domains and settle to one, or never establish two
domains?

- Lengths: `L = 1, 3, ..., 39 µm`.
- Input: reflected cosine caps over the outer `0.1L` at both ends, maximum
  activation rate `0.2 s⁻¹`, applied from `t=0` to `40 s` with a half-cosine
  taper during the last `20 s`.
- Symmetry breaking: the left cue is multiplied by `1.01` and the right by
  `0.99`, a 2% left/right difference.
- Initial conditions: wave-pinning, Otsuji, and Goryachev receive conservative
  zero-mean perturbations with concentration-scale amplitude `0.01`, using
  seeds 1–5. LEGI is deterministic and has one run per length.
- Run end: `max(2000, 500(L/10)²) s`.
- Primary readout: competition classification and time from sustained
  two-front establishment to sustained one-front resolution.

The cue width is a fixed fraction of length, not a fixed physical width. Thus,
this scan changes both cell length and the physical width of each stimulus.

### J2 — gradient reversal

Question: if a linear input first polarizes the model toward one end and is
then reversed, does a sustained front form at the new end, and how long does it
take?

- Lengths: `L = 1, 3, ..., 39 µm`; one deterministic run per length.
- Input before reversal: a linear gradient with maximum activation rate
  `0.1 s⁻¹` at the left end and zero at the right.
- Conditioning time: `max(200, 100(L/10)²) s`.
- Reversal: linear cross-fade to the reflected gradient over
  `min(20, 0.1 × conditioning time) s`.
- Post-reversal observation: `max(400, 200(L/10)²) s`.
- Primary readout: whether a qualifying front existed immediately before
  reversal, whether a sustained opposite front subsequently formed, and its
  latency.

A run without a polarized, directionally oriented state before the cue changed
is reported as `not_conditioned`; it is not counted as a failed reversal.
Multiple detected local maxima do not by themselves invalidate an old side.

### J3 — alternating-gradient tracking

Question: how strongly and how quickly does left/right model orientation follow
a periodically reversing gradient?

- Lengths: `5, 15, 25, 39 µm`.
- Periods: `50, 100, 200, 400 s`.
- Input: a nonnegative linear gradient with fixed peak `0.1 s⁻¹`; its direction
  varies sinusoidally.
- Duration: six periods; the first period is excluded from scoring.
- Replication: one deterministic run per length/period cell.
- Primary readouts: normalized orientation amplitude, phase lag, correlation
  with cue direction, and a declared pass/fail rule.

The pass rule requires amplitude at least `0.05`, absolute lag no more than
`0.25` periods, and positive cue/response correlation. The continuous amplitude
and lag should be read alongside the thresholded label.

## De Belly experiments

The De Belly observable used for common polarity classifications is Rac. Every
condition is one deterministic run; the campaign does not estimate biological
or numerical variability for this model.

Lengths are `12, 16, 20, 25.1984, 32, 40, 50, 64 µm`. The `25.1984-µm` point
is `2^(1/3) × 20 µm`, included as a 1.26× linear-size anchor motivated by whole
genome doubling. The mesh uses `max(60, round(L/0.4))` cells and the default
fixed step is `0.5 s`.

### Coupling perturbations

The screen crosses every De Belly assay with these six parameter states:

| Manifest name | Rac→mechanics→Rho arm | Rho→mechanics→Rac arm |
|---|---:|---:|
| `intact` | 1.0 | 1.0 |
| `rac_to_rho_half` | 0.5 | 1.0 |
| `rho_to_rac_half` | 1.0 | 0.5 |
| `rac_to_rho_off` | 0.0 | 1.0 |
| `rho_to_rac_off` | 1.0 | 0.0 |
| `both_off` | 0.0 | 0.0 |

The Rac→Rho multiplier scales Rac-driven protrusion and tension-dependent Rho
activation. The Rho→Rac multiplier scales Rho active stress and MCA-dependent
Rac activation. Basal Rac/Rho mutual antagonism remains present when either or
both mechanochemical arms are disabled.

### D1 — single-site propagation

Question: after local Rac activation at the left edge, what response reaches
the unstimulated right side, and how is it changed by the two mechanochemical
arms?

- Input: a left-centered Gaussian Rac-activation pulse, amplitude `1`, standard
  deviation `0.2L`, active from `t=50` to `150 s`.
- Run end: `700 s`.
- Proximal and distal regions: the leftmost and rightmost 20% of the domain.
- Baseline: the last saved frame before `t=50 s` (`t=48 s`).
- Primary readouts: signed peak changes in proximal/distal Rac and Rho,
  distal/proximal absolute response ratios, peak times, and maximum absolute
  tension, velocity, and displacement.

Because the Gaussian width is `0.2L`, physical stimulus width grows with cell
length. The size curve cannot distinguish a pure domain-size effect from a
cue-width effect.

### D2 — two-site challenge

Question: after both ends receive Rac activation, is there a sustained
two-front phase, is it followed by a sustained one-front phase, and what is the
final Rac endpoint?

- Input: Gaussians centered on both ends, standard deviation `0.08L`, amplitude
  `1`, active from `t=50` to `150 s`.
- Symmetry breaking: left and right inputs are multiplied by `1.01` and `0.99`.
- Run end: `700 s`.
- Primary readouts: competition classification, resolution latency, and final
  Rac endpoint class.

### D3 — gradient reversal

Question: after a linear Rac-activation gradient is reversed, does Rac
orientation form at the new end?

- Input peak: `1` at one end and zero at the other.
- Reversal begins at `t=250 s` and cross-fades over `25 s`.
- Run end: `750 s`.
- Primary readouts: side-orientation conditioning, sustained side-switch status,
  and latency. The one-detected-maximum version is retained as a stricter
  morphology diagnostic.

### D4 — alternating-gradient tracking

Question: does Rac orientation follow a periodically alternating activation
gradient, and do the coupling perturbations change amplitude or lag?

- Lengths: `20, 25.1984, 40, 64 µm`.
- Periods: `100, 200, 400 s`.
- Input: sinusoidally alternating gradient with peak `1`.
- Duration and scoring: six periods, first period excluded; the same amplitude,
  lag, and correlation definitions as J3.

## Common phenotype definitions

For a saved observable profile, spatial contrast is
`(maximum - minimum) / (|maximum| + |minimum|)`. A frame is polarized when this
is at least `0.1`. Peaks are counted with prominence at least 8% of the profile
range and separation at least 12% of the grid, including eligible boundary
peaks.

Each frame is classified as:

- `unpolarized`: contrast below 0.1;
- `single_front`: polarized with one detected domain;
- `multifront`: polarized with at least two detected domains;
- `polarized_no_resolved_peak`: polarized but no peak passes the domain rule.

Competition is defined as an event sequence, not just a final frame:

1. `not_established`: no multifront state persists for the assay's declared
   sustained duration;
2. `resolved`: after sustained multifront establishment, a single-front state
   later persists for the sustained duration;
3. `right_censored`: multifront establishment occurred, but no sustained
   single-front resolution was observed before the run ended.

Importantly, `right_censored` does not imply that the final frame is
multifront. A run can establish two fronts, never enter a sustained one-front
state, and later become unpolarized. The endpoint class is therefore always
reported separately.

Orientation is the normalized difference between the mean observable in the
rightmost and leftmost 20% of the domain. Positive values point right and
negative values point left.

For reversal, the primary precondition is evaluated at the last saved frame
before reversal begins. A run is orientation-conditioned when that frame is
polarized and `|orientation| ≥ 0.05`; the number and locations of detected
local maxima are reported but do not determine whether an old side exists. A
reversal succeeds at the first post-reversal frame from which the run remains
polarized with orientation at least 0.05 toward the opposite side for the
declared sustained duration. Latency is measured from the start of the cue
cross-fade. A conditioned run without such an event before run end is
`right_censored`.

The parallel strict fields (`reversal_conditioned`, `reversal_success`, and
`reversal_latency_s`) additionally require no more than one detected maximum
both before and after the switch. They are useful when the question is
specifically one-domain morphology, but they are not the primary measure of
which side Rac occupies. The primary side-switch fields are prefixed
`orientation_reversal_`.

## Interpretation limits

- Only the wave-pinning, Otsuji, and Goryachev two-site scans have seeded
  replicates. LEGI and all other assay cells are single deterministic runs.
- The seed scan probes sensitivity to a declared small initial perturbation;
  five seeds are not a biological sample size.
- Fractional stimulus widths confound physical cue width with cell length in
  both localized size scans.
- Thresholded endpoint and tracking labels are convenience summaries. The
  continuous trajectories, amplitudes, lags, and contrasts remain the primary
  evidence.
- The Jilkine-model results are properties of the current implementations and
  protocols, not a reproduction of undocumented published branch selection.
- The De Belly results are properties of the accepted 1D screening port, not a
  native Gridap validation.

## Reproduce

Use the existing `scaling_sandbox` Conda environment:

```bash
conda run -n scaling_sandbox python polarity1d/scripts/create_phase1_campaign.py

mkdir -p polarity1d/results/campaigns/phase1_production_v1/logs
sbatch --array=0-687 polarity1d/workflow/run_phase1_campaign.slurm

conda run -n scaling_sandbox python polarity1d/scripts/collect_phase1_campaign.py \
  polarity1d/campaigns/phase1_production_v1.json \
  polarity1d/results/campaigns/phase1_production_v1

conda run -n scaling_sandbox python polarity1d/scripts/analyze_phase1_campaign.py \
  polarity1d/results/campaigns/phase1_production_v1
```

Do not add an array concurrency throttle unless scheduler, filesystem, or
resource constraints justify one. An existing task result is reused only when
its configuration digest matches the manifest.
