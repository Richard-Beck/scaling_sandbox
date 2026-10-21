# Figure 4 E1–E4 analysis adapter

The implementation is
`polarity1d/manuscript/code/mechanisms/signature_figure4_stage_analysis_20261020.py`.
It reuses the tested A1–A4/A7 measurements and direct independent audits from
the E0 checkpoint. It does not integrate native models. Root owns A5/A6 and
the final integrated hypothesis decisions; stage-local decision rows are not
unconditional H1/H2/H3 conclusions.

The corrected E0 authority is
`polarity1d/manuscript/data/runs/signature_figure4_e0_corrected_basals_20261020`.
The superseded compact-profile version is rejected. E1/E2 analysis preparation
requires an explicit delivered-interim release with hashes of its concrete
artifacts. The adapter verifies the corrected basal inventory and all 12
native-initial-profile identity checks.

## Inputs and immutable preparation

Each completed native campaign supplies `manifest.json`,
`results/status.json` and `results/case_census.json`. Missing or failed native
work is an error. Normalized rows retain exact scalar/profile/status/basal
paths and hashes, actual initial circumference, exact basal contrast,
observation horizon and both sampling cadences.

Run `prepare` with `--stage E1`, `E2`, `E3` or `E4`, `--e0` pointing to the
corrected authority, one or more completed `--native-runs`, a fresh `--out`
beneath manuscript `data/runs`, and `--gate` pointing to the root's delivered
E0 release. E1 also requires
`--forced-rule docs/tasks/briefs/20261020_figure4_e1_forced_rule.json`.
Include the relevant E2 native run when preparing E1: it supplies the missing
Goryachev .02 sigma-one small/mid companions. These are data dependencies,
not substituted or interpolated conditions.

The complete E1 analysis requires both disjoint native roots (792 ready and
72 later sigma outcomes) plus its sigma-one companions. The canonical
reviewed selection defines all 24 model–size–input conditions; a ready-only
22-condition analysis cannot pass as complete. E2 requires all 36 new
conditions, E3 all 40 intermediate-size additions, and E4 the literal 52
trigger conditions at the declared fourfold horizon. Every native manifest
job must be represented once in its completed census with identical labels
and observation schedule. Every new condition requires twelve outcomes.

The producer freezes all five executable analysis modules, its inventory,
release digest, measurement definitions and the reviewed E1 rule. Run its
frozen `verify --out ...` command on SLURM; one CPU, 4 GiB and a suitable
`small` QoS are adequate starting resources. Verification rehashes native
scientific files and checks every retained native initial profile against its
exact compact basal. It then creates sigma-arm manifests and explicit
`submit_analysis.slurm`, `submit_audit.slurm` and `submit_collect.slurm` scripts
under each arm. No array concurrency throttle is added.

Run the analysis arrays before their dependent audit/export jobs. Collect
each arm's independent audit with the frozen audit module after all its audit
workers finish. Then run the frozen stage module's `collect --out ...` on
SLURM; allocate memory for the combined CSV population (start with 8–16 GiB
for large E2/E4 exports). Every final table is bound to completed worker,
independent-audit and source digests. No partial table becomes censoring.

## Observation populations

- E1 has separate sigma 0.25, 0.5, 1 and 2 analysis arms. Every model–size–input
  group must have exactly the same twelve latent seed blocks across these
  arms, and the underlying Gaussian-knot hashes must agree. Matching integer
  seeds alone is insufficient. Whole-seed bootstrap weights are shared across
  these paired arms. Sigma-zero mean controls remain deterministic controls.
- E2 augments the corrected E0 sigma-one grid with completed dense-input
  cases. Explicit finer-profile observation companions replace the same
  realization only after all shared native scalar arrays agree exactly. They
  are never counted as additional independent seeds.
- E3 augments E0 with qualified intermediate-size cases and their exact
  controls. The actual weakest-input baseline remains 0.0001. Missing weak
  baselines at new intermediate sizes make noise tolerance unavailable; the
  baseline is not silently replaced by 0.001.
- E4 is a separate long-horizon arm. Short E0 observations are not pooled
  into its persistence curves. Its means and no-input controls must match the
  long horizon. Comparisons retain the actual available matched strengths;
  unobserved conditions do not become lower bounds.

## S7/S8 and retained table contracts

The seven E0 tidy tables keep their measurement names. Added metadata include
`experiment`, `sigma`, `analysis_sigma`, `observation_arm`,
`native_horizon`, sampling cadence, threshold fields, `n_seeds` where a single
condition is identified, and the replication unit. The inherited `forced`
column remains the matched-mean union flag.

S7 additionally receives:

- `control_metrics.csv`: deterministic mean/zero persistence status, point or
  bounds, coverage and regime fractions at sigma zero. It has no invented
  bootstrap interval; `bootstrap_resamples=0`.
- `sigma_comparisons.csv`: paired whole-seed bootstrap intervals for
  persistence and short-lag mobility ratios relative to sigma one, and paired
  differences in regime fractions. Censored persistence ratios stay bounded
  or unidentified.
- `E1_forced_flags.csv`: the reviewed additional flag, low-sigma abnormal
  persistence, observed mean agreement and rule status. Condition metrics
  separately expose `forced_mean`, `forced_E1` and `forced_combined`.

S8 uses the same `condition_metrics`, `heading_curves`, `site_curves`,
`size_ratios` and `noise_tolerance` tables on the augmented sigma-one grid.
`source_experiments` distinguishes original from new conditions. A5 threshold,
optical and tracking statistics are separate root-owned outputs; these raw
heading/site tables do not change their definitions to mimic A5.

## Additional E1 forced-state rule

The authority is
`docs/tasks/briefs/20261020_figure4_e1_forced_rule.json`, reviewed by root before
E1 analysis. The sigma-one majority must be multiple-front or low-contrast;
the same abnormal majority must occur at both lower sigmas and in the exact
matched mean control. An unchanged healthy single front never triggers it.
Persistence at low sigma without mean agreement is reported separately and
is not proof of a forced state. A missing required arm/control leaves the
additional flag unavailable. The original matched-mean magnitude/contrast
union flag stays independent.

The tests include exact blocked Spring basal extraction, missing/mismatched
controls, false independent duplicates, changed latent fields despite matching
seed integers, unavailable weak-input baselines, explicit interim gating and
a synthetic worker → export → independent-audit roundtrip. No native E1–E4
scientific result was analyzed while implementing or testing this adapter.
