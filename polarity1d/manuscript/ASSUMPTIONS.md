# Assumptions for the corrected five-model comparison

Updated 2026-10-03. [PROTOCOL.md](PROTOCOL.md) defines the authorized matched
rerun: nominal amplitude factor 1, three current widths, existing lengths,
qualifying corrected establishment checkpoints, mirrored finite right pulses.
Wave-pinning, Goryachev, deBelly, Holmes 3 and Spring are active. Both LEGI and
Otsuji are excluded. Old campaigns supply only the personal field explorer.

## 1. Selected regimes, not complete model families

The adapters in [code/core/models](code/core/models/) retain specific native
regimes rather than exhausting each model family's parameter dependence.

| Model | Declared convention | Interpretation |
|---|---|---|
| Wave-pinning | Mori/Jilkine conserved adapter; initial active 0.2683, inactive 2 | Concentrations stay fixed with L, so total available material increases |
| Goryachev | Reduced Jilkine adapter; ec = 1 and u0 = 1 fill unspecified cytoplasmic Cdc24-Bem1 / initial active abundance | Explicit reproduction assumptions, not unique published defaults; inactive abundance follows the fixed point |
| Holmes 3 | No PI feedback; original Eqs. 5/6 factor 2; selected Rac activation 3.5; lowest stable homogeneous Rac root at rest | The publication supplies no unique Model 3 activation default; 3.5 is a selected tested regime |
| deBelly | Finite-volume port of staggered finite-element model; reflection-symmetric two-end extension; algebraic fixed-point start; cell-centred input; Rac/Rho basal activation 1; step 0.5 | These mechanics extend the original one-sided assay; state extension and discretization |
| Spring | Zmurchok author-model-derived moving domain; active start 0.05645 + 0.001*cos(pi*x); native L0 and time | Initial asymmetry is part of native rest; spontaneous polarization is allowed; L0 = 1 is not calibrated microns |

Current and frozen retained deBelly adapters explicitly choose cell-centred cue
sampling. The alternate `julia_element_right_edge` branch is not selected.
Intrinsic basal reactions remain active when extra stimulation is off; Spring's
basal activation 4 is distinct from extra cue amplitude 4.

## 2. Inputs, length and preparation

Holding concentration fixed as length changes is not fixed-total-mass scaling.
Spring includes dilution from domain motion. deBelly cues use nominal solver
coordinates with displacement monitored separately; Spring relative cues follow
physical length. Common width fractions retain these coordinate differences.
The current fixed supports are 1.5/3 or Spring 0.1/0.2; oversized support is
truncated by the domain without rescaling its native peak or prescribed width.
Equal native amplitude or spatial area does not imply equivalent receptor
occupancy, biochemical activation or force. Equal-integral controls are absent.

Only qualified nominal-factor-1 establishment states can start reversal.
Preserve full checkpoints and hashes; copy-based reproducibility checks must not
silently replace starting states. Spring's initial perturbation can seed no-cue
polarization: zero selected establishment duration means no cue needed, not
instantaneous physical establishment. Zero right duration is a no-reversal
control; passing it flags unstable preparation rather than a valid zero minimum.
LEGI's corrected positive basal withdrawal removes persistent polarity; Otsuji
is excluded by the user's decision. Their old results do not enter new fits.

## 3. Scalar readouts and field stability

Orientation is normalized active/Rac first moment, not complete spatial shape.
Contrast guards against large normalized scores of vanishing signal.
Scalar strength replacement compares the reversed first-moment magnitude with
initial strength within 1%; it is not reflected full-profile identity.
The final 500-unit scalar gate and independent 0.1% field-tail flag answer
separate questions. Preserve both and state fit masks explicitly. Neither scalar
success nor a stable tail alone proves long-time basin membership or the full
preparation-style continuation/holdout convergence of the response.

Frozen field scales use the per-field maximum of prior conditioning/release
transients, fixed across current width, amplitude and duration. Percent errors
are on those scales, not every tiny instantaneous value. deBelly displacement is
centered to remove rigid translation; physical length remains monitored.
The 1% scalar strength/arrival, duration-bracket and complete-profile clustering
rules are distinct. Apparent endpoints are sampled profile clusters, not attractors.
Zero sampled arrival means already within tolerance, not exact equilibrium.

## 4. Bounded search and descriptive length trends

Both scalar gates use exactly the same trials and all original final-500 samples.
Refine local fail-to-pass transitions; retain nonmonotonic intervals and unresolved
lower/upper bounds. Report shortest passing **tested** duration. A finite dense
plus log grid can miss narrow disconnected intervals; no universal monotone
threshold is claimed. Censored within tested range, preparation-ineligible and
missing execution are different statuses. Recovery observation is excluded from
the pulse-duration readout.

New length summaries preserve gaps, preparation exclusions, censoring and
field-tail flags. A fit conditional on successful eligible/stable-field cases
must show that membership, length range, sample count, raw slope and R-squared.
Native times and input amplitudes differ across models. Neither numerical
completion nor a success-only fit establishes a universal length law.
The retained length grid uses 201 spatial samples per model/length case. Physical
spacing therefore changes with nominal length; Spring's physical spacing also
changes as its domain moves. Chemical cues are evaluated on each solver's native
grid. This is not a fixed-physical-spacing comparison; a mesh sensitivity study
would be needed before attributing every small length effect to dynamics alone.

Retained checks address recorded calculations; this cleanup adds no systematic
mesh/time-step study, perturbation ensemble, parameter robustness screen or
physical unit calibration.

## 5. Historical diagnostic isolation

The personal explorer retains old model-specific widths, linear-conditioned
states and half-domain restriction. Old LEGI total-zero withdrawal arrests
response after activator/inhibitor decay; it is kinetic freezing, not basal
memory. Old Otsuji starting-reference ambiguity remains diagnostic. Neither
provides a current manuscript response or fit. Keep only required diagnostic
inputs/frozen source provenance and follow cleanup manifests.

The explorer displays native concentrations and actual time from pulse onset,
separate cutoffs and interpolation only within recorded histories. Its x50 is
half spatial integral for nonnegative fields, not peak/half-height; signed or
vanishing fields have undefined x50. Spring feedback
`80 * (physical_length - L0)` and restoring term
`0.01 * (physical_length - L0)` are distinct source-backed global proxies,
not calibrated spatial force fields. They derive from saved length histories.
