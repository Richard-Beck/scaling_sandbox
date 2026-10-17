# Main-text candidate 1, revision 2: incumbent substrate capture

![Six-panel revised candidate: native basal profiles and activation laws, actual substrate-removal damage and recovery, and returned-state allocation of net recovery to changed active field, inactive surplus and diffusion.](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/main_candidate_1_20261017/main_candidate_1.png)

[Vector PDF](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/main_candidate_1_20261017/supporting/vectors/main_candidate_1.pdf) · [SVG](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/main_candidate_1_20261017/supporting/vectors/main_candidate_1.svg) · [Technical caption and mathematical definitions](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/main_candidate_1_20261017/supporting/metadata/caption.md) · [Measurements](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/main_candidate_1_20261017/supporting/data/length_measurements.csv)

## Proposed message

Enlargement accompanies a taller Goryachev incumbent and a broader Wave-pinning incumbent, with the Goryachev state leaving less inactive substrate available for a remote cue. After a small inactive-pool disturbance and exact return, Goryachev's strongly length-dependent rebound is supplied by capture of the endogenous inactive surplus, while its active-field contribution at unchanged sham inactive concentration is slightly negative.

The candidate concerns **native substrate capture and remote substrate availability** in these chemical models. The returned-state allocation is conditional on this specified disturbance.

## Why this design serves the message

The model columns keep each starting state beside its own functional readouts, and length colors retain one meaning throughout. The magnified Goryachev incumbent makes concentration growth visible; Wave-pinning shows broadening. Goryachev's native conversion factor grows strongly with its taller incumbent. Wave-pinning's sampled plateau concentration changes little, so its local factor changes little over this range; its broader active extent still matters for integrated response. The exact local cue identity links lower remote inactive concentration to reduced initial cue-driven activation from that same pool.

The middle row shows the actual disturbance and response: 0.5% of the original whole inactive profile is removed at −0.1 native time units and returned pointwise at time zero, while native reactions and diffusion continue. The measured core is each length's original boundary-connected half-contrast region, held fixed in both perturbed and matched-sham histories. The first unit of recovery shows real active-material loss and subsequent conversion back to active material.

The final row identifies the chemical mediator. It evaluates the active-field change first at the sham inactive field, then the inactive-field change at the perturbed active field; diffusion completes the actual net rate. Under this exact ordering, the positive Goryachev chemical rebound comes entirely from the inactive-surplus contribution. Wave-pinning has positive contributions from both fields. Linear axes preserve their relative magnitudes; Goryachev's small negative active-field contribution is stated directly.

“Capture” means conversion sensitivity to extra inactive material at fixed active concentration, $A(u)=\partial f/\partial v$. From L9 to L39, Goryachev's surplus-weighted coefficient grows **16.81-fold**, while core surplus per core damage changes **0.913-fold**, giving **15.35-fold** substrate allocation; Wave-pinning's corresponding product is **1.093 × 1.060 = 1.159-fold**. The [R3 factorization data](reviews/incumbent_r3_capture_factorization.csv) and [exact definition](reviews/incumbent_r3.md) identify both factors at the actual coupled states.

The surplus is endogenous: **after exact return, global inactive surplus balances global active loss** relative to sham. The experiment returns exactly the material withheld; the subsequent rebound reflects native conversion of the redistributed pool.

## What the retained evidence shows

Each entry below compares **L9 → L39**, in model-native units. The gain panels include all 18 retained lengths within each model.

| Observable or returned-state allocation | Goryachev | Wave-pinning |
| --- | --- | --- |
| Original half-contrast incumbent width | 0.712 → 0.799 | 2.821 → 11.919 |
| Peak active concentration | 27.76 → 117.39 | 1.101 → 1.125 |
| Opposite-end inactive concentration, hence initial added reaction per unit local cue | 0.2633 → 0.08771; **0.333×** | 1.80782 → 1.80764; **1.000×** |
| Fractional core active damage at return | 0.001002% → 0.000488% | 0.04915% → 0.04821% |
| Active-field change at sham inactive concentration, per damage | **−0.02438 → −0.02346** | 0.09699 → 0.18724 |
| Inactive-surplus capture at perturbed active concentration, per damage | **0.99747 → 15.30885** | 0.52035 → 0.60295 |
| Native diffusion difference, per damage | 0.00456 → 0.00410 | 0.01728 → 0.00367 |
| Immediate net rate per damage, including diffusion | 0.978 → 15.289; **15.64×** | 0.635 → 0.794; **1.25×** |

All allocated rates are integrated over the original core and divided by its active loss at return, with units of inverse native time. Their sum is the net rate. The ordering assigns the nonlinear cross term to the inactive contribution; it is an algebraic allocation at coupled states, not separately evolved ablations. Goryachev's active-field term remains negative at all 18 lengths, including after adding its actual diffusion difference.

Mean conserved-pool concentration and local kinetic parameters are fixed within each model, so total available material grows with length by construction. Absolute native rates are not calibrated between models.

## Sources, checks and remaining gaps

The source is the current [length-mechanism PDF](../polarity1d/manuscript/figures/length_mechanisms_20261008/length_mechanisms.pdf), Figures 3–4, with original checkpoints and the [retained inactive-profile perturbation evidence](../polarity1d/manuscript/data/runs/inactive_pool_figure_exports_20261008/README.md). This bounded reanalysis uses no new integrations. It reconstructs all 36 exact-return budgets, the ordered split and global conservation, and checks both phases of all six displayed matched history pairs against the original observables. [Verification](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/main_candidate_1_20261017/supporting/metadata/verification.json) and [source/artifact hashes](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/main_candidate_1_20261017/supporting/metadata/manifest.json) accompany the exports.

Goryachev's displayed profile is an explicitly labelled **x = 0–3 incumbent-end view**. The opposite-end measurements come from its full saved state. [Full basal profiles](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/main_candidate_1_20261017/supporting/data/basal_profiles.csv) and [complete retained damage/recovery series](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/main_candidate_1_20261017/supporting/data/core_damage_and_recovery.csv) are exported, including all 5000 native recovery units.

The native law and measured allocation give a physical substrate-capture interpretation, but the experiment does not independently measure a general capture capacity or isolate the architectural change from the rest of the state. The same fractional intervention produces different field perturbations and active damage. It creates no rival, so successful competition and takeover remain unresolved, as do global stability and a connection to WGD-associated migration persistence.

## Response to fresh reviews

[R2](reviews/incumbent_r2.md) changed the interpretation. [Version 1](main_candidate_1_v1_20261017.md) showed that reaction-only and net rates nearly coincide, but left the source of rebound unidentified. The exact returned-state split now exposes substrate-surplus capture as the positive Goryachev chemical contribution; the figure replaces the redundant overplot with this allocation and states the native activation laws. Its message is shared-pool substrate capture and remote availability, rather than generic intrinsic incumbent maintenance. Version 1's figure, data, captions, renderer and original Markdown are preserved in its linked archive.

[R3](reviews/incumbent_r3.md) accepts this six-panel candidate and asks for one numerical clarification. Its supplied coefficient-times-surplus factorization is now stated above and in the caption, making the growing native property explicit without another panel. Wave-pinning's description now refers to its nearly unchanged sampled plateau rather than implying deep saturation.

These remain two chemical examples for discussion, with no final story selection. The separate distance hypothesis is not tested here.
