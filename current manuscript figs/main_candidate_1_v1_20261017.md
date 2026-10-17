# Main-text candidate 1: incumbent organization and local chemical response

![Six-panel candidate in Goryachev and Wave-pinning columns: basal profiles, actual substrate-removal damage and recovery, and immediate reaction-only versus net recovery per damage.](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/archive_main_candidate_1_v1_20261017/main_candidate_1.png)

[Vector PDF](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/archive_main_candidate_1_v1_20261017/supporting/vectors/main_candidate_1.pdf) · [SVG](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/archive_main_candidate_1_v1_20261017/supporting/vectors/main_candidate_1.svg) · [Technical caption and mathematical definitions](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/archive_main_candidate_1_v1_20261017/supporting/metadata/caption.md) · [Measurements](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/archive_main_candidate_1_v1_20261017/supporting/data/length_measurements.csv)

## Proposed message

Enlargement accompanies a taller, substrate-depleting incumbent in Goryachev and a broader incumbent in Wave-pinning. In the specified small inactive-pool removal/return experiment, longer Goryachev cells sustain less fractional incumbent damage and have a substantially greater immediate net core-repair rate per damage; their depleted inactive pool also reduces the initial source driven by an equal local opposing cue.

This is a provisional figure about **architecture and explicitly local chemical responses**. It does not establish incumbent scale as an independently manipulated cause of global stability.

## Why this design serves the message

The model columns keep each starting state beside its own functional readouts, and length colors retain one meaning throughout. The magnified Goryachev incumbent makes concentration growth visible; the Wave-pinning profiles expose broadening. The exact local source identity gives inactive depletion a functional interpretation without adding a curve that merely repeats inactive concentration.

The middle row shows the actual disturbance and response: 0.5% of the original whole inactive profile is removed at −0.1 native time units and returned pointwise at time zero, while native reactions and diffusion continue. The measured active core is each length's original boundary-connected half-contrast region, held fixed during both perturbed and matched-sham histories. The plotted first unit of recovery makes the perturbation concrete; fractional damage is shown directly rather than hidden by normalization to the damage itself.

The final row asks whether the earlier reaction-only diagnostic represents actual immediate repair of that core deficit. Independently adding native diffusion gives the net rate; both are divided by damage at return. Their agreement shows that the strong Goryachev trend persists in the **net original-core response** under this intervention. Separate axes retain native rate units; absolute rates are not calibrated between models.

## What the retained evidence shows

Each entry below compares **L9 → L39**, in model-native units. The gain panels include all 18 retained lengths within each model.

| Observable | Goryachev | Wave-pinning |
| --- | --- | --- |
| Original half-contrast incumbent width | 0.712 → 0.799 | 2.821 → 11.919 |
| Peak active concentration | 27.76 → 117.39 | 1.101 → 1.125 |
| Opposite-end inactive concentration, hence initial added reaction per unit local cue | 0.2633 → 0.08771; **0.333×** | 1.80782 → 1.80764; **1.000×** |
| Fractional core active damage at return | 0.001002% → 0.000488% | 0.04915% → 0.04821% |
| Immediate reaction-only rate per damage | 0.973 → 15.285; **15.71×** | 0.617 → 0.790; **1.28×** |
| Immediate net rate per damage, including diffusion | 0.978 → 15.289; **15.64×** | 0.635 → 0.794; **1.25×** |

Mean conserved-pool concentration and local kinetic parameters are fixed within each model. Total available material therefore grows with length by construction. The informative comparison concerns how that material is organized and how the resulting states respond; increased amount alone does not predict the displayed response trends.

## Sources, checks and remaining gaps

The source is the current [length-mechanism PDF](../polarity1d/manuscript/figures/length_mechanisms_20261008/length_mechanisms.pdf), Figures 3–4, with original basal checkpoints and the [retained inactive-profile perturbation evidence](../polarity1d/manuscript/data/runs/inactive_pool_figure_exports_20261008/README.md). This candidate uses no new integrations. The reanalysis reconstructs the exact-return reaction, native finite-volume diffusion and net rate at all 36 model–length states; it also checks both phases of all six displayed matched history pairs against their original core observables. [Verification](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/archive_main_candidate_1_v1_20261017/supporting/metadata/verification.json) and [source/artifact hashes](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/archive_main_candidate_1_v1_20261017/supporting/metadata/manifest.json) accompany the exports.

Goryachev's displayed profile is an explicitly labelled **x = 0–3 incumbent-end view**. The opposite-end measurements come from its full saved state. [Full basal profiles](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/archive_main_candidate_1_v1_20261017/supporting/data/basal_profiles.csv) and [complete retained damage/recovery series](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/archive_main_candidate_1_v1_20261017/supporting/data/core_damage_and_recovery.csv) are exported, including all 5000 native recovery units.

Architecture, inactive depletion and response covary; the figure does not isolate their causal contributions. The same fractional inactive intervention causes different active damage across lengths and models. Immediate net repair is one response direction within one original core, rather than a global stability measure, a force or a general late recovery rate. The perturbation creates no rival, and the initial opposing-cue source does not establish challenger formation, survival or takeover. These polarity-model diagnostics do not measure spontaneous migration persistence or calibrate a biological WGD enlargement.

The editorial decision remains whether these two chemical examples carry the desired main-text insight, or whether the retained rapid challenger activation/persistent incumbent contrast is the more important story. No final story selection is made here; the separate distance hypothesis is not tested by this figure.
