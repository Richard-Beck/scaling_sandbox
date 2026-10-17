# Inactive-pool perturbation: revised mechanism figures

The [Goryachev mechanism figure](../../../figures/Goryachev_length_mechanism/inactive_pool_perturbation/Goryachev_length_mechanism.png) and [Wave-pinning comparison](../../../figures/Goryachev_wave_pinning_trial/inactive_pool_perturbation/goryachev_wave_pinning_trial.png) use an alternative perturbation to the earlier incumbent activation-suppression tests. The [Goryachev report](../../../figures/Goryachev_length_mechanism/inactive_pool_perturbation/inactive_pool_report.pdf) has four pages (main figure and recovery for each duration). The [comparison report](../../../figures/Goryachev_wave_pinning_trial/inactive_pool_perturbation/inactive_pool_comparison_report.pdf) has both durations.

## Protocol

From each exact qualifying native steady checkpoint, remove 0.1%, 0.5% or 1% of the original inactive field at every spatial cell. Hold that material in a nonreacting external reservoir while the two fields evolve normally for 0.1 native time units (primary) or 1 unit (sensitivity). Then add back exactly the removed original spatial profile to the evolved inactive field. Active does not jump at removal or return. Native reactions and both diffusion coefficients stay unchanged throughout. This is a bolus removal/return, without a continuing clamp. An exactly zero-duration remove/return leaves the state unchanged, so it cannot give a nonzero active damage or a restoring reaction per unit damage.

Both models cover the same 18 lengths, 9–39. Each model/length/duration has a sham and three nonzero amplitudes: 72 workers, 288 native histories. No reversal cue is added. The original activation-suppression evidence and figures are preserved separately.

## Measurements

Damage is the sham-subtracted active amount lost inside the fixed original half-height incumbent core when inactive material is returned. Fractional loss divides by the corresponding sham core amount. Net reaction is native activation minus deactivation, integrated over that same fixed core; the excess subtracts the sham and its normalized value divides by the active amount lost at return. Active diffusion across the core is recorded separately, along with total material recovery rate. Fractional edge-cell finite-volume integration is used throughout.

The immediate reaction reflects both altered active and inactive fields after material return. This perturbation tests restoration after a global substrate disturbance; it does not isolate an autocatalytic coefficient or a global stability eigenvalue. Absolute amounts removed differ across models and lengths despite equal fractional inactive reduction. Native time and length units have not been calibrated between models.

Full steady concentrations retain physical positions and original amplitudes. Shape panels use L=21 and 0.5% removal, at return and the recorded sample near 9.93 native time units afterward. Width is the first interpolated crossing of each profile's own minimum plus half its contrast. This relative threshold separates falling height from relative edge motion; negative width loss means broadening. Original 201-cell grids are retained; no new mesh-convergence result is claimed.

## Findings under the alternative protocol

Longer Goryachev cells again lose a smaller fraction of their fixed-core active material and have stronger immediate net reaction per unit damage. Wave-pinning's remaining active deficit becomes concentrated near its incumbent edge during recovery, while Goryachev's is largely peak-centered. The shape contrast is time-dependent: both initially lose height. Small late peak overshoots and width changes remain visible.

The recovery tails belong to this protocol. The earlier activation-suppression statement that longer Goryachev cells recover more slowly late is not carried over: sustained 90% recovery after the 0.1-unit pool perturbation is faster at L=39 than L=9.

For 0.5% inactive removal:

| Model | Hold time | L | Core active loss (%) | Immediate reaction / damage | Sustained 90% recovery time |
|---|---:|---:|---:|---:|---:|
| Goryachev | 0.1 | 9 | 0.00100242 | 0.973092 | 18.8717 |
| Goryachev | 0.1 | 21 | 0.000790069 | 5.39402 | 8.80211 |
| Goryachev | 0.1 | 39 | 0.000487608 | 15.2854 | 2.63999 |
| Goryachev | 1 | 9 | 0.00828518 | 0.407071 | 19.6446 |
| Goryachev | 1 | 21 | 0.0038771 | 2.11645 | 21.2867 |
| Goryachev | 1 | 39 | 0.00161227 | 5.87305 | 21.2867 |
| Wave-pinning | 0.1 | 9 | 0.0491507 | 0.617337 | 8.12309 |
| Wave-pinning | 0.1 | 21 | 0.0484319 | 0.748702 | 5.66012 |
| Wave-pinning | 0.1 | 39 | 0.0482058 | 0.790194 | 3.94394 |
| Wave-pinning | 1 | 9 | 0.391922 | 0.485613 | 8.80211 |
| Wave-pinning | 1 | 21 | 0.36032 | 0.676465 | 6.91815 |
| Wave-pinning | 1 | 39 | 0.349967 | 0.748554 | 4.44866 |

## Evidence and reproduction

Original histories, solver metadata, frozen scientific sources, source checkpoint hashes, independent observable/return/material-reservoir audits and summary measurements live in `../inactive_pool_restoring_response_20261008/`. The native campaign ran on SLURM array 20760480 and finalizer 20760481; all 288 trials completed and passed independent audits.

`evidence_D0p1.npz` and `evidence_D1.npz` freeze compact steady fields, spatial reaction/deficit fields and matched-sham pulse/recovery observables. `response_vs_length.csv` includes both models, both durations and all three amplitudes. `morphology.csv` gives representative shape changes. `provenance.json` hashes all compact inputs and their sources. `render_provenance.json` hashes the renderer and figures.

The smallest-amplitude primary pulse was repeated at L=9 and L=39 in both models with rtol=2e-12, atol=2e-14 and max_step=0.01. Maximum relative difference in fractional loss or immediate reaction per damage was 0.000174409%. See the source campaign’s `results/accuracy_check.json`.

From the repository root, regenerate all alternative figures with:

```
/home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/inactive_pool_mechanism_figures.py render
```

For a fresh simulation campaign, use `inactive_pool_restoring_response.py init --durations 0.1 1 --out <fresh-directory>`. Original evidence is immutable. All producers are finished; source and figure generation times appear in the HTML reports.
