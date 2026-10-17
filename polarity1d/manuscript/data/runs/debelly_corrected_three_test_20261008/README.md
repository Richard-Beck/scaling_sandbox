# Corrected independent-endpoint deBelly replacement: completed

Generated 2026-10-08; producer finished. The flawed centered extension is deprecated. No old deBelly checkpoint or trial contributes to this replacement.

## Protocol and coverage

At all 22 original deBelly lengths (1, 3, 5, 7, 9, 11, 12, 13, 15, 17, 18, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39), independently start at homogeneous algebraic rest, apply nominal left relative-20% cosine patches (peak 3.6), and return to verified basal release. The model retains local Rac/tension protrusion and endpoint relaxation with `center_boundary_motion=False`. Endpoint commands/displacements are independent; inactive ends remain prescribed constraints. This is not a stress-free free-migration model.

The preparation search contains 688 native trials. Twelve preparations qualify (18, 19, 21–39 on the original grid), and ten shorter lengths do not. Twelve individual preparation trials hit the convergence cap; their flags are retained, and none supplies a selected checkpoint. All selected states have the original convergence windows and 10,000-unit holdout. The inherited eight-step preparation search left some lower brackets unresolved; a separately frozen bounded 56-step continuation resolves them, preserving `initial_eight_step_result.json`. These very small preparation pulse durations are numerical branch-selection outcomes under the finite convergence rules, not calibrated biological thresholds.

Each of the three reversal challenges restarts its length's exact selected full checkpoint: fixed support 3; relative support 20%; constant slope 0.24. Patch peak stays 3.6. Both scalar gates, 10,000-unit recovery, final-500 observation, separate native-field flag, dense/log duration grid and local nonmonotonic refinement retain their definitions. The duration grid and four local refinement rounds (final round index 5) resolve every sampled fail-to-pass transition within 1%.

The other four models retain their exact original inputs and immutable native trial archives. All 528 quantitative readout rows are unchanged; only copied archive paths differ. The full five-model suite has 330 cases / 660 readouts: 536 scalar successes, 46 censored and 78 preparation-ineligible readouts, with zero missing work or execution errors. Of 74,183 audited trials, 71,297 are exactly unchanged four-model evidence and 2,886 are newly integrated deBelly trials. All 72 eligible deBelly readouts succeed.

## What the corrected reversal curves show

Orientation-reversal minimum tested pulse durations (native seconds):

| Length | Fixed support 3 | Relative support 20% | Constant slope 0.24 |
| --- | ---: | ---: | ---: |
| 18 | 565.612 | 472.168 | 487.321 |
| 19 | 306.531 | 254.974 | 203.010 |
| 21 | 280.230 | 218.355 | 148.598 |
| 23 | 303.265 | 228.457 | 146.282 |
| 25 | 333.265 | 247.781 | 159.005 |
| 27 | 348.418 | 267.602 | 179.209 |
| 29 | 340.842 | 282.755 | 200.000 |
| 31 | 308.163 | 282.755 | 217.092 |
| 33 | 309.796 | 275.179 | 229.719 |
| 35 | 361.046 | 290.332 | 251.487 |
| 37 | 428.304 | 355.995 | 280.230 |
| 39 | 530.255 | 809.082 | 292.857 |

Above the slow-onset regime, the constant-gradient minimum duration increases from about 146 s at L=23 to 293 s at L=39. Patch curves contain nonmonotonic structure, so there is no universal monotonic length claim across all three geometries. These minimum successful pulses are distinct from continuous-cue MCA-zero times (289 and 686 s at L=23 and 39).

Native-field flags remain separate from scalar success. deBelly has only 3/3, 2/3 and 1/1 field-stable endpoints for orientation/strength in fixed/relative/gradient tests respectively. All six deBelly slope cells therefore show **No fit** under the existing minimum-four-points rule. The 24 fits for the other models pass independent closed-form OLS checks. Open markers preserve the scalar successes whose full fields have not passed the finite-window flag.

## Figures and evidence

- [Five-model × three-test duration figure](../../../figures/supplement_width_response.png), with PDF/SVG beside it.
- [Updated main figure and heatmap](../../../figures/main_polarity.png).
- [Corrected MCA-zero / Rac-relaxation / spontaneous-polarizability figure](../../../figures/deBelly_mechanism/debelly_zero_crossing.png).
- [Corrected advection-mechanism figure](../../../figures/deBelly_mechanism/debelly_advection.png).

`results/independent_audit.json` independently checks every imported and fresh native trial, inputs, no-cue controls, gates, field flags, coverage and endpoint selection. `preparation/refinement_independent_audit.json` verifies the extended preparation source/provenance and selections. `results/replacement_comparison_audit.json` verifies unchanged four-model quantitative results; `comparison_to_deprecated_centered.csv` records changed deBelly readouts. `results/fit_validation.json` independently verifies all 24 available fits.

Mechanism arrays and independent audits are in the separate sibling runs `debelly_corrected_zero_crossing_20261008`, `debelly_corrected_advection_20261008` and `debelly_corrected_advection_revision_20261008`. Both figures use the same fresh selected preparation checkpoints as this suite. They retain the original panel structure and the corrected prescribed-endpoint/direct-Rho/membrane-drag instantaneous decomposition. MCA centroid speeds and spatial profiles use nominal model coordinates; they are not displacement-mapped laboratory measurements.

Earlier canonical comparison files are archived under `figures/archive/centered_debelly_comparison_20261008/`. The old centered mechanism figure directories contain explicit deprecation notices. Original campaigns and native scientific solver files remain preserved.

## Reproduction

Use `/home/4473331/.conda/envs/scaling_sandbox/bin/python`, with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`. Fresh SLURM campaigns use `code/debelly_corrected_launch.py launch --out /absolute/path/to/FRESH_RUN`.

Completed-run regeneration:

```
python code/debelly_corrected_reversal_audit.py --campaign data/runs/debelly_corrected_three_test_20261008
python code/audit_debelly_preparation_refinement.py
python code/three_test_figures.py --campaign data/runs/debelly_corrected_three_test_20261008
python code/mechanisms/audit_debelly_corrected_advection.py
python code/mechanisms/audit_debelly_corrected_zero_crossing.py
python code/mechanisms/debelly_corrected_advection_figure.py
python code/mechanisms/debelly_corrected_zero_crossing_figure.py
```

Major completed jobs: preparation 20757828; extended brackets 20757915; preparation/reversal controller 20757916; grid 20758032; final controller 20758430; seeds 20757851; withdrawal 20757917; advection 20758022; mechanism audit/render 20758321. An initial unlaunched/canceled submission attempt and a stale completed-job dependency rejection are retained as execution history; corrected execution has no missing cases or errors.
