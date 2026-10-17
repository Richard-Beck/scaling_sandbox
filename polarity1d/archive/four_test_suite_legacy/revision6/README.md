# Steady-state four-test suite rerun (revision 6)

Authorized by the user on 2026-09-29. **Spring competition remains paused** by the user's explicit clarification. Periodic following retains the existing preparation; the new fixed-cue settling requirement does not apply to periodic following.

## Campaign

- 770 assay records across the complete revision-5 model/length manifest (201 cells, selected fixed regimes).
- 441 preparation jobs: competition 133 (spring excluded), gradient reversal 154, finite countercue 154.
- Each constant preparation phase has a maximum of 10,000 native time units. These are seconds for deBelly and Holmes, uncalibrated model time for the four original RD models and spring. The cap applies separately to each cue/no-cue phase, not to the response/recovery horizon.
- Constant conditioning cues preserve the accepted amplitudes and spatial shapes. Duration is replaced by detected settling time. Competition cues are held at full amplitude while settling; the original 20-unit withdrawal taper is retained for the four original RD models. Holmes/deBelly withdraw abruptly. Spring competition is never simulated.
- Both reversal tests are included separately. Gradient reversal uses the settled initial left gradient and its original 20-unit crossfade. Finite countercue settles under its own initial conditioning cue, then under no cue before challenging. Spring finite countercue starts from the settled spontaneous state and reflects its orientation left if needed, retaining physical length. Spring gradient preparation additionally settles under the sustained left gradient.
- Qualification requires all evolving chemistry and mechanical state to settle. Over a 500-unit window, maximum pointwise field range and final drift projected across that window must each be <=0.1% of the field's maximum absolute value (floor 1e-8). RD/spring use the full final RHS; deBelly uses the last saved finite difference. Physical length and displacement shape are normalized to reference L. Rigid translation is treated as a coordinate gauge. deBelly finite autonomous preparation uses a 1000-unit window to retain its existing Rac/MCA half-mass stability checks.
- Detection occurs after the confirmation window, then at 100-unit intervals. Reported times include the confirmation window and have 100-unit resolution. This is an operational numerical criterion, not a proof of asymptotic stability.
- Failure to settle disqualifies that model/test/length combination. A settled preparation can still fail the accepted polarity/peak eligibility conditions. These outcomes remain distinct in data; figures use grayscale for unavailable or unusable values without distinct failure symbols.
- Full finite-countercue duration grids are copied unchanged from revision 5. Every duration restarts independently from the same hashed complete prepared checkpoint. No bisection, monotonicity assumption, or partial-sweep minimum is used.
- Periodic following is rerun with the accepted revision-3 initialization, cue, sampling, harmonic fit and reliability gates, over the full revision-5 length list. Both P=100 and P=400 are included.

## Complete state restart

`debelly_checkpoint_solver.py` copies the canonical specialized solver, adds a full-state checkpoint input and records the raw boundary variables needed to restart. Its kinetics, coefficients and update order are unchanged. `check_checkpoint.py` verified exact equality with the unmodified default solver and exact agreement between uninterrupted and split runs, including after a cue change, at 201 cells. All seven models' restart interfaces were checked at the production mesh. `smoke/` retains the checks.

RD checkpoints include every species. Spring checkpoints include both species and physical length (translation gauge reset), not a restored resting length. deBelly checkpoints include Rac, Rho, bound/unbound MCA, displacement, velocity, tension, centered and raw boundary state. Original simulation sources and accepted prior results are preserved.

## Launch and outputs

`launch.json` records SLURM IDs and stages. `frozen_provenance.json` freezes all executable inputs, model sources, manifest and duration grids. Workers verify sources before/after computation. Preparation checkpoints and response manifests are hashed.

1. `preparation.slurm` runs the 441 preparations.
2. `controller.slurm`, after the preparation array finishes, plans qualified responses, disqualifies settling/eligibility failures, flags missing execution outputs separately and submits the response array (10 durations per finite chunk).
3. `report.slurm`, after response completion, collects full duration coverage and generates figures.
4. `watch.slurm` regenerates timestamped report files as outputs become available. The browser never refreshes itself. Missing outputs and incomplete sweeps remain pending until final collection; execution failures are listed separately.

Figures and CSVs are in [results](results/README.md): `four_test_all_models.{png,pdf,svg}`, `time_to_settle.{png,pdf,svg}`, `endpoints.csv`, `settling_times.csv`, and `status.json`. These are revision-6 results only. The HTML report is a timestamped view of currently generated files. Until the campaign finishes, figures are preliminary. The accepted assay explanation v2 remains unchanged and is a keeper. This campaign does not generate the postponed length-effect heatmap.

To reconstruct the campaign from a fresh revision-6 directory: run `campaign.py setup`, then `launch.py` using the existing `scaling_sandbox` Conda environment. Do not rerun setup or change frozen sources while this campaign is active. Do not invoke launch twice. Failed execution jobs require diagnosis; they are not scientific evidence of failure to settle.

### Figure orchestration correction

Legacy adapters modify Python's module search path, so a bare `import build_report` in the original frozen launcher resolved to revision 3. This was discovered after the response array was submitted. The report/progress jobs were replaced with jobs using [the explicit revision-6 report service](../revision6_report_service/service.py), without changing or canceling scientific simulation jobs or frozen simulation inputs. Job replacements and supplementary driver hashes are recorded in `report_orchestration_repair.json`. The earlier revision-3 report was inadvertently regenerated from its own unchanged data; revision-6 outputs now use the correct report implementation. For manual regeneration, use `revision6_report_service/service.py report`; do not invoke the frozen watch/finalize modules directly. Fix the import resolution when adapting this snapshot into any future campaign.
