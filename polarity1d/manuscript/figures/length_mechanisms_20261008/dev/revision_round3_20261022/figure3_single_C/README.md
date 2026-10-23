# Figure 3 round 3

[Main PDF](figure3_round3.pdf) · [legend](legend.md) · [complete supplement](figure3_round3_supplement.pdf) · [reader](index.html).

Generated 2026-10-22T20:03:17.235122+00:00; four relaxation protocols complete. Original evidence is immutable. No kinetic retuning or new response campaign. Readout source: `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261020/figure3_distance/supporting/data/A1_heading_target_time_brackets.csv`, target_fraction = 0.5, coherence_floor = 0.05, non-control placements. Coincident cue/heading axes have no nonzero turn and are not log-axis observations. All other placements remain in the coordinate CSV, including undefined headings.

Time coordinates are roots of the linear native-field interpolant within the original saved brackets, rather than new integrations. Complete bracket endpoints and reciprocal speed intervals are exported. Topology is classified by the existing two-region watershed/saddle detector at that same interpolated state; endpoint agreement is checked. Undefined-heading cells extend halfway to adjacent tested gaps, separately for each size. Grey bands show their union across sizes; only the affected size's curve is interrupted. Regime lines report the first tested saddle at 50%, not an inferred continuous bifurcation. No-saddle sizes have no line. The 50% source contains no positive-target censored observations at the selected gate; triangle meanings remain defined for reproducibility. There are 1 crossings whose bracket endpoints disagree on two-region topology.

A's cue-on profiles reuse the original full fields, reconstructed at 50/90% inside saved brackets. Only its relaxation row is new. Four cue-on replays start from the exact immutable basal input, using the original frozen model and BDF settings (rtol 10⁻⁹, atol 10⁻¹¹, max step 10). They stop at the dense 95% heading event or 20,000. Cue-off starts from exactly that full state. Existing active-field scaled-RHS tolerance is 10⁻⁴ over 1,000 time units and active scale 2.2683, hence max |∂u/∂t| ≤ 2.2683 × 10⁻⁷. The original all-species basal qualification gate is distinct from this requested active-species relaxation stop. Full cue-on/cue-off states, accepted-step mass/positivity checks and stop conditions are in `polarity1d/manuscript/data/runs/figure3_round3_relaxation_v2_20261022`.

Every round-2 main and supplement page is appended to the new supplement, after 90/99%-heading pages. Original files are untouched. Native clocks have no cross-model calibration; observed migration and conditional OU persistence remain separate evidence from these simulated polarity readouts.

Rebuild using the existing scaling_sandbox environment:

```sh
python -B polarity1d/manuscript/code/mechanisms/revision_round3_figure3_20261022.py freeze
# Submit the exact sbatch_arguments in the generated submission.json once.
python -B polarity1d/manuscript/code/mechanisms/revision_round3_figure3_single_C_20261022.py
```

This layout revision places the eight C = 30 profiles on the left and B above C on the right. Only i–ii are highlighted in the main. The previous four-placement main and its legend are appended to the supplement, retaining the C = 47.6 profiles and their horizon qualification. All original source CSVs and four completed relaxation protocols are reused; no simulations were run for this layout revision.
