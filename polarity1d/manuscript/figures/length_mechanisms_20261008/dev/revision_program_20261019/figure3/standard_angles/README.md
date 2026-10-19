# Standardized angle–duration checkpoint

This is the user-requested stationary-cue experiment at π, π/2, π/4 and π/8. It preserves the previous Figure 3 sheets while the experiment and interpretation are reviewed.

- [Orientation response](standard_angles_orientation.png) · [PDF](supporting/vectors/standard_angles_orientation.pdf)
- [Activity-moment magnitude](standard_angles_magnitude.png)
- [Initial basal profiles and cue footprints](standard_angles_initial_context.png)
- [Activity-vector concept sketch](vector_concept/stationary_vector_concept.png)
- [Actual signed spatial changes in the selected histories](vector_concept/stationary_vector_native_changes.png)
- [Compact endpoint table](supporting/data/standard_angles_endpoints.csv)
- [Complete endpoints, native provenance and field flags](supporting/data/standard_angles_summary.csv)
- [Input anchor and overlap diagnostics](supporting/data/standard_angles_initial_context.csv)
- [Caption](supporting/metadata/CAPTION.md) · [Evidence and definitions](supporting/metadata/EVIDENCE.md) · [Validation](supporting/metadata/validation.json)

All 80 conditions are complete: 48 new native branches and 32 proved exact archive reuses. All 35 initial-overlap conditions remain. Chemical tests use fixed C = 40 and native durations 1/10/100/300/1000; Spring uses resting C₀ = 2 and native durations 0.1/1/10/100/1000. Every endpoint follows 13,000 cue-free units. Cue-centre angles are shared; native physical widths, peaks and clocks remain model-specific.

All 60 non-antipodal final primary profiles have one resolved peak and match a rotated copy of their exact basal profile to within 2.21×10⁻¹² of its contrast. Their graded retained angles describe one translated full-strength activity region. Antipodal Holmes instead loses first-moment magnitude: at T = 1000, R ≈ 0.040 and direction is masked. Native profiles and magnitude distinguish these outcomes.

Raw angle and R are primary. Target error and fractional requested turn remain exported; no refined minimum-duration or arbitrary small-angle alignment claim is made. Overlapping inputs perturb the original activity region, so fixed-region loss does not identify a separate challenger or site lineage. Activity orientation is not measured cell movement.

The [frozen native manifest](../../../../../../data/runs/revision_standard_angles_duration_20261019/manifest.json) and [collection counts](../../../../../../data/runs/revision_standard_angles_duration_20261019/results/collection.json) record the full experiment. SLURM array: 21187059. New native allocated storage: 334.7 MiB.

To rebuild the three plotted sheets from retained compact tables, run `revision_standard_angles_report_20261019.py --render-only` in the existing environment. The default command independently replays the native archives before generating tables, plots and validation; that replay requires the preserved native data checkout.
