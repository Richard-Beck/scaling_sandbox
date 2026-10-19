# Common-initial-arc stationary-cue arm

Generated: 2026-10-19T18:23:35.591976+00:00. Complete: 240/240 requested conditions. No overlap is excluded.

Chemical rings use C = 20/30/40 and initial arc offsets 1.25/2.5/5/10. Spring uses resting C₀ = 1/2/3 and the same four physical initial offsets defined by the exact qualified C₀=1 basal circumference times 1/16, 1/8, 1/4 and 1/2. Cue centres remain stationary in material coordinates; their physical offset subsequently scales with Spring's instantaneous circumference. Native supports, peaks, durations, basals and 13,000-unit cue-free continuation retain the first arm's definitions. Distance and angle remain linked at each fixed circumference.

- [Current 240 requested conditions and available readouts](supporting/data/stationary_arclength_summary.csv)
- [Initial context, overlap and exact cue shapes](supporting/data/stationary_arclength_initial_context.csv)
- [Initial native profiles and cue samples](supporting/data/stationary_arclength_initial_profiles_and_cues.csv)
- [Saved regional first-entry and retained-entry brackets](supporting/data/stationary_arclength_crossings.csv)
- [Final native primary profiles](supporting/data/stationary_arclength_final_profiles.csv)
- [Both arms: 480 rows](supporting/data/stationary_two_arms_all480.csv) and [355 unique protocols](supporting/data/stationary_two_arms_unique355.csv)
- [Paired design levels](supporting/data/stationary_two_arms_reference_level_pairs.csv): levels coincide at chemical C20 / Spring C0=1, then compare fixed angle with fixed initial arc across size. Shared protocol IDs separately capture all 125 exact reused conditions, including matches at different design-level indices.
- [Latest collection status](supporting/metadata/latest_status.json), [validation](supporting/metadata/validation.json) and [source hashes](supporting/metadata/source_sha256.json)

Raw M₁ components use fixed physical ring axes. P = M₁/M₁(0) aligns and scales the initial vector once; R = |P|, physical M₀ is total active amount, and q₁ = |M₁|/M₀. Moment quadrature is instantaneous C times the native-cell mean. Endpoint angles are masked at R < 0.1; raw phases remain diagnostics. Compare alignment using wrapped target_error_final_deg and the actual target_angle_in_P_coordinates_deg, which includes the initial moment's measured offset from the preparation axis. Fraction of requested turn uses the raw signed endpoint angle only when R ≥ 0.1; its branch discontinuity near ±180° makes it a raw diagnostic rather than a universal progress or success measure. The fixed original material region subtracts its transported initial background amount and does not track a site's ancestry. Event times are saved-sample brackets, separate for cue-on and release; sparse duration levels do not provide refined minimum cue durations. Full-field tail flags are empirical finite-window excursions; retained primary-only archives keep the full-field flag unavailable.

All compact files replace atomically; latest_status.json is written last and hashes that generation's files. Missing native branches remain pending/error rows. These exports contain no new integration or figure layout. The old angle-arm files are unchanged.
