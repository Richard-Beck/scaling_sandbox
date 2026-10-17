# Full data catalog for the `migration1d` experiments

## Purpose

This catalog describes the data available beyond the descriptive tables in `EXPERIMENT_SUMMARIES.xlsx`. Paths are relative to the repository root. The catalog distinguishes final merged scientific outputs from task-level and scheduler-oriented intermediates. It does not rank results or recommend how they should be interpreted or visualized.

## Data tiers

| Tier | Contents | Use | Included in writer bundle? |
|---|---|---|---|
| 1 | Workbook summaries | Initial descriptive review | Yes |
| 2 | Final merged CSV and JSON files | Condition-level reanalysis and verification | Described here; retained in the working repository |
| 3 | Moving-interval task JSON and compressed NPZ trajectories | Time-resolved and spatial reanalysis | Described here; retained in the working repository |
| 4 | WGD shards and extended per-task JSON | Merge validation, recovery, and condition-level audit | Described here; normally unnecessary after merge |
| 5 | Generated HTML reports and image files | Existing presentation layer | Described here; not required to reproduce workbook summaries |

The active organized outputs under `migration1d/results/single_species` and `migration1d/results/three_species`, together with the top-level directional atlas CSV, are the sources used for this bundle. Legacy/pre-reorganization outputs elsewhere in the repository are not cataloged as active results.

## Final condition-level CSV files

These are the preferred inputs for analyses that can be performed from terminal outcomes or precomputed trajectory summaries.

| Path | Rows | Columns | Size | Experimental grain | Description |
|---|---:|---:|---:|---|---|
| `migration1d/results/focused_motile_length_directional_atlas.csv` | 960 | 15 | 171.86 KiB | one moving-interval trajectory | Derived cumulative leftward/rightward center travel for all moving-interval trajectories. |
| `migration1d/results/single_species/calibration/local_bistability/roots.csv` | 353 | 3 | 11.11 KiB | one reaction root at one inactive concentration | Roots of the local single-species reaction at each scanned inactive concentration. |
| `migration1d/results/single_species/calibration/regime_map/regime_map.csv` | 867 | 11 | 95.94 KiB | one 300 s trajectory/parameter condition | Single-cue single-species calibration over mu, Lambda, and h. |
| `migration1d/results/single_species/testing/dual_pole_map/dual_pole_map.csv` | 2,028 | 8 | 96.49 KiB | one 300 s trajectory/parameter condition | Two-seeded-domain test over mu, Lambda, h, and epsilon. |
| `migration1d/results/single_species/testing/focused_motile_length/focused_motile_length.csv` | 240 | 46 | 183.76 KiB | one 500 s trajectory/design condition | Single-species moving-interval sweep with global elastic-tension feedback. |
| `migration1d/results/single_species/testing/focused_motile_length_spatial_tension/focused_motile_length.csv` | 240 | 54 | 204.39 KiB | one 500 s trajectory/design condition | Single-species moving-interval sweep with a spatial tension signal. |
| `migration1d/results/single_species/testing/polarity_reversal/polarity_reversal_map.csv` | 867 | 11 | 95.79 KiB | one 300 s trajectory/parameter condition | Sequential left/right cue test over mu, Lambda, and h. |
| `migration1d/results/single_species/testing/wgd_null/double_stimulus/two_seeded_domains.csv` | 6,084 | 17 | 1.05 MiB | one baseline condition, epsilon, and volume multiplier | Single-species WGD-null two-seeded-domain extension. |
| `migration1d/results/single_species/testing/wgd_null/polarity_reversal/polarity_reversal.csv` | 2,601 | 20 | 627.33 KiB | one baseline condition at one volume multiplier | Single-species WGD-null polarity-reversal extension. |
| `migration1d/results/single_species/testing/wgd_null/wgd_null_screen.csv` | 2,601 | 17 | 531.90 KiB | one baseline condition at one volume multiplier | Single-species single-cue WGD-null endpoint screen. |
| `migration1d/results/three_species/calibration/tier1_regime_map/tier1_regime_map.csv` | 507 | 41 | 301.05 KiB | one 300 s trajectory/parameter condition | Holmes six-field tier-1 calibration over Lambda, kappa_R, and n; exported from the 1x WGD slice. |
| `migration1d/results/three_species/testing/focused_motile_length/focused_motile_length.csv` | 240 | 58 | 232.15 KiB | one 500 s trajectory/design condition | Three-species moving-interval sweep with global elastic-tension feedback. |
| `migration1d/results/three_species/testing/focused_motile_length_spatial_tension/focused_motile_length.csv` | 240 | 66 | 246.15 KiB | one 500 s trajectory/design condition | Three-species moving-interval sweep with a spatial tension signal. |
| `migration1d/results/three_species/testing/wgd_null/double_stimulus/two_seeded_domains.csv` | 6,084 | 21 | 1.54 MiB | one baseline condition, epsilon, and volume multiplier | Three-species WGD-null two-seeded-active-Rac extension. |
| `migration1d/results/three_species/testing/wgd_null/polarity_reversal/polarity_reversal.csv` | 2,601 | 22 | 786.97 KiB | one baseline condition at one volume multiplier | Three-species WGD-null active-Rac reversal extension. |
| `migration1d/results/three_species/testing/wgd_null/wgd_null_screen.csv` | 1,521 | 39 | 871.65 KiB | one baseline condition at one volume multiplier | Three-species single-Rac-cue WGD-null endpoint screen. |

## CSV schemas and units

Every final CSV schema is listed below in file order. Units are inferred from the implemented field definition; fields in model units are identified as such rather than assigned physical units that were not calibrated.

### `migration1d/results/focused_motile_length_directional_atlas.csv`

Rows: 960. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `condition_index` | categorical, Boolean, or model units; see METHODS.md |
| `condition` | categorical, Boolean, or model units; see METHODS.md |
| `chemistry_class` | categorical, Boolean, or model units; see METHODS.md |
| `tension_mode` | categorical, Boolean, or model units; see METHODS.md |
| `task_index` | categorical, Boolean, or model units; see METHODS.md |
| `reference_length_um` | um |
| `mechanical_case` | categorical, Boolean, or model units; see METHODS.md |
| `protocol` | categorical, Boolean, or model units; see METHODS.md |
| `tension_feedback_alpha` | dimensionless/model control |
| `leftward_travel_um` | um |
| `rightward_travel_um` | um |
| `total_travel_um` | um |
| `leftward_fraction` | dimensionless fraction |
| `rightward_fraction` | dimensionless fraction |
| `direction_class` | categorical, Boolean, or model units; see METHODS.md |

### `migration1d/results/single_species/calibration/local_bistability/roots.csv`

Rows: 353. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `inactive` | categorical, Boolean, or model units; see METHODS.md |
| `active` | categorical, Boolean, or model units; see METHODS.md |
| `stable` | categorical, Boolean, or model units; see METHODS.md |

### `migration1d/results/single_species/calibration/regime_map/regime_map.csv`

Rows: 867. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `mu` | dimensionless/model control |
| `Lambda` | dimensionless/model control |
| `h` | dimensionless/model control |
| `terminal_polarity_status` | categorical, Boolean, or model units; see METHODS.md |
| `time_to_polarization` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_polarization_amplitude` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_pole_count` | count |
| `terminal_high_width_fraction` | dimensionless fraction |
| `terminal_polarity_contrast` | dimensionless |
| `late_maximum_normalized_change` | categorical, Boolean, or model units; see METHODS.md |
| `maximum_conservation_error` | dimensionless relative error |

### `migration1d/results/single_species/testing/dual_pole_map/dual_pole_map.csv`

Rows: 2,028. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `mu` | dimensionless/model control |
| `Lambda` | dimensionless/model control |
| `h` | dimensionless/model control |
| `epsilon` | dimensionless/model control |
| `polarized_domain_count_t0` | count |
| `polarized_domain_count_t100` | count |
| `polarized_domain_count_t300` | count |
| `maximum_conservation_error` | dimensionless relative error |

### `migration1d/results/single_species/testing/focused_motile_length/focused_motile_length.csv`

Rows: 240. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `chemistry_class` | categorical, Boolean, or model units; see METHODS.md |
| `dashpot_viscosity` | categorical, Boolean, or model units; see METHODS.md |
| `duration_seconds` | s |
| `final_length_um` | um |
| `fraction_time_under_tension` | s unless stated as a sampled index |
| `integrated_tension_feedback_exposure` | categorical, Boolean, or model units; see METHODS.md |
| `maximum_absolute_strain` | dimensionless |
| `maximum_active_gtpase_conservation_error` | dimensionless relative error |
| `maximum_elastic_tension` | categorical, Boolean, or model units; see METHODS.md |
| `maximum_length_um` | um |
| `maximum_tension_added_deactivation_rate` | model rate units; generally s^-1 for chemical rates |
| `mean_absolute_speed_cue` | um/s |
| `mean_absolute_speed_early` | um/s |
| `mean_absolute_speed_established` | um/s |
| `mean_absolute_speed_late` | um/s |
| `mean_absolute_speed_reorientation` | um/s |
| `mean_absolute_speed_reversal_cue` | um/s |
| `mean_polarity_projected_speed_cue` | um/s |
| `mean_polarity_projected_speed_early` | um/s |
| `mean_polarity_projected_speed_established` | um/s |
| `mean_polarity_projected_speed_late` | um/s |
| `mean_polarity_projected_speed_reorientation` | um/s |
| `mean_polarity_projected_speed_reversal_cue` | um/s |
| `mean_signed_velocity_cue` | um/s |
| `mean_signed_velocity_early` | um/s |
| `mean_signed_velocity_established` | um/s |
| `mean_signed_velocity_late` | um/s |
| `mean_signed_velocity_reorientation` | um/s |
| `mean_signed_velocity_reversal_cue` | um/s |
| `mechanical_case` | categorical, Boolean, or model units; see METHODS.md |
| `median_elastic_tension` | categorical, Boolean, or model units; see METHODS.md |
| `minimum_length_um` | um |
| `net_displacement_um` | um |
| `nodes` | count |
| `p90_elastic_tension` | categorical, Boolean, or model units; see METHODS.md |
| `p90_tension_added_deactivation_rate` | model rate units; generally s^-1 for chemical rates |
| `protocol` | categorical, Boolean, or model units; see METHODS.md |
| `reference_length_um` | um |
| `spring_stiffness` | categorical, Boolean, or model units; see METHODS.md |
| `tension_feedback_alpha` | dimensionless/model control |
| `terminal_absolute_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_polarity_status` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_polarization_amplitude` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_pole_count` | count |
| `terminal_signed_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `time_to_polarization` | categorical, Boolean, or model units; see METHODS.md |

### `migration1d/results/single_species/testing/focused_motile_length_spatial_tension/focused_motile_length.csv`

Rows: 240. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `chemistry_class` | categorical, Boolean, or model units; see METHODS.md |
| `dashpot_viscosity` | categorical, Boolean, or model units; see METHODS.md |
| `duration_seconds` | s |
| `final_length_um` | um |
| `fraction_time_under_tension` | s unless stated as a sampled index |
| `integrated_tension_feedback_exposure` | categorical, Boolean, or model units; see METHODS.md |
| `maximum_absolute_strain` | dimensionless |
| `maximum_active_gtpase_conservation_error` | dimensionless relative error |
| `maximum_elastic_tension` | categorical, Boolean, or model units; see METHODS.md |
| `maximum_length_um` | um |
| `maximum_tension_added_deactivation_rate` | model rate units; generally s^-1 for chemical rates |
| `maximum_tension_feedback_signal` | categorical, Boolean, or model units; see METHODS.md |
| `mean_absolute_speed_cue` | um/s |
| `mean_absolute_speed_early` | um/s |
| `mean_absolute_speed_established` | um/s |
| `mean_absolute_speed_late` | um/s |
| `mean_absolute_speed_reorientation` | um/s |
| `mean_absolute_speed_reversal_cue` | um/s |
| `mean_polarity_projected_speed_cue` | um/s |
| `mean_polarity_projected_speed_early` | um/s |
| `mean_polarity_projected_speed_established` | um/s |
| `mean_polarity_projected_speed_late` | um/s |
| `mean_polarity_projected_speed_reorientation` | um/s |
| `mean_polarity_projected_speed_reversal_cue` | um/s |
| `mean_signed_velocity_cue` | um/s |
| `mean_signed_velocity_early` | um/s |
| `mean_signed_velocity_established` | um/s |
| `mean_signed_velocity_late` | um/s |
| `mean_signed_velocity_reorientation` | um/s |
| `mean_signed_velocity_reversal_cue` | um/s |
| `mechanical_case` | categorical, Boolean, or model units; see METHODS.md |
| `median_elastic_tension` | categorical, Boolean, or model units; see METHODS.md |
| `median_tension_feedback_signal` | categorical, Boolean, or model units; see METHODS.md |
| `minimum_length_um` | um |
| `net_displacement_um` | um |
| `nodes` | count |
| `p90_elastic_tension` | categorical, Boolean, or model units; see METHODS.md |
| `p90_tension_added_deactivation_rate` | model rate units; generally s^-1 for chemical rates |
| `p90_tension_feedback_signal` | categorical, Boolean, or model units; see METHODS.md |
| `protocol` | categorical, Boolean, or model units; see METHODS.md |
| `reference_length_um` | um |
| `spatial_tension_characteristic_length` | um when field name states physical/reference/scaled length; Lambda is dimensionless |
| `spatial_tension_diffusivity` | categorical, Boolean, or model units; see METHODS.md |
| `spatial_tension_generation_rate` | model rate units; generally s^-1 for chemical rates |
| `spatial_tension_relaxation_time` | s unless stated as a sampled index |
| `spring_stiffness` | categorical, Boolean, or model units; see METHODS.md |
| `tension_feedback_alpha` | dimensionless/model control |
| `tension_feedback_mode` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_absolute_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_polarity_status` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_polarization_amplitude` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_pole_count` | count |
| `terminal_signed_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `time_to_polarization` | categorical, Boolean, or model units; see METHODS.md |

### `migration1d/results/single_species/testing/polarity_reversal/polarity_reversal_map.csv`

Rows: 867. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `mu` | dimensionless/model control |
| `Lambda` | dimensionless/model control |
| `h` | dimensionless/model control |
| `prechallenge_time` | s unless stated as a sampled index |
| `prechallenge_polarized_domain_count` | categorical, Boolean, or model units; see METHODS.md |
| `prechallenge_endpoint_contrast` | dimensionless |
| `terminal_polarized_domain_count` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_endpoint_contrast` | dimensionless |
| `reversal_outcome` | categorical, Boolean, or model units; see METHODS.md |
| `reversal_time_after_challenge` | s unless stated as a sampled index |
| `maximum_conservation_error` | dimensionless relative error |

### `migration1d/results/single_species/testing/wgd_null/double_stimulus/two_seeded_domains.csv`

Rows: 6,084. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `baseline_Lambda` | categorical, Boolean, or model units; see METHODS.md |
| `baseline_h` | categorical, Boolean, or model units; see METHODS.md |
| `baseline_mu` | categorical, Boolean, or model units; see METHODS.md |
| `epsilon` | dimensionless/model control |
| `experiment` | categorical, Boolean, or model units; see METHODS.md |
| `linear_scale` | fold |
| `maximum_conservation_error` | dimensionless relative error |
| `polarized_domain_count_t0` | count |
| `polarized_domain_count_t100` | count |
| `polarized_domain_count_t300` | count |
| `resting_active_u` | categorical, Boolean, or model units; see METHODS.md |
| `resting_inactive_v` | categorical, Boolean, or model units; see METHODS.md |
| `scaled_Lambda` | categorical, Boolean, or model units; see METHODS.md |
| `scaled_h` | categorical, Boolean, or model units; see METHODS.md |
| `scaled_length` | um when field name states physical/reference/scaled length; Lambda is dimensionless |
| `scaled_mu` | categorical, Boolean, or model units; see METHODS.md |
| `volume_multiplier` | fold |

### `migration1d/results/single_species/testing/wgd_null/polarity_reversal/polarity_reversal.csv`

Rows: 2,601. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `baseline_Lambda` | categorical, Boolean, or model units; see METHODS.md |
| `baseline_h` | categorical, Boolean, or model units; see METHODS.md |
| `baseline_mu` | categorical, Boolean, or model units; see METHODS.md |
| `experiment` | categorical, Boolean, or model units; see METHODS.md |
| `linear_scale` | fold |
| `maximum_conservation_error` | dimensionless relative error |
| `prechallenge_endpoint_contrast` | dimensionless |
| `prechallenge_polarized_domain_count` | categorical, Boolean, or model units; see METHODS.md |
| `prechallenge_time` | s unless stated as a sampled index |
| `resting_active_u` | categorical, Boolean, or model units; see METHODS.md |
| `resting_inactive_v` | categorical, Boolean, or model units; see METHODS.md |
| `reversal_outcome` | categorical, Boolean, or model units; see METHODS.md |
| `reversal_time_after_challenge` | s unless stated as a sampled index |
| `scaled_Lambda` | categorical, Boolean, or model units; see METHODS.md |
| `scaled_h` | categorical, Boolean, or model units; see METHODS.md |
| `scaled_length` | um when field name states physical/reference/scaled length; Lambda is dimensionless |
| `scaled_mu` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_endpoint_contrast` | dimensionless |
| `terminal_polarized_domain_count` | categorical, Boolean, or model units; see METHODS.md |
| `volume_multiplier` | fold |

### `migration1d/results/single_species/testing/wgd_null/wgd_null_screen.csv`

Rows: 2,601. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `volume_multiplier` | fold |
| `linear_scale` | fold |
| `baseline_mu` | categorical, Boolean, or model units; see METHODS.md |
| `baseline_Lambda` | categorical, Boolean, or model units; see METHODS.md |
| `baseline_h` | categorical, Boolean, or model units; see METHODS.md |
| `resting_active_u` | categorical, Boolean, or model units; see METHODS.md |
| `resting_inactive_v` | categorical, Boolean, or model units; see METHODS.md |
| `scaled_mu` | categorical, Boolean, or model units; see METHODS.md |
| `scaled_Lambda` | categorical, Boolean, or model units; see METHODS.md |
| `scaled_h` | categorical, Boolean, or model units; see METHODS.md |
| `scaled_length` | um when field name states physical/reference/scaled length; Lambda is dimensionless |
| `terminal_outcome` | categorical, Boolean, or model units; see METHODS.md |
| `polarized_domain_count_t100` | count |
| `polarized_domain_count_t300` | count |
| `terminal_polarization_amplitude` | categorical, Boolean, or model units; see METHODS.md |
| `late_maximum_normalized_change` | categorical, Boolean, or model units; see METHODS.md |
| `maximum_conservation_error` | dimensionless relative error |

### `migration1d/results/three_species/calibration/tier1_regime_map/tier1_regime_map.csv`

Rows: 507. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `volume_multiplier` | fold |
| `linear_scale` | fold |
| `baseline_Lambda` | categorical, Boolean, or model units; see METHODS.md |
| `scaled_Lambda` | categorical, Boolean, or model units; see METHODS.md |
| `baseline_physical_length_um` | um |
| `scaled_physical_length_um` | um |
| `kappa_R` | dimensionless/model control |
| `I_R1_uM_per_s` | categorical, Boolean, or model units; see METHODS.md |
| `n` | dimensionless/model control |
| `homogeneous_growth_rate` | model rate units; generally s^-1 for chemical rates |
| `maximum_spatial_growth_rate` | model rate units; generally s^-1 for chemical rates |
| `fastest_spatial_mode` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_polarity_status` | categorical, Boolean, or model units; see METHODS.md |
| `time_to_polarization` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_polarization_amplitude` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_pole_count` | count |
| `terminal_high_width_fraction` | dimensionless fraction |
| `late_maximum_normalized_change` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_composite_response_status` | categorical, Boolean, or model units; see METHODS.md |
| `composite_late_maximum_normalized_change` | categorical, Boolean, or model units; see METHODS.md |
| `time_to_coherent_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `maximum_coherent_polarity_duration` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_rac_active_fraction` | dimensionless fraction |
| `cdc42_range_contrast` | dimensionless |
| `rac_range_contrast` | dimensionless |
| `rho_range_contrast` | dimensionless |
| `cdc42_rac_correlation` | dimensionless |
| `cdc42_rho_correlation` | dimensionless |
| `cdc42_moment` | categorical, Boolean, or model units; see METHODS.md |
| `rac_moment` | categorical, Boolean, or model units; see METHODS.md |
| `rho_moment` | categorical, Boolean, or model units; see METHODS.md |
| `rac_pole_count` | count |
| `coherent_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_cdc42_endpoint_contrast` | dimensionless |
| `terminal_rac_endpoint_contrast` | dimensionless |
| `terminal_rho_endpoint_contrast` | dimensionless |
| `maximum_cdc42_conservation_error` | dimensionless relative error |
| `maximum_rac_conservation_error` | dimensionless relative error |
| `maximum_rho_conservation_error` | dimensionless relative error |
| `Lambda` | dimensionless/model control |
| `physical_length_um` | um |

### `migration1d/results/three_species/testing/focused_motile_length/focused_motile_length.csv`

Rows: 240. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `chemistry_class` | categorical, Boolean, or model units; see METHODS.md |
| `dashpot_viscosity` | categorical, Boolean, or model units; see METHODS.md |
| `duration_seconds` | s |
| `final_length_um` | um |
| `fraction_time_under_tension` | s unless stated as a sampled index |
| `integrated_tension_feedback_exposure` | categorical, Boolean, or model units; see METHODS.md |
| `maximum_absolute_strain` | dimensionless |
| `maximum_cdc42_conservation_error` | dimensionless relative error |
| `maximum_elastic_tension` | categorical, Boolean, or model units; see METHODS.md |
| `maximum_length_um` | um |
| `maximum_rac_conservation_error` | dimensionless relative error |
| `maximum_rho_conservation_error` | dimensionless relative error |
| `maximum_tension_added_deactivation_rate` | model rate units; generally s^-1 for chemical rates |
| `mean_absolute_speed_cue` | um/s |
| `mean_absolute_speed_early` | um/s |
| `mean_absolute_speed_established` | um/s |
| `mean_absolute_speed_late` | um/s |
| `mean_absolute_speed_reorientation` | um/s |
| `mean_absolute_speed_reversal_cue` | um/s |
| `mean_polarity_projected_speed_cue` | um/s |
| `mean_polarity_projected_speed_early` | um/s |
| `mean_polarity_projected_speed_established` | um/s |
| `mean_polarity_projected_speed_late` | um/s |
| `mean_polarity_projected_speed_reorientation` | um/s |
| `mean_polarity_projected_speed_reversal_cue` | um/s |
| `mean_signed_velocity_cue` | um/s |
| `mean_signed_velocity_early` | um/s |
| `mean_signed_velocity_established` | um/s |
| `mean_signed_velocity_late` | um/s |
| `mean_signed_velocity_reorientation` | um/s |
| `mean_signed_velocity_reversal_cue` | um/s |
| `mechanical_case` | categorical, Boolean, or model units; see METHODS.md |
| `median_elastic_tension` | categorical, Boolean, or model units; see METHODS.md |
| `minimum_length_um` | um |
| `net_displacement_um` | um |
| `nodes` | count |
| `p90_elastic_tension` | categorical, Boolean, or model units; see METHODS.md |
| `p90_tension_added_deactivation_rate` | model rate units; generally s^-1 for chemical rates |
| `protocol` | categorical, Boolean, or model units; see METHODS.md |
| `reference_length_um` | um |
| `spring_stiffness` | categorical, Boolean, or model units; see METHODS.md |
| `tension_feedback_alpha` | dimensionless/model control |
| `terminal_absolute_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_cdc42_moment` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_cdc42_rac_correlation` | dimensionless |
| `terminal_cdc42_range_contrast` | dimensionless |
| `terminal_cdc42_rho_correlation` | dimensionless |
| `terminal_coherent_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_polarity_status` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_polarization_amplitude` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_pole_count` | count |
| `terminal_rac_moment` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_rac_pole_count` | count |
| `terminal_rac_range_contrast` | dimensionless |
| `terminal_rho_moment` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_rho_range_contrast` | dimensionless |
| `terminal_signed_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `time_to_polarization` | categorical, Boolean, or model units; see METHODS.md |

### `migration1d/results/three_species/testing/focused_motile_length_spatial_tension/focused_motile_length.csv`

Rows: 240. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `chemistry_class` | categorical, Boolean, or model units; see METHODS.md |
| `dashpot_viscosity` | categorical, Boolean, or model units; see METHODS.md |
| `duration_seconds` | s |
| `final_length_um` | um |
| `fraction_time_under_tension` | s unless stated as a sampled index |
| `integrated_tension_feedback_exposure` | categorical, Boolean, or model units; see METHODS.md |
| `maximum_absolute_strain` | dimensionless |
| `maximum_cdc42_conservation_error` | dimensionless relative error |
| `maximum_elastic_tension` | categorical, Boolean, or model units; see METHODS.md |
| `maximum_length_um` | um |
| `maximum_rac_conservation_error` | dimensionless relative error |
| `maximum_rho_conservation_error` | dimensionless relative error |
| `maximum_tension_added_deactivation_rate` | model rate units; generally s^-1 for chemical rates |
| `maximum_tension_feedback_signal` | categorical, Boolean, or model units; see METHODS.md |
| `mean_absolute_speed_cue` | um/s |
| `mean_absolute_speed_early` | um/s |
| `mean_absolute_speed_established` | um/s |
| `mean_absolute_speed_late` | um/s |
| `mean_absolute_speed_reorientation` | um/s |
| `mean_absolute_speed_reversal_cue` | um/s |
| `mean_polarity_projected_speed_cue` | um/s |
| `mean_polarity_projected_speed_early` | um/s |
| `mean_polarity_projected_speed_established` | um/s |
| `mean_polarity_projected_speed_late` | um/s |
| `mean_polarity_projected_speed_reorientation` | um/s |
| `mean_polarity_projected_speed_reversal_cue` | um/s |
| `mean_signed_velocity_cue` | um/s |
| `mean_signed_velocity_early` | um/s |
| `mean_signed_velocity_established` | um/s |
| `mean_signed_velocity_late` | um/s |
| `mean_signed_velocity_reorientation` | um/s |
| `mean_signed_velocity_reversal_cue` | um/s |
| `mechanical_case` | categorical, Boolean, or model units; see METHODS.md |
| `median_elastic_tension` | categorical, Boolean, or model units; see METHODS.md |
| `median_tension_feedback_signal` | categorical, Boolean, or model units; see METHODS.md |
| `minimum_length_um` | um |
| `net_displacement_um` | um |
| `nodes` | count |
| `p90_elastic_tension` | categorical, Boolean, or model units; see METHODS.md |
| `p90_tension_added_deactivation_rate` | model rate units; generally s^-1 for chemical rates |
| `p90_tension_feedback_signal` | categorical, Boolean, or model units; see METHODS.md |
| `protocol` | categorical, Boolean, or model units; see METHODS.md |
| `reference_length_um` | um |
| `spatial_tension_characteristic_length` | um when field name states physical/reference/scaled length; Lambda is dimensionless |
| `spatial_tension_diffusivity` | categorical, Boolean, or model units; see METHODS.md |
| `spatial_tension_generation_rate` | model rate units; generally s^-1 for chemical rates |
| `spatial_tension_relaxation_time` | s unless stated as a sampled index |
| `spring_stiffness` | categorical, Boolean, or model units; see METHODS.md |
| `tension_feedback_alpha` | dimensionless/model control |
| `tension_feedback_mode` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_absolute_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_cdc42_moment` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_cdc42_rac_correlation` | dimensionless |
| `terminal_cdc42_range_contrast` | dimensionless |
| `terminal_cdc42_rho_correlation` | dimensionless |
| `terminal_coherent_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_polarity_status` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_polarization_amplitude` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_pole_count` | count |
| `terminal_rac_moment` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_rac_pole_count` | count |
| `terminal_rac_range_contrast` | dimensionless |
| `terminal_rho_moment` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_rho_range_contrast` | dimensionless |
| `terminal_signed_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `time_to_polarization` | categorical, Boolean, or model units; see METHODS.md |

### `migration1d/results/three_species/testing/wgd_null/double_stimulus/two_seeded_domains.csv`

Rows: 6,084. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `baseline_Lambda` | categorical, Boolean, or model units; see METHODS.md |
| `epsilon` | dimensionless/model control |
| `experiment` | categorical, Boolean, or model units; see METHODS.md |
| `kappa_R` | dimensionless/model control |
| `linear_scale` | fold |
| `maximum_cdc42_conservation_error` | dimensionless relative error |
| `maximum_rac_conservation_error` | dimensionless relative error |
| `maximum_rho_conservation_error` | dimensionless relative error |
| `n` | dimensionless/model control |
| `polarized_domain_count_t0` | count |
| `polarized_domain_count_t100` | count |
| `polarized_domain_count_t300` | count |
| `resting_cdc42` | categorical, Boolean, or model units; see METHODS.md |
| `resting_rac` | categorical, Boolean, or model units; see METHODS.md |
| `resting_rho` | categorical, Boolean, or model units; see METHODS.md |
| `scaled_Lambda` | categorical, Boolean, or model units; see METHODS.md |
| `scaled_length_um` | um |
| `terminal_cdc42_rac_correlation` | dimensionless |
| `terminal_cdc42_rho_correlation` | dimensionless |
| `terminal_coherent_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `volume_multiplier` | fold |

### `migration1d/results/three_species/testing/wgd_null/polarity_reversal/polarity_reversal.csv`

Rows: 2,601. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `baseline_Lambda` | categorical, Boolean, or model units; see METHODS.md |
| `experiment` | categorical, Boolean, or model units; see METHODS.md |
| `kappa_R` | dimensionless/model control |
| `linear_scale` | fold |
| `maximum_cdc42_conservation_error` | dimensionless relative error |
| `maximum_rac_conservation_error` | dimensionless relative error |
| `maximum_rho_conservation_error` | dimensionless relative error |
| `n` | dimensionless/model control |
| `prechallenge_endpoint_contrast` | dimensionless |
| `prechallenge_polarized_domain_count` | categorical, Boolean, or model units; see METHODS.md |
| `prechallenge_time` | s unless stated as a sampled index |
| `resting_cdc42` | categorical, Boolean, or model units; see METHODS.md |
| `resting_rac` | categorical, Boolean, or model units; see METHODS.md |
| `resting_rho` | categorical, Boolean, or model units; see METHODS.md |
| `reversal_outcome` | categorical, Boolean, or model units; see METHODS.md |
| `reversal_time_after_challenge` | s unless stated as a sampled index |
| `scaled_Lambda` | categorical, Boolean, or model units; see METHODS.md |
| `scaled_length_um` | um |
| `terminal_coherent_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_endpoint_contrast` | dimensionless |
| `terminal_polarized_domain_count` | categorical, Boolean, or model units; see METHODS.md |
| `volume_multiplier` | fold |

### `migration1d/results/three_species/testing/wgd_null/wgd_null_screen.csv`

Rows: 1,521. Missing numeric entries are stored as empty CSV fields or `None`-like values and generally mean that a quantity is undefined for that case, such as reversal time without a qualifying successful reversal.

| Field | Units/type note |
|---|---|
| `volume_multiplier` | fold |
| `linear_scale` | fold |
| `baseline_Lambda` | categorical, Boolean, or model units; see METHODS.md |
| `scaled_Lambda` | categorical, Boolean, or model units; see METHODS.md |
| `baseline_physical_length_um` | um |
| `scaled_physical_length_um` | um |
| `kappa_R` | dimensionless/model control |
| `I_R1_uM_per_s` | categorical, Boolean, or model units; see METHODS.md |
| `n` | dimensionless/model control |
| `homogeneous_growth_rate` | model rate units; generally s^-1 for chemical rates |
| `maximum_spatial_growth_rate` | model rate units; generally s^-1 for chemical rates |
| `fastest_spatial_mode` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_polarity_status` | categorical, Boolean, or model units; see METHODS.md |
| `time_to_polarization` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_polarization_amplitude` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_pole_count` | count |
| `terminal_high_width_fraction` | dimensionless fraction |
| `late_maximum_normalized_change` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_composite_response_status` | categorical, Boolean, or model units; see METHODS.md |
| `composite_late_maximum_normalized_change` | categorical, Boolean, or model units; see METHODS.md |
| `time_to_coherent_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `maximum_coherent_polarity_duration` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_rac_active_fraction` | dimensionless fraction |
| `cdc42_range_contrast` | dimensionless |
| `rac_range_contrast` | dimensionless |
| `rho_range_contrast` | dimensionless |
| `cdc42_rac_correlation` | dimensionless |
| `cdc42_rho_correlation` | dimensionless |
| `cdc42_moment` | categorical, Boolean, or model units; see METHODS.md |
| `rac_moment` | categorical, Boolean, or model units; see METHODS.md |
| `rho_moment` | categorical, Boolean, or model units; see METHODS.md |
| `rac_pole_count` | count |
| `coherent_polarity` | categorical, Boolean, or model units; see METHODS.md |
| `terminal_cdc42_endpoint_contrast` | dimensionless |
| `terminal_rac_endpoint_contrast` | dimensionless |
| `terminal_rho_endpoint_contrast` | dimensionless |
| `maximum_cdc42_conservation_error` | dimensionless relative error |
| `maximum_rac_conservation_error` | dimensionless relative error |
| `maximum_rho_conservation_error` | dimensionless relative error |

## Summary JSON and manifests

Each completed analysis directory generally contains `summary.json`. These files record design axes, fixed settings, aggregate counts, numerical maxima, and—where relevant—exact matched transition counts. They are small metadata companions to the CSVs rather than substitutes for the condition-level data.

| Type | Files | Total logical size |
|---|---:|---:|
| `summary.json` | 16 | 36.84 KiB |
| `figure_manifest.json` | 5 | 5.43 KiB |

Summary JSON paths:

- `migration1d/results/single_species/calibration/local_bistability/summary.json`
- `migration1d/results/single_species/calibration/regime_map/summary.json`
- `migration1d/results/single_species/exploratory/single_pulse_dual_example/summary.json`
- `migration1d/results/single_species/testing/dual_pole_map/summary.json`
- `migration1d/results/single_species/testing/focused_motile_length/summary.json`
- `migration1d/results/single_species/testing/focused_motile_length_spatial_tension/summary.json`
- `migration1d/results/single_species/testing/polarity_reversal/summary.json`
- `migration1d/results/single_species/testing/wgd_null/double_stimulus/summary.json`
- `migration1d/results/single_species/testing/wgd_null/polarity_reversal/summary.json`
- `migration1d/results/single_species/testing/wgd_null/summary.json`
- `migration1d/results/three_species/calibration/tier1_regime_map/summary.json`
- `migration1d/results/three_species/testing/focused_motile_length/summary.json`
- `migration1d/results/three_species/testing/focused_motile_length_spatial_tension/summary.json`
- `migration1d/results/three_species/testing/wgd_null/double_stimulus/summary.json`
- `migration1d/results/three_species/testing/wgd_null/polarity_reversal/summary.json`
- `migration1d/results/three_species/testing/wgd_null/summary.json`

Figure manifests exist only for the fixed-strip single-species calibration and testing sections and the three-species tier-1 calibration. They contain titles, captions, and associated data filenames. Moving-length summaries do not have equivalent manifests.

## Moving-interval trajectory tier

There are 960 compressed NPZ files, one for each moving-interval trajectory, totaling 841.27 MiB in logical file size. They are organized under the four `focused_motile_length` task directories. Each trajectory has 501 saved times (0 through 500 s at 1 s intervals). Node count is `max(81, ceil(L0/0.25)+1)`, giving 81, 101, 121, 161, or 201 nodes across the tested lengths.

Common NPZ arrays:

| Array | Shape | Meaning |
|---|---|---|
| `times` | `(501,)` | saved times in seconds |
| `position` | `(501, nodes)` | laboratory-frame node positions |
| `membrane` | `(501, species, nodes)` | active fields; species=1 or 3 |
| `center` | `(501,)` | interval center |
| `length` | `(501,)` | interval length |
| `center_velocity` | `(501,)` | mean endpoint velocity |
| `polarity` | `(501,)` | signed polarity definition from METHODS.md |
| `elastic_tension` | `(501,)` | positive elastic spring force |
| `tension_feedback_signal` | `(501,)` or `(501,nodes)` when present | global elastic signal or spatial excess-tension signal |
| `tension_added_deactivation_rate` | `(501,)` or `(501,nodes)` | signal multiplied by alpha_t |

The realized files from the global campaigns may omit the separate `tension_feedback_signal` array while retaining elastic tension and the added deactivation rate. The spatial files retain the full spatial signal and rate.

The corresponding 960 task JSON files store the task index, full design key, merged-row metrics, and NPZ filename. The merged `focused_motile_length.csv` files contain the same row-level summary metrics, so task JSON is mainly useful for linking a summary row to its trajectory.

## WGD intermediate data

The three-species endpoint WGD array retains 1,521 one-condition task JSON files. The single-species endpoint screen uses CSV shards. Across both WGD endpoint screens there are 128 shard files (CSV data plus JSON metadata), totaling 1.30 MiB. The merged endpoint CSVs validate and supersede these for ordinary scientific analysis.

The two extended WGD campaigns retain 17,370 per-condition JSON records. These contain scalar outcomes only, not time-resolved fields. They are merged into the four competition/reversal CSVs cataloged above. The large on-disk footprint of the extended work directories is partly due to filesystem allocation overhead from thousands of small files; their logical content is 11.49 MiB.

Intermediate paths:

- `migration1d/results/single_species/testing/wgd_null/shards/`
- `migration1d/results/three_species/testing/wgd_null/shards/`
- `migration1d/results/three_species/testing/wgd_null/tasks/`
- `migration1d/results/single_species/testing/wgd_null_extended_work/extended_tasks/`
- `migration1d/results/three_species/testing/wgd_null_extended_work/extended_tasks/`

## Existing reports and rendered outputs

| Type | Count | Total logical size | Notes |
|---|---:|---:|---|
| HTML | 5 | 76.21 KiB | Three assembled reports plus two embedded report fragments |
| PNG | 77 | 21.78 MiB | Generated analysis plots and kymograph montages |
| PDF | 1 | 66.67 KiB | Directional atlas PDF |

HTML files:

- `migration1d/results/single_species/report.html`
- `migration1d/results/single_species/testing/wgd_null/extended_report_fragment.html`
- `migration1d/results/single_species/testing/wgd_null/report.html`
- `migration1d/results/three_species/testing/wgd_null/extended_report_fragment.html`
- `migration1d/results/three_species/testing/wgd_null/report.html`

These presentation files are not inputs to the workbook and are not required to access the numerical results.

## Source code and exact experimental definitions

The following modules are the primary audit trail for equations, metrics, and campaign transformations:

- `migration1d/src/migration1d/model.py`: moving finite-volume solver and endpoint mechanics;
- `migration1d/src/migration1d/forces.py`: active endpoint-force mappings;
- `migration1d/src/migration1d/tension.py`: global and spatial tension feedback;
- `migration1d/src/migration1d/single_species/model.py`: wave-pinning chemistry;
- `migration1d/src/migration1d/single_species/observables.py`: active-field classifications;
- `migration1d/src/migration1d/single_species/calibration.py`: local roots and single-cue map;
- `migration1d/src/migration1d/single_species/testing.py`: dual-domain and reversal tests;
- `migration1d/src/migration1d/holmes.py`: six-field Holmes chemistry;
- `migration1d/src/migration1d/three_species/model.py`: mapped Holmes parameters and stability;
- `migration1d/src/migration1d/three_species/observables.py`: composite metrics;
- `migration1d/src/migration1d/single_species/wgd.py` and `wgd_extended.py`: single-species WGD null;
- `migration1d/src/migration1d/three_species/wgd.py` and `wgd_extended.py`: three-species WGD null;
- `migration1d/src/migration1d/focused_motile_length.py`: moving-length design and row metrics;
- `migration1d/src/migration1d/directional_motility_atlas.py`: cumulative directional travel;
- `migration1d/configs/`: campaign selections and scheduler resources.

## Selecting data beyond the workbook

Use the merged CSVs when a requested analysis depends only on parameter conditions, terminal classifications, recorded time-point domain counts, precomputed interval summaries, or numerical diagnostics. Use moving-task NPZ files when the analysis requires a different time window, a spatial profile, laboratory-frame positions, time-resolved length or velocity, or recomputation under a different outcome threshold. WGD extended runs cannot be reconstructed as full trajectories from retained task JSON; rerunning the recorded campaign would be required for unrecorded time-resolved quantities.
