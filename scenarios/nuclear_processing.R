# Scenario: nuclear production with saturable processing in the cytoplasm.

nuclear_processing_scenario <- function() {
  list(
    name = "Nuclear export with saturable consumption and time-varying degradation",
    description = paste(
      "A nuclear source exports a volume-scaled amount to a cytoplasm with",
      "saturable consumption and a separate first-order degradation term."
    ),
    geometry = list(cell_radius_um = 10, nucleus_radius_fraction = 0.45),
    boundary = list(type = "no_flux"),
    parameters = list(
      diffusion_um2_s = 20,
      # Export per cytoplasmic volume (uM/s); total export is proportional to volume.
      source_rate = 0.03,
      source_exponent = 1,
      consumption_vmax_uM_s = 0.06,
      consumption_km_uM = 0.05,
      capacity_exponent = 1,
      initial_degradation_s = 0.005,
      final_degradation_s = 0.04,
      degradation_transition_s = 75
    )
  )
}
