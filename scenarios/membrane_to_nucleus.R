# Scenario: membrane-generated cAMP-like signal competing for two sinks.

membrane_to_nucleus_scenario <- function() {
  list(
    name = "Membrane-to-nucleus signalling with distributed degradation",
    description = paste(
      "A plasma-membrane source supplies a cAMP-like signal to a spherical",
      "cytoplasmic shell. Distributed Michaelis-Menten PDE degradation",
      "competes with productive capture at a reactive nuclear envelope."
    ),
    geometry = list(
      cell_radius_um = 9.34,
      nucleus_radius_um = 5.26
    ),
    boundary = list(
      membrane_flux_uM_um_s = 0.3611,
      nuclear_capture_permeability_um_s = 0.03
    ),
    parameters = list(
      initial_concentration_uM = 0.05,
      diffusion_um2_s = c(300, 30, 3),
      pde_vmax_uM_s = 0.295,
      pde_km_uM = 2,
      duration_s = 300,
      nuclear_capture_sweep_um_s = 10 ^ seq(log10(0.003), log10(0.3), length.out = 7)
    )
  )
}
