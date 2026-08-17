# Scenario: extracellular substrate supplied through the plasma membrane.

external_supply_scenario <- function() {
  list(
    name = "External supply with a fixed influx per cell volume",
    description = paste(
      "Total substrate influx scales with cell volume; the required membrane",
      "flux density therefore increases as the cell gets larger."
    ),
    geometry = list(cell_radius_um = 10, nucleus_radius_fraction = 0.45),
    boundary = list(type = "fixed_inward_flux"),
    parameters = list(
      diffusion_um2_s = 20,
      # Input per cell volume (uM/s); total influx is proportional to volume.
      input_rates_uM_s = c(low = 0.01, high = 0.05),
      consumers = list(
        first_order = list(label = "First-order consumption", loss_s = 0.05),
        michaelis_menten = list(label = "Saturable consumption", vmax_uM_s = 0.10, km_uM = 0.05)
      ),
      external_medium_radius_um = 100,
      external_bulk_concentration_uM = 1,
      source_exponent = 1,
      capacity_exponent = 1
    )
  )
}
