# Helpers for the Sedlack-inspired, homogenized single-cell oxygen model.

#' Convert oxygen partial pressure (mmHg) to dissolved concentration (mol/m^3).
oxygen_concentration <- function(po2_mmhg, solubility_mol_m3_mmhg) {
  po2_mmhg * solubility_mol_m3_mmhg
}

#' Convert dissolved oxygen concentration (mol/m^3) to partial pressure (mmHg).
oxygen_partial_pressure <- function(concentration_mol_m3, solubility_mol_m3_mmhg) {
  concentration_mol_m3 / solubility_mol_m3_mmhg
}

#' Homogenize a mitochondrial maximum consumption rate over cell volume.
homogenized_vmax <- function(mitochondrial_volume_fraction, vmax_mito_mol_m3_s) {
  stopifnot(mitochondrial_volume_fraction >= 0, vmax_mito_mol_m3_s >= 0)
  mitochondrial_volume_fraction * vmax_mito_mol_m3_s
}

#' Uniform Michaelis-Menten oxygen-consumption field (mol/m^3/s).
oxygen_consumption_mm <- function(concentration_mol_m3, vmax_mol_m3_s, km_mol_m3) {
  stopifnot(vmax_mol_m3_s >= 0, km_mol_m3 > 0)
  vmax_mol_m3_s * concentration_mol_m3 / (km_mol_m3 + concentration_mol_m3)
}

#' Solve a coupled spherical cell/extracellular oxygen model to steady state.
#'
#' This is a fast radial screening model, not the final adherent-cell (r,z)
#' geometry. The membrane flux is coupled to the concentrations on both sides
#' through a finite permeability and includes the half-cell diffusion resistance
#' of each adjacent finite-volume cell.
simulate_coupled_spherical_oxygen <- function(cell_radius_um, outer_radius_um,
                                               bulk_po2_mmhg, parameters,
                                               n_cell = 100L, n_medium = 160L,
                                               duration_s = 1000,
                                               method = "lsodes") {
  if (!exists("radial_grid", mode = "function") ||
      !requireNamespace("ReacTran", quietly = TRUE) ||
      !requireNamespace("deSolve", quietly = TRUE)) {
    stop("Source R/diffusion.R and install ReacTran and deSolve before running this model.")
  }
  stopifnot(outer_radius_um > cell_radius_um, cell_radius_um > 0)
  cell_grid <- radial_grid(cell_radius_um, n_cell)
  medium_grid <- radial_grid(outer_radius_um, n_medium, inner_radius_um = cell_radius_um)
  transport_parts <- function(grid) {
    n <- length(grid$centers_um)
    list(
      grid = list(dx = rep(grid$dr_um, n),
                  dx.aux = c(grid$dr_um / 2, rep(grid$dr_um, n - 1L), grid$dr_um / 2)),
      area = list(int = grid$face_areas_um2, mid = grid$volumes_um3 / grid$dr_um)
    )
  }
  cell_transport <- transport_parts(cell_grid)
  medium_transport <- transport_parts(medium_grid)
  d_cell <- parameters$transport$diffusion_cytoplasm_um2_s
  d_medium <- parameters$transport$diffusion_medium_um2_s
  s_cell <- parameters$transport$solubility_cytoplasm_mol_m3_mmhg
  s_medium <- parameters$transport$solubility_medium_mol_m3_mmhg
  permeability <- parameters$transport$membrane_permeability_um_s
  c_bulk <- oxygen_concentration(bulk_po2_mmhg, s_medium)
  n_cell <- length(cell_grid$centers_um)
  n_medium <- length(medium_grid$centers_um)

  membrane_flux <- function(cell_concentration, medium_concentration) {
    # Positive J is transport from medium into cell. Surface concentrations are
    # reconstructed from adjacent cell-centre concentrations and flux.
    membrane_resistance <- (cell_grid$dr_um / 2) / (d_cell * s_cell) +
      (medium_grid$dr_um / 2) / (d_medium * s_medium)
    driving_po2 <- medium_concentration[1] / s_medium -
      cell_concentration[n_cell] / s_cell
    permeability * s_medium * driving_po2 /
      (1 + permeability * s_medium * membrane_resistance)
  }

  rhs <- function(time_s, state, parameters_unused) {
    cell_concentration <- state[seq_len(n_cell)]
    medium_concentration <- state[n_cell + seq_len(n_medium)]
    influx <- membrane_flux(cell_concentration, medium_concentration)
    cell_transport_term <- ReacTran::tran.1D(
      C = cell_concentration, D = d_cell, flux.down = -influx,
      A = cell_transport$area, dx = cell_transport$grid
    )$dC
    medium_transport_term <- ReacTran::tran.1D(
      C = medium_concentration, D = d_medium, flux.up = -influx,
      # The outer edge of this screening shell represents Sedlack's maintained
      # medium boundary, placed one medium-column height from the cell surface.
      C.down = c_bulk,
      A = medium_transport$area, dx = medium_transport$grid
    )$dC
    consumption <- oxygen_consumption_mm(
      pmax(cell_concentration, 0), parameters$consumption$vmax_effective_mol_m3_s,
      parameters$consumption$km_mol_m3
    )
    list(c(cell_transport_term - consumption, medium_transport_term))
  }

  solution <- deSolve::ode(
    y = rep(c_bulk, n_cell + n_medium), times = c(0, duration_s),
    func = rhs, parms = NULL, method = method, rtol = 1e-7, atol = 1e-10
  )
  final_state <- solution[nrow(solution), -1]
  final_cell <- final_state[seq_len(n_cell)]
  final_medium <- final_state[n_cell + seq_len(n_medium)]
  consumption <- oxygen_consumption_mm(
    pmax(final_cell, 0), parameters$consumption$vmax_effective_mol_m3_s,
    parameters$consumption$km_mol_m3
  )
  list(
    cell_grid = cell_grid, medium_grid = medium_grid,
    cell_concentration_mol_m3 = final_cell,
    medium_concentration_mol_m3 = final_medium,
    cell_po2_mmhg = oxygen_partial_pressure(final_cell, s_cell),
    medium_po2_mmhg = oxygen_partial_pressure(final_medium, s_medium),
    consumption_mol_m3_s = consumption,
    consumption_per_cell_volume_mol_m3_s = volume_weighted_mean(consumption, cell_grid),
    membrane_flux_mol_m2_s = membrane_flux(final_cell, final_medium) * 1e-6
  )
}
