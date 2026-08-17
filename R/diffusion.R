# One-dimensional radial finite-volume helpers for spherical cells.

#' Make a cell-centred radial grid for a sphere or cytoplasmic shell.
radial_grid <- function(radius_um, n_cells = 100L, inner_radius_um = 0) {
  stopifnot(radius_um > inner_radius_um, inner_radius_um >= 0, n_cells >= 3)
  edges <- seq(inner_radius_um, radius_um, length.out = n_cells + 1L)
  centers <- (edges[-1] + edges[-length(edges)]) / 2
  list(
    edges_um = edges,
    centers_um = centers,
    dr_um = (radius_um - inner_radius_um) / n_cells,
    volumes_um3 = 4 / 3 * pi * diff(edges ^ 3),
    face_areas_um2 = 4 * pi * edges ^ 2
  )
}

#' Integrate radial diffusion with ReacTran and deSolve.
#'
#' This is a method-of-lines model: `ReacTran::tran.1D()` performs the
#' finite-volume transport calculation and `deSolve::ode()` chooses adaptive,
#' stiff-capable time steps. Boundary fluxes are expressed per unit boundary area. An outer
#' `fixed_inward_flux` adds material to the cell; an inner `fixed_outward_flux`
#' exports material from the nucleus to the cytoplasm. An inner `reactive_sink`
#' removes material with permeability `value` (um/s). Its Robin condition is
#' evaluated at the boundary face, including the half-cell diffusion resistance.
#' `loss_s` may be a scalar or a function of time in seconds. Units must be
#' mutually consistent, e.g. um, seconds, and uM.
simulate_radial_diffusion <- function(grid, diffusion_um2_s, initial,
                                      duration_s, source = 0, demand = 0,
                                      loss_s = 0, reaction = NULL,
                                      inner_boundary = list(type = "no_flux"),
                                      outer_boundary = list(type = "no_flux"),
                                      save_times = NULL,
                                      method = "lsodes") {
  stopifnot(diffusion_um2_s > 0, duration_s >= 0)
  if (!requireNamespace("ReacTran", quietly = TRUE) ||
      !requireNamespace("deSolve", quietly = TRUE)) {
    stop("simulate_radial_diffusion() requires the ReacTran and deSolve packages.")
  }
  n <- length(grid$centers_um)
  initial <- rep_len(initial, n)
  source <- rep_len(source, n)
  demand <- rep_len(demand, n)
  return_trajectory <- !is.null(save_times)
  times <- sort(unique(c(0, save_times, duration_s)))
  dx_aux <- c(grid$dr_um / 2, rep(grid$dr_um, n - 1L), grid$dr_um / 2)
  transport_grid <- list(dx = rep(grid$dr_um, n), dx.aux = dx_aux)
  transport_area <- list(
    int = grid$face_areas_um2,
    # A_mid * dx is exactly the volume of each spherical shell.
    mid = grid$volumes_um3 / grid$dr_um
  )

  derivative <- function(time_s, concentration, parameters) {
    flux_up <- if (identical(inner_boundary$type, "fixed_outward_flux")) {
      inner_boundary$value
    } else if (identical(inner_boundary$type, "reactive_sink")) {
      permeability_um_s <- inner_boundary$value
      if (length(permeability_um_s) != 1L || !is.finite(permeability_um_s) ||
          permeability_um_s < 0) {
        stop("reactive_sink permeability must be one non-negative finite value")
      }
      boundary_concentration <- concentration[1] /
        (1 + permeability_um_s * grid$dr_um / (2 * diffusion_um2_s))
      -permeability_um_s * boundary_concentration
    } else if (identical(inner_boundary$type, "no_flux")) {
      NULL
    } else {
      stop("Unsupported inner boundary type: ", inner_boundary$type)
    }
    flux_down <- if (identical(outer_boundary$type, "fixed_inward_flux")) {
      -outer_boundary$value
    } else if (identical(outer_boundary$type, "no_flux") ||
               identical(outer_boundary$type, "infinite_reservoir")) {
      NULL
    } else {
      NULL
    }
    boundary_value <- if (identical(outer_boundary$type, "fixed_concentration") ||
                          identical(outer_boundary$type, "infinite_reservoir")) {
      outer_boundary$value
    } else if (!identical(outer_boundary$type, "fixed_inward_flux") &&
               !identical(outer_boundary$type, "no_flux")) {
      stop("Unsupported outer boundary type: ", outer_boundary$type)
    } else {
      concentration[n]
    }
    transport <- ReacTran::tran.1D(
      C = concentration, C.down = boundary_value,
      flux.up = flux_up, flux.down = flux_down, D = diffusion_um2_s,
      A = transport_area, dx = transport_grid,
      # For an exterior sphere, this implements C -> C_bulk as r -> infinity.
      a.bl.down = if (identical(outer_boundary$type, "infinite_reservoir")) {
        diffusion_um2_s / tail(grid$edges_um, 1)
      } else NULL
    )
    degradation <- if (is.function(loss_s)) loss_s(time_s) else loss_s
    degradation <- rep_len(degradation, n)
    if (any(degradation < 0)) stop("loss_s must be non-negative")
    reaction_rate <- if (is.null(reaction)) 0 else reaction(time_s, concentration)
    reaction_rate <- rep_len(reaction_rate, n)
    if (any(reaction_rate < 0)) stop("reaction rates must be non-negative")
    list(transport$dC + source - demand - degradation * concentration - reaction_rate)
  }
  solution <- deSolve::ode(
    y = initial, times = times, func = derivative, parms = NULL, method = method,
    rtol = 1e-7, atol = 1e-9
  )
  concentrations <- solution[, -1, drop = FALSE]
  if (!return_trajectory) return(as.numeric(concentrations[nrow(concentrations), ]))
  list(time_s = solution[, 1], concentration_uM = concentrations, grid = grid)
}

#' Volume-weighted mean concentration.
volume_weighted_mean <- function(values, grid) {
  weighted.mean(values, grid$volumes_um3)
}

#' Concentration and removal flux at a reactive inner boundary.
#'
#' Converts the first cell-centre concentration to the boundary-face value for
#' the same Robin discretization used by `simulate_radial_diffusion()`.
reactive_inner_boundary <- function(first_cell_concentration_uM, grid,
                                    diffusion_um2_s, permeability_um_s) {
  stopifnot(diffusion_um2_s > 0, permeability_um_s >= 0)
  boundary_concentration <- first_cell_concentration_uM /
    (1 + permeability_um_s * grid$dr_um / (2 * diffusion_um2_s))
  list(
    concentration_uM = boundary_concentration,
    flux_uM_um_s = permeability_um_s * boundary_concentration
  )
}
