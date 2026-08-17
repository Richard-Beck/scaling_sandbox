# Minimal 3D Cartesian oxygen model for an adherent half-ellipsoidal cell.

#' Build a uniform Cartesian medium box containing a half-ellipsoidal cell.
build_oxygen_grid_3d <- function(spacing_um, medium_height_um = 150,
                                 cell_axes_um = c(x = 30, y = 30, z = 12),
                                 dx_um = 5, dy_um = dx_um, dz_um = 5) {
  stopifnot(spacing_um > 2 * max(cell_axes_um[c("x", "y")]),
            medium_height_um > cell_axes_um[["z"]])
  nx <- round(spacing_um / dx_um)
  ny <- round(spacing_um / dy_um)
  nz <- round(medium_height_um / dz_um)
  dx_um <- spacing_um / nx
  dy_um <- spacing_um / ny
  dz_um <- medium_height_um / nz
  x <- (seq_len(nx) - 0.5) * dx_um - spacing_um / 2
  y <- (seq_len(ny) - 0.5) * dy_um - spacing_um / 2
  z <- (seq_len(nz) - 0.5) * dz_um
  coordinates <- expand.grid(x_um = x, y_um = y, z_um = z)
  is_cell <- with(coordinates,
    (x_um / cell_axes_um[["x"]]) ^ 2 +
    (y_um / cell_axes_um[["y"]]) ^ 2 +
    (z_um / cell_axes_um[["z"]]) ^ 2 <= 1
  )
  index <- array(seq_len(nx * ny * nz), dim = c(nx, ny, nz))
  list(
    nx = nx, ny = ny, nz = nz, x_um = x, y_um = y, z_um = z,
    dx_um = dx_um, dy_um = dy_um, dz_um = dz_um,
    coordinates = coordinates, is_cell = is_cell, index = index,
    voxel_volume_um3 = dx_um * dy_um * dz_um,
    represented_cell_volume_um3 = sum(is_cell) * dx_um * dy_um * dz_um,
    analytic_cell_volume_um3 = 2 / 3 * pi * prod(cell_axes_um),
    spacing_um = spacing_um, medium_height_um = medium_height_um,
    cell_axes_um = cell_axes_um
  )
}

#' Assemble the sparse finite-volume diffusion operator.
#'
#' The current Sedlack baseline uses equal O2 solubility in medium and cell.
#' At cell-medium faces, membrane resistance 1/P is added in series with the
#' two half-voxel diffusion resistances. All other exterior faces are no-flux,
#' except the top z face, which is fixed at the applied pO2.
assemble_oxygen_transport_3d <- function(grid, parameters, applied_po2_mmhg) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The 3D oxygen model requires the Matrix package.")
  }
  s_medium <- parameters$transport$solubility_medium_mol_m3_mmhg
  s_cell <- parameters$transport$solubility_cytoplasm_mol_m3_mmhg
  if (!isTRUE(all.equal(s_medium, s_cell))) {
    stop("The initial 3D implementation currently requires equal medium and cell solubility.")
  }
  n <- nrow(grid$coordinates)
  diffusivity <- ifelse(
    grid$is_cell,
    parameters$transport$diffusion_cytoplasm_um2_s,
    parameters$transport$diffusion_medium_um2_s
  )
  permeability <- parameters$transport$membrane_permeability_um_s
  all_i <- integer()
  all_j <- integer()
  all_w <- numeric()

  add_faces <- function(i, j, distance_um) {
    i <- as.integer(i)
    j <- as.integer(j)
    membrane <- xor(grid$is_cell[i], grid$is_cell[j])
    conductance_um_s <- 1 / (
      distance_um / (2 * diffusivity[i]) +
      distance_um / (2 * diffusivity[j]) +
      ifelse(membrane, 1 / permeability, 0)
    )
    list(i = i, j = j, w = conductance_um_s / distance_um)
  }
  if (grid$nx > 1) {
    face <- add_faces(grid$index[seq_len(grid$nx - 1), , , drop = FALSE],
                      grid$index[2:grid$nx, , , drop = FALSE], grid$dx_um)
    all_i <- c(all_i, face$i); all_j <- c(all_j, face$j); all_w <- c(all_w, face$w)
  }
  if (grid$ny > 1) {
    face <- add_faces(grid$index[, seq_len(grid$ny - 1), , drop = FALSE],
                      grid$index[, 2:grid$ny, , drop = FALSE], grid$dy_um)
    all_i <- c(all_i, face$i); all_j <- c(all_j, face$j); all_w <- c(all_w, face$w)
  }
  if (grid$nz > 1) {
    face <- add_faces(grid$index[, , seq_len(grid$nz - 1), drop = FALSE],
                      grid$index[, , 2:grid$nz, drop = FALSE], grid$dz_um)
    all_i <- c(all_i, face$i); all_j <- c(all_j, face$j); all_w <- c(all_w, face$w)
  }
  diagonal <- -as.numeric(Matrix::sparseMatrix(
    i = c(all_i, all_j), j = rep(1L, 2 * length(all_i)),
    x = c(all_w, all_w), dims = c(n, 1)
  ))
  operator <- Matrix::sparseMatrix(
    i = c(all_i, all_j, seq_len(n)), j = c(all_j, all_i, seq_len(n)),
    x = c(all_w, all_w, diagonal), dims = c(n, n)
  )
  top <- as.integer(grid$index[, , grid$nz])
  top_weight <- 2 * diffusivity[top] / grid$dz_um ^ 2
  operator <- operator + Matrix::sparseMatrix(
    i = top, j = top, x = -top_weight, dims = c(n, n)
  )
  bulk_concentration <- oxygen_concentration(applied_po2_mmhg, s_medium)
  boundary_source <- numeric(n)
  boundary_source[top] <- top_weight * bulk_concentration
  diffusion_factor <- Matrix::Cholesky(Matrix::forceSymmetric(-operator))
  list(operator = operator, diffusion_factor = diffusion_factor,
       boundary_source = boundary_source,
       diffusivity_um2_s = diffusivity, bulk_concentration_mol_m3 = bulk_concentration)
}

#' Solve the nonlinear 3D steady state by damped Picard iteration.
solve_oxygen_steady_3d <- function(grid, transport, parameters,
                                   vmax_effective_mol_m3_s =
                                     parameters$consumption$vmax_effective_mol_m3_s,
                                   initial = NULL, tolerance = 1e-9,
                                   max_iterations = 100L) {
  n <- nrow(grid$coordinates)
  if (is.null(initial)) initial <- rep(transport$bulk_concentration_mol_m3, n)
  concentration <- pmax(as.numeric(initial), 0)
  cell_indices <- which(grid$is_cell)
  km <- parameters$consumption$km_mol_m3
  residual_at <- function(current) {
    positive_cell <- pmax(current[cell_indices], 0)
    consumption <- numeric(n)
    consumption[cell_indices] <- vmax_effective_mol_m3_s * positive_cell / (km + positive_cell)
    as.numeric(transport$operator %*% current) + transport$boundary_source - consumption
  }
  converged <- FALSE
  for (iteration in seq_len(max_iterations)) {
    residual <- residual_at(concentration)
    if (max(abs(residual)) < tolerance) {
      converged <- TRUE
      break
    }
    positive_cell <- pmax(concentration[cell_indices], 0)
    consumption <- numeric(n)
    consumption[cell_indices] <- vmax_effective_mol_m3_s * positive_cell /
      (km + positive_cell)
    proposed <- as.numeric(Matrix::solve(
      transport$diffusion_factor, transport$boundary_source - consumption
    ))
    step <- proposed - concentration
    alpha <- 1
    if (any(step < 0)) {
      alpha <- min(alpha, 0.95 * min(-concentration[step < 0] / step[step < 0]))
    }
    baseline_norm <- max(abs(residual))
    while (alpha > 1e-6 &&
           max(abs(residual_at(concentration + alpha * step))) >= baseline_norm) {
      alpha <- alpha / 2
    }
    concentration <- pmax(concentration + alpha * step, 0)
  }
  newton_iteration <- 0L
  if (!converged) {
    # Picard iteration can stall close to an anoxic front. Finish with sparse
    # Newton steps only in that harder regime.
    for (newton_iteration in seq_len(20L)) {
      residual <- residual_at(concentration)
      if (max(abs(residual)) < tolerance) {
        converged <- TRUE
        break
      }
      derivative <- numeric(n)
      positive <- concentration[cell_indices] > 0
      derivative[cell_indices[positive]] <- vmax_effective_mol_m3_s * km /
        (km + concentration[cell_indices[positive]]) ^ 2
      jacobian <- transport$operator - Matrix::Diagonal(n, derivative)
      step <- as.numeric(Matrix::solve(-jacobian, residual))
      alpha <- 1
      if (any(step < 0)) {
        alpha <- min(alpha, 0.95 * min(-concentration[step < 0] / step[step < 0]))
      }
      baseline_norm <- max(abs(residual))
      while (alpha > 1e-8 &&
             max(abs(residual_at(concentration + alpha * step))) >= baseline_norm) {
        alpha <- alpha / 2
      }
      concentration <- pmax(concentration + alpha * step, 0)
    }
    residual <- residual_at(concentration)
    converged <- max(abs(residual)) < tolerance
  }
  cell_concentration <- concentration[cell_indices]
  cell_consumption <- oxygen_consumption_mm(
    cell_concentration, vmax_effective_mol_m3_s, km
  )
  list(
    concentration_mol_m3 = concentration,
    po2_mmhg = oxygen_partial_pressure(
      concentration, parameters$transport$solubility_medium_mol_m3_mmhg
    ),
    consumption_mol_m3_s = replace(numeric(n), cell_indices, cell_consumption),
    mean_cell_po2_mmhg = mean(oxygen_partial_pressure(
      cell_concentration, parameters$transport$solubility_cytoplasm_mol_m3_mmhg
    )),
    mean_cell_consumption_mol_m3_s = mean(cell_consumption),
    total_ocr_fmol_s = mean(cell_consumption) * grid$represented_cell_volume_um3 * 1e-3,
    converged = converged,
    iterations = iteration + newton_iteration,
    max_residual_mol_m3_s = max(abs(residual_at(concentration))), grid = grid
  )
}

#' Calibrate effective Vmax to a target mean cytoplasmic pO2 in the 3D model.
calibrate_oxygen_vmax_3d <- function(grid, transport, parameters,
                                     target_mean_po2_mmhg,
                                     bracket_mol_m3_s = c(0, 2),
                                     tolerance_mol_m3_s = 1e-4) {
  objective <- function(vmax) {
    fit <- solve_oxygen_steady_3d(grid, transport, parameters,
                                  vmax_effective_mol_m3_s = vmax)
    fit$mean_cell_po2_mmhg - target_mean_po2_mmhg
  }
  root <- stats::uniroot(objective, interval = bracket_mol_m3_s,
                         tol = tolerance_mol_m3_s)
  root$root
}

#' Reach difficult hypoxic steady states by continuation in effective Vmax.
solve_oxygen_continuation_3d <- function(grid, transport, parameters,
                                         vmax_effective_mol_m3_s,
                                         steps = 20L) {
  initial <- NULL
  fit <- NULL
  for (vmax in seq(vmax_effective_mol_m3_s / steps,
                   vmax_effective_mol_m3_s, length.out = steps)) {
    fit <- solve_oxygen_steady_3d(
      grid, transport, parameters, vmax_effective_mol_m3_s = vmax,
      initial = initial
    )
    if (!fit$converged) stop("3D oxygen continuation failed at Vmax = ", signif(vmax, 4))
    initial <- fit$concentration_mol_m3
  }
  fit
}
