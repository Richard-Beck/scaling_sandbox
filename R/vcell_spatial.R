source("R/vcell_vcml.R")

build_vcell_system <- function(path = "scenarios/wgd_spatial_scale.vcml",
                               simulation = "simulation2",
                               linear_scale = 1, physical_extent = NULL,
                               mesh_size = NULL) {
  stopifnot(requireNamespace("Matrix", quietly = TRUE),
            requireNamespace("Rcpp", quietly = TRUE),
            is.numeric(linear_scale), length(linear_scale) == 1L,
            is.finite(linear_scale), linear_scale > 0)
  model <- read_vcml_model(path, simulation)
  if (!is.null(mesh_size)) {
    stopifnot(is.numeric(mesh_size), length(mesh_size) %in% c(1L, 2L),
              all(is.finite(mesh_size)), all(mesh_size >= 3),
              all(mesh_size == as.integer(mesh_size)))
    if (length(mesh_size) == 1L) mesh_size <- rep(mesh_size, 2L)
    model$mesh <- setNames(as.integer(mesh_size), c("x", "y"))
  }
  baseline_extent <- model$geometry$extent
  if (is.null(physical_extent)) {
    model$geometry$extent <- baseline_extent * linear_scale
  } else {
    stopifnot(is.numeric(physical_extent), length(physical_extent) == 2L,
              all(is.finite(physical_extent)), all(physical_extent > 0))
    if (is.null(names(physical_extent))) names(physical_extent) <- c("x", "y")
    model$geometry$extent <- physical_extent[c("x", "y")]
  }
  axis_scales <- model$geometry$extent / baseline_extent
  mask <- resample_vcell_mask(model)
  nx <- nrow(mask); ny <- ncol(mask)
  dx <- model$geometry$extent[["x"]] / nx
  dy <- model$geometry$extent[["y"]] / ny
  cell_index <- matrix(NA_integer_, nx, ny)
  cell_index[mask] <- seq_len(sum(mask))
  coords <- which(mask, arr.ind = TRUE)
  x <- model$geometry$origin[["x"]] + (coords[, 1] - 0.5) * dx
  y <- model$geometry$origin[["y"]] + (coords[, 2] - 0.5) * dy

  rows <- integer(); cols <- integer(); vals <- numeric()
  add_edges <- function(offset, weight) {
    for (k in seq_len(nrow(coords))) {
      neighbor <- coords[k, ] + offset
      if (neighbor[1] >= 1L && neighbor[1] <= nx &&
          neighbor[2] >= 1L && neighbor[2] <= ny && mask[neighbor[1], neighbor[2]]) {
        rows <<- c(rows, k, k)
        cols <<- c(cols, cell_index[neighbor[1], neighbor[2]], k)
        vals <<- c(vals, weight, -weight)
      }
    }
  }
  add_edges(c(1L, 0L), 1 / dx^2)
  add_edges(c(-1L, 0L), 1 / dx^2)
  add_edges(c(0L, 1L), 1 / dy^2)
  add_edges(c(0L, -1L), 1 / dy^2)
  laplacian <- Matrix::sparseMatrix(i = rows, j = cols, x = vals,
                                    dims = rep(sum(mask), 2L))

  constants <- model$constants
  parameter_order <- c("n", "a1", "a2", "alpha", "Beta", "cdc42_tot",
    "Rac_tot", "Rho_tot", "d_cdc42", "d_rac", "d_rho", "f", "Ip", "Ip1",
    "Ir", "k21", "k_PI3K", "k_PI5K", "k_PTEN", "P3b", "Rb", "Rho_b",
    "delta_P1")
  reaction_parameters <- unname(constants[parameter_order])
  if (anyNA(reaction_parameters)) stop("A required VCML parameter is missing")
  state <- matrix(rep(model$initial, each = sum(mask)), nrow = sum(mask),
                  dimnames = list(NULL, names(model$initial)))
  stimulus_x <- model$geometry$origin[["x"]] +
    (x - model$geometry$origin[["x"]]) / axis_scales[["x"]]
  list(model = model, mask = mask, cell_index = cell_index,
       coords = data.frame(index = seq_len(sum(mask)), x = x, y = y,
                           ix = coords[, 1], iy = coords[, 2]),
       dx = dx, dy = dy, laplacian = laplacian, initial_state = state,
       reaction_parameters = reaction_parameters, stimulus_x = stimulus_x,
       linear_scale = if (is.null(physical_extent)) linear_scale else NA_real_,
       axis_scales = axis_scales, baseline_extent = baseline_extent)
}

compile_vcell_kernel <- function(rebuild = FALSE) {
  result <- tryCatch({
    Rcpp::sourceCpp("R/vcell_model.cpp", rebuild = rebuild, showOutput = FALSE,
                    verbose = FALSE)
    TRUE
  }, error = function(error) {
    message("Rcpp compilation unavailable; using vectorized R reaction kernel.")
    FALSE
  })
  assign(".vcell_kernel_attempted", TRUE, envir = .GlobalEnv)
  result
}

reaction_rates_r <- function(state, x, time, p) {
  n <- p[1]; a1 <- p[2]; a2 <- p[3]; alpha <- p[4]; beta <- p[5]
  ctot <- p[6]; ractot <- p[7]; rhotot <- p[8]
  dcdc <- p[9]; drac <- p[10]; drho <- p[11]; f <- p[12]
  ip <- p[13]; ip1 <- p[14]; ir <- p[15]; k21 <- p[16]
  kpi3k <- p[17]; kpi5k <- p[18]; kpten <- p[19]
  p3b <- p[20]; rb <- p[21]; rhob <- p[22]; delta_p1 <- p[23]
  ic <- if (time <= 10) 2.6 + 0.05 * x else rep(2.95, length(x))
  feedback <- (1 - f) + f * state[, 9] / p3b
  j0 <- ic / (1 + (state[, 6] / a1)^n) * (state[, 1] / ctot) * feedback -
    dcdc * state[, 2]
  j1 <- (ir + alpha * state[, 2]) * (state[, 3] / ractot) * feedback -
    drac * state[, 4]
  j2 <- (ip + beta * state[, 4]) / (1 + (state[, 2] / a2)^n) *
    (state[, 5] / rhotot) - drho * state[, 6]
  j3 <- -k21 * state[, 8] + kpi5k / 2 * (1 + state[, 4] / rb) * state[, 7]
  j4 <- kpi3k / 2 * (1 + state[, 4] / rb) * state[, 8] -
    kpten / 2 * (1 + state[, 6] / rhob) * state[, 9]
  j5 <- ip1 - delta_p1 * state[, 7]
  list(rate = cbind(-j0, j0, -j1, j1, -j2, j2,
                    -j3 + j5, j3 - j4, j4), mean_j5 = mean(j5))
}

reaction_step_r <- function(state, x, time, dt, parameters) {
  k1 <- reaction_rates_r(state, x, time, parameters)
  k2 <- reaction_rates_r(state + 0.5 * dt * k1$rate, x, time + 0.5 * dt,
                         parameters)
  k3 <- reaction_rates_r(state + 0.5 * dt * k2$rate, x, time + 0.5 * dt,
                         parameters)
  k4 <- reaction_rates_r(state + dt * k3$rate, x, time + dt, parameters)
  list(state = state + dt / 6 * (k1$rate + 2*k2$rate + 2*k3$rate + k4$rate),
       pi_delta = -dt / 6 * (k1$mean_j5 + 2*k2$mean_j5 +
                              2*k3$mean_j5 + k4$mean_j5))
}

prepare_diffusion <- function(system, dt) {
  identity <- Matrix::Diagonal(nrow(system$laplacian))
  unique_d <- unique(unname(system$model$diffusion))
  operators <- lapply(unique_d, function(d) {
    left <- identity - 0.5 * dt * d * system$laplacian
    right <- identity + 0.5 * dt * d * system$laplacian
    list(factor = Matrix::Cholesky(left, LDL = FALSE), right = right)
  })
  names(operators) <- format(unique_d, scientific = FALSE, trim = TRUE)
  operators
}

simulate_vcell <- function(system, end_time = system$model$time[["end"]],
                           dt = 0.02, save_times = c(0, 5, 10, 20, 50, 100),
                           progress = TRUE, start_time = 0,
                           initial_state = NULL, initial_PI = NULL,
                           reaction_substeps = 1L) {
  if (!exists("reaction_step_cpp", mode = "function") &&
      !exists(".vcell_kernel_attempted", envir = .GlobalEnv, inherits = FALSE))
    compile_vcell_kernel()
  reaction_step <- if (exists("reaction_step_cpp", mode = "function"))
    reaction_step_cpp else reaction_step_r
  stopifnot(dt > 0, end_time > start_time, start_time >= 0,
            length(reaction_substeps) == 1L, reaction_substeps >= 1,
            reaction_substeps == as.integer(reaction_substeps))
  reaction_substeps <- as.integer(reaction_substeps)
  interval <- end_time - start_time
  n_steps <- round(interval / dt)
  if (abs(n_steps * dt - interval) > 1e-10)
    stop("simulation interval must be divisible by dt")
  selected_save_times <- save_times[save_times >= start_time & save_times <= end_time]
  save_steps <- unique(round((selected_save_times - start_time) / dt))
  save_steps <- sort(unique(c(0L, save_steps, n_steps)))
  operators <- prepare_diffusion(system, dt)
  state <- if (is.null(initial_state)) system$initial_state else initial_state
  pi_value <- if (is.null(initial_PI))
    unname(system$model$constants[["PI_init_uM"]]) else initial_PI
  snapshots <- setNames(list(state), format(start_time, scientific = FALSE,
                                             trim = TRUE))
  pi_trace <- data.frame(time = start_time, PI = pi_value)
  diffusion <- unname(system$model$diffusion)
  groups <- split(seq_along(diffusion), format(diffusion, scientific = FALSE, trim = TRUE))

  for (step in seq_len(n_steps)) {
    time <- start_time + (step - 1) * dt
    reaction_x <- if (!is.null(system$stimulus_x)) system$stimulus_x else system$coords$x
    reaction_dt <- dt / (2 * reaction_substeps)
    for (substep in seq_len(reaction_substeps)) {
      first <- reaction_step(state, reaction_x,
        time + (substep - 1) * reaction_dt, reaction_dt,
        system$reaction_parameters)
      state <- first$state
      pi_value <- pi_value + first$pi_delta
    }
    for (key in names(groups)) {
      columns <- groups[[key]]
      rhs <- operators[[key]]$right %*% state[, columns, drop = FALSE]
      state[, columns] <- as.matrix(Matrix::solve(operators[[key]]$factor, rhs))
    }
    for (substep in seq_len(reaction_substeps)) {
      second <- reaction_step(state, reaction_x,
        time + dt / 2 + (substep - 1) * reaction_dt, reaction_dt,
        system$reaction_parameters)
      state <- second$state
      pi_value <- pi_value + second$pi_delta
    }
    if (step %in% save_steps) {
      saved_state <- state
      colnames(saved_state) <- names(system$model$initial)
      actual_time <- start_time + step * dt
      snapshots[[format(actual_time, scientific = FALSE, trim = TRUE)]] <- saved_state
      pi_trace <- rbind(pi_trace, data.frame(time = actual_time, PI = pi_value))
    }
    if (progress && step %% max(1L, floor(n_steps / 10L)) == 0L)
      message(sprintf("%3d%%  t = %.2f s", round(100 * step / n_steps),
                      start_time + step * dt))
  }
  list(system = system, dt = dt, start_time = start_time, end_time = end_time,
       snapshots = snapshots,
       PI = pi_trace, PI5K = system$model$constants[["PI5K_init_uM"]],
       PI3K = system$model$constants[["PI3K_init_uM"]],
       PTEN = system$model$constants[["PTEN_init_uM"]])
}

snapshot_long <- function(result, times = NULL, species = NULL) {
  available_times <- as.numeric(names(result$snapshots))
  if (is.null(times)) times <- available_times
  if (is.null(species)) species <- colnames(result$snapshots[[1]])
  rows <- lapply(times, function(time) {
    index <- which.min(abs(available_times - time))
    actual <- available_times[index]
    state <- result$snapshots[[index]]
    do.call(rbind, lapply(species, function(name) data.frame(
      x = result$system$coords$x, y = result$system$coords$y,
      time = actual, species = name, value = state[, name])))
  })
  do.call(rbind, rows)
}
