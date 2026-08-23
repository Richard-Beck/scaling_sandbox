# Chiou et al. mass-conserved activator-substrate (MCAS) competition model.
# Defaults transcribe Main1D_color_plot.m; the alternate set is Parameters.m.

chiou_fig6_parameters <- function() {
  list(D_u = 0.01, D_v = 1, a = 1, b = 1, k = 0,
       baseline_length = 10, baseline_mass = 20, n = 500L,
       reaction_error = 0.05, steady_error = 2e-4,
       initial_dt = 1e-4, dt_growth = 1.5,
       equilibration_time = 2000, competition_limit = 1e5,
       output_interval = 1, winner_fraction = 0.99)
}

chiou_general_parameters <- function() {
  p <- chiou_fig6_parameters(); p$k <- 0.01; p
}

chiou_reaction <- function(u, v, parameters = chiou_fig6_parameters()) {
  with(parameters, a * u^2 * v / (1 + k * u^2) - b * u)
}

chiou_saturation_values <- function(parameters = chiou_fig6_parameters()) {
  with(parameters, list(u_sat = sqrt(2 * b * D_v / (a * D_u)),
                        q_sat = sqrt(9 * b * D_u / (2 * a * D_v))))
}

# Read the preserved MATLAB state. Use R.matlab if available, otherwise scipy.
read_chiou_initial_state <- function(
    path = file.path("dev", "chiou", "MATLAB_code", "initial_state_1D.mat")) {
  if (!file.exists(path)) stop("Chiou initial state not found: ", path)
  if (requireNamespace("R.matlab", quietly = TRUE)) {
    state <- R.matlab::readMat(path)
    return(list(u = as.numeric(state$u), v = as.numeric(state$v)))
  }
  python <- Sys.which(c("python", "python3")); python <- unname(python[nzchar(python)][1])
  if (is.na(python) || !nzchar(python))
    stop("Reading the original MAT file requires R.matlab or Python with scipy.")
  bridge <- file.path("scripts", "read_chiou_mat.py")
  if (!file.exists(bridge)) stop("MAT import bridge not found: ", bridge)
  raw <- system2(python, c(shQuote(normalizePath(bridge)), shQuote(normalizePath(path))),
                 stdout = TRUE, stderr = TRUE)
  status <- attr(raw, "status")
  if (!is.null(status) && status != 0L) stop(paste(raw, collapse = "\n"))
  raw <- raw[nzchar(raw)]
  values <- suppressWarnings(as.numeric(raw))
  values <- values[is.finite(values)]
  if (!length(values) || length(values) %% 2L)
    stop("Unexpected output while reading ", path)
  h <- length(values) %/% 2L
  list(u = values[seq_len(h)], v = values[h + seq_len(h)])
}

chiou_mass <- function(u, v, length) mean(u + v) * length

.chiou_resample_periodic <- function(y, n) {
  old_n <- length(y); x_old <- (seq_len(old_n) - 1) / old_n
  approx(c(x_old, 1), c(y, y[1]), xout = (seq_len(n) - 1) / n, rule = 2)$y
}

.chiou_diffuse_implicit <- function(y, D, length, dt) {
  n <- length(y); dx <- length / n; mode <- 0:(n - 1)
  eigenvalue <- -4 * sin(pi * mode / n)^2 / dx^2
  Re(fft(fft(y) / (1 - dt * D * eigenvalue), inverse = TRUE) / n)
}

.chiou_error_max <- function(x) {
  x <- x[is.finite(x)]; if (length(x)) max(x) else Inf
}

# Adaptive explicit-reaction/implicit-diffusion integration used by MATLAB.
# With separate=TRUE, periodic diffusion is applied independently to each half.
.chiou_integrate <- function(u, v, length, end_time, parameters,
                             separate = FALSE, record = FALSE,
                             stop_fraction = NULL) {
  n <- length(u)
  if (length(v) != n || n %% 2L) stop("u and v need equal, even lengths.")
  t <- 0; dt <- parameters$initial_dt; next_record <- 0; steady <- FALSE
  times <- numeric(); us <- list(); vs <- list(); accepted <- 0L
  save_state <- function() {
    times <<- c(times, t); us[[length(us) + 1L]] <<- u; vs[[length(vs) + 1L]] <<- v
  }
  if (record) save_state()
  diffuse <- function(y, D, step) {
    if (!separate) return(.chiou_diffuse_implicit(y, D, length, step))
    h <- n %/% 2L
    c(.chiou_diffuse_implicit(y[seq_len(h)], D, length / 2, step),
      .chiou_diffuse_implicit(y[h + seq_len(h)], D, length / 2, step))
  }
  while (t < end_time - 1e-12) {
    # MATLAB forces a step onto every integer-second boundary even when that
    # state is not retained; keep the same checkpoint behavior here.
    target <- min(end_time, next_record + parameters$output_interval)
    remaining <- target - t
    boundary_step <- dt >= remaining
    step <- min(dt, remaining)
    if (step <= 1e-14) {
      next_record <- target; if (record) save_state(); next
    }
    f1 <- chiou_reaction(u, v, parameters)
    fhalf <- chiou_reaction(u + step * f1 / 2, v - step * f1 / 2, parameters)
    err <- .chiou_error_max(abs((fhalf - f1) / fhalf))
    if (err < parameters$reaction_error || steady || boundary_step || step < 1e-13) {
      un <- diffuse(u + step * f1, parameters$D_u, step)
      vn <- diffuse(v - step * f1, parameters$D_v, step)
      if (any(!is.finite(c(un, vn))) || any(c(un, vn) < -1e-7))
        stop("Non-finite or materially negative concentration during integration.")
      ss_error <- .chiou_error_max(abs((un - u) / un) + abs((vn - v) / vn))
      steady <- ss_error < parameters$steady_error
      u <- pmax(un, 0); v <- pmax(vn, 0); t <- t + step
      accepted <- accepted + 1L; dt <- step * parameters$dt_growth
      if (abs(t - target) < 1e-9) {
        next_record <- target
        if (record) save_state()
        if (record && !is.null(stop_fraction)) {
          h <- n %/% 2L
          background <- min(u)
          masses <- c(sum(pmax(u[seq_len(h)] - background, 0)),
                      sum(pmax(u[h + seq_len(h)] - background, 0))) * length / n
          if (max(masses) / sum(masses) >= stop_fraction) break
        }
      }
    } else dt <- step / parameters$dt_growth
  }
  list(u = u, v = v, time = t, times = times,
       u_history = if (record) do.call(rbind, us) else NULL,
       v_history = if (record) do.call(rbind, vs) else NULL,
       accepted_steps = accepted)
}

load_chiou_fast_solver <- local({
  loaded <- FALSE
  function() {
    if (!requireNamespace("reticulate", quietly = TRUE))
      stop("The optimized Chiou solver requires reticulate and Python with numba.")
    if (!loaded || !exists("chiou_integrate_numba", mode = "function")) {
      python <- unname(Sys.which(c("python", "python3")))
      python <- python[nzchar(python)][1]
      if (is.na(python) || !nzchar(python)) stop("Python was not found.")
      reticulate::use_python(python, required = TRUE)
      environment <- new.env(parent = globalenv())
      reticulate::source_python(file.path("scripts", "chiou_solver_fast.py"),
                                envir = environment)
      assign("chiou_integrate_numba", environment$integrate, envir = globalenv())
      loaded <<- TRUE
    }
    invisible(TRUE)
  }
})

.chiou_integrate_fast <- function(u, v, length, end_time, parameters,
                                  separate = FALSE, record_profiles = FALSE,
                                  record_masses = FALSE, stop_fraction = NULL) {
  if (record_profiles)
    stop("The optimized metric solver stores initial/final profiles only.")
  load_chiou_fast_solver()
  chiou_integrate_numba(u, v, length, end_time, parameters, separate,
                        record_masses, stop_fraction)
}

chiou_peak_metrics <- function(u, length, saturation = NULL) {
  n <- length(u); dx <- length / n
  peaks <- which(u >= c(u[n], u[-n]) & u > c(u[-1], u[1]))
  amplitude <- max(u) - min(u)
  peaks <- peaks[u[peaks] >= min(u) + 0.1 * amplitude]
  if (!length(peaks)) peaks <- which.max(u)
  widths <- vapply(peaks, function(i) {
    level <- min(u) + (u[i] - min(u)) / 2; count <- 1L
    j <- if (i == 1L) n else i - 1L
    while (count < n && u[j] >= level) { count <- count + 1L; j <- if (j == 1L) n else j - 1L }
    j <- if (i == n) 1L else i + 1L
    while (count < n && u[j] >= level) { count <- count + 1L; j <- if (j == n) 1L else j + 1L }
    count * dx
  }, numeric(1))
  separation <- NA_real_
  if (length(peaks) >= 2L) {
    top <- peaks[order(u[peaks], decreasing = TRUE)[1:2]]
    raw <- abs(diff(top)) * dx; separation <- min(raw, length - raw)
  }
  if (is.null(saturation)) saturation <- chiou_saturation_values()$u_sat
  list(number = length(peaks), positions = (peaks - 0.5) * dx,
       heights = u[peaks], widths = widths, separation = separation,
       u_max = max(u), saturation_index = max(u) / saturation)
}

# Faithful unequal competitors: give isolated half-domains 60% and 40% of M,
# equilibrate them, and only then permit diffusion across the interface.
chiou_competition_initial_state <- function(length = 10, total_mass = 20,
                                            parameters = chiou_fig6_parameters(),
                                            template = read_chiou_initial_state(),
                                            equilibration_time = parameters$equilibration_time,
                                            engine = c("optimized", "reference")) {
  engine <- match.arg(engine)
  n <- as.integer(parameters$n)
  if (n %% 2L) stop("parameters$n must be even.")
  h <- n %/% 2L; old_h <- length(template$u) %/% 2L
  tu <- .chiou_resample_periodic(template$u[seq_len(old_h)], h)
  tv <- .chiou_resample_periodic(template$v[seq_len(old_h)], h)
  current <- chiou_mass(tu, tv, length / 2); target <- total_mass * c(0.6, 0.4)
  u <- c(tu * target[1] / current, tu * target[2] / current)
  v <- c(tv * target[1] / current, tv * target[2] / current)
  eq <- if (engine == "optimized") {
    .chiou_integrate_fast(u, v, length, equilibration_time, parameters,
                          separate = TRUE)
  } else {
    .chiou_integrate(u, v, length, equilibration_time, parameters,
                     separate = TRUE, record = FALSE)
  }
  list(u = eq$u, v = eq$v, length = length, total_mass = total_mass,
       target_half_fractions = c(0.6, 0.4), accepted_steps = eq$accepted_steps)
}

run_chiou_competition <- function(length = 10, total_mass = 20,
                                  parameters = chiou_fig6_parameters(),
                                  template = read_chiou_initial_state(),
                                  equilibration_time = parameters$equilibration_time,
                                  competition_limit = parameters$competition_limit,
                                  engine = c("optimized", "reference"),
                                  initial_state = NULL,
                                  record_profiles = FALSE) {
  engine <- match.arg(engine)
  initial <- initial_state
  if (is.null(initial))
    initial <- chiou_competition_initial_state(length, total_mass, parameters,
                                                template, equilibration_time, engine)
  if (engine == "optimized") {
    run <- .chiou_integrate_fast(initial$u, initial$v, length, competition_limit,
                                 parameters, record_profiles = record_profiles,
                                 record_masses = TRUE,
                                 stop_fraction = parameters$winner_fraction)
    mass_history <- cbind(run$mass1, run$mass2)
  } else {
    run <- .chiou_integrate(initial$u, initial$v, length, competition_limit,
                            parameters, separate = FALSE, record = TRUE,
                            stop_fraction = parameters$winner_fraction)
    h <- length(run$u) %/% 2L
    mass_history <- t(vapply(seq_len(nrow(run$u_history)), function(i) {
      background <- min(run$u_history[i, ])
      c(sum(pmax(run$u_history[i, seq_len(h)] - background, 0)),
        sum(pmax(run$u_history[i, h + seq_len(h)] - background, 0))) * length / (2 * h)
    }, numeric(2)))
  }
  fractions <- mass_history / rowSums(mass_history)
  reached <- max(fractions[nrow(fractions), ]) >= parameters$winner_fraction
  sat <- chiou_saturation_values(parameters)$u_sat
  pre <- chiou_peak_metrics(initial$u, length, sat)
  final <- chiou_peak_metrics(run$u, length, sat)
  list(competition_time = if (reached) tail(run$times, 1) else NA_real_,
       censored_at = if (reached) NA_real_ else tail(run$times, 1),
       time = run$times, winner_mass = apply(mass_history, 1, max),
       loser_mass = apply(mass_history, 1, min), half_fractions = fractions,
       initial_u = initial$u, initial_v = initial$v, final_u = run$u, final_v = run$v,
       peak_heights = pre$heights, peak_widths = pre$widths,
       peak_separation = pre$separation, u_max = pre$u_max,
       saturation_index = pre$saturation_index, total_mass = total_mass,
       mean_concentration = total_mass / length,
       final_number_of_peaks = final$number,
       accepted_steps = run$accepted_steps,
       rejected_steps = if (!is.null(run$rejected_steps)) run$rejected_steps else NA_integer_,
       u_history = if (record_profiles) run$u_history else NULL,
       v_history = if (record_profiles) run$v_history else NULL,
       parameters = parameters)
}

chiou_wgd_adapter <- function(volume_ratio, baseline_length = 10,
                              baseline_mass = 20, abundance_exponent = 1) {
  stopifnot(all(volume_ratio > 0), length(abundance_exponent) == 1L)
  data.frame(volume_ratio = volume_ratio,
             length_ratio = volume_ratio^(1 / 3),
             length = baseline_length * volume_ratio^(1 / 3),
             total_mass = baseline_mass * volume_ratio^abundance_exponent,
             mean_concentration = baseline_mass * volume_ratio^abundance_exponent /
               (baseline_length * volume_ratio^(1 / 3)))
}
