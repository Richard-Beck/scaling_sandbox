source(file.path("R", "chiou_mcas.R"))

p <- chiou_fig6_parameters()
p$output_interval <- 0.25
p$initial_dt <- 1e-3
template <- read_chiou_initial_state()

# Warm the JIT before timing.
invisible(.chiou_integrate_fast(rep(1, 4), rep(1, 4), 1, 0.01, p))

cases <- expand.grid(n = c(40L, 80L), separate = c(FALSE, TRUE))
results <- lapply(seq_len(nrow(cases)), function(i) {
  n <- cases$n[i]
  p$n <- n
  u <- .chiou_resample_periodic(template$u, n)
  v <- .chiou_resample_periodic(template$v, n)
  # Move away from the archived steady profile so adaptive decisions are tested.
  u <- u * rep(c(1.03, 0.97), each = n / 2)
  v <- v * rep(c(1.03, 0.97), each = n / 2)
  mass_before <- chiou_mass(u, v, 10)

  fast_time <- system.time(fast <- .chiou_integrate_fast(
    u, v, 10, 1, p, separate = cases$separate[i], record_masses = TRUE))
  reference_time <- system.time(reference <- .chiou_integrate(
    u, v, 10, 1, p, separate = cases$separate[i], record = TRUE))

  h <- n / 2
  reference_peak_mass <- t(vapply(seq_len(nrow(reference$u_history)), function(j) {
    background <- min(reference$u_history[j, ])
    c(sum(pmax(reference$u_history[j, seq_len(h)] - background, 0)),
      sum(pmax(reference$u_history[j, h + seq_len(h)] - background, 0))) * 10 / n
  }, numeric(2)))

  data.frame(
    n = n, separate = cases$separate[i],
    max_abs_u = max(abs(fast$u - reference$u)),
    max_abs_v = max(abs(fast$v - reference$v)),
    max_abs_peak_mass = max(abs(cbind(fast$mass1, fast$mass2) - reference_peak_mass)),
    conservation_error = abs(chiou_mass(fast$u, fast$v, 10) - mass_before),
    fast_steps = fast$accepted_steps,
    reference_steps = reference$accepted_steps,
    fast_elapsed_s = unname(fast_time["elapsed"]),
    reference_elapsed_s = unname(reference_time["elapsed"])
  )
})
verification <- do.call(rbind, results)
verification$speedup <- verification$reference_elapsed_s / verification$fast_elapsed_s
print(verification)

stopifnot(max(verification$max_abs_u) < 1e-10)
stopifnot(max(verification$max_abs_v) < 1e-10)
stopifnot(max(verification$max_abs_peak_mass) < 1e-10)
stopifnot(max(verification$conservation_error) < 1e-10)
stopifnot(all(verification$fast_steps == verification$reference_steps))
message("Optimized Chiou solver agrees with the independent R reference.")
