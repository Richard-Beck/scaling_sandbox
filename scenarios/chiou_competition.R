source(file.path("R", "chiou_mcas.R"))

chiou_published_size_design <- function(
    lengths = c(5, 7.5, 10, 12.5, 15, 20), baseline_length = 10,
    baseline_mass = 20, parameters = chiou_fig6_parameters()) {
  sat <- chiou_saturation_values(parameters)
  rbind(
    data.frame(control = "constant overall concentration", length = lengths,
               total_mass = baseline_mass * lengths / baseline_length),
    data.frame(control = "constant total M", length = lengths, total_mass = baseline_mass),
    data.frame(control = "constant protein content in peaks", length = lengths,
               total_mass = baseline_mass + sat$q_sat * (lengths - baseline_length)))
}

run_chiou_published_size_experiment <- function(
    design = chiou_published_size_design(), parameters = chiou_fig6_parameters(), ...) {
  initial_cache <- new.env(parent = emptyenv())
  runs <- lapply(seq_len(nrow(design)), function(i) {
    key <- sprintf("n%d_L%.12g_M%.12g", parameters$n, design$length[i],
                   design$total_mass[i])
    if (!exists(key, envir = initial_cache, inherits = FALSE)) {
      dots <- list(...)
      initial <- do.call(chiou_competition_initial_state,
                         c(list(length = design$length[i],
                                total_mass = design$total_mass[i],
                                parameters = parameters),
                           dots[names(dots) %in% c("template", "equilibration_time", "engine")]))
      assign(key, initial, envir = initial_cache)
    }
    run_chiou_competition(design$length[i], design$total_mass[i], parameters,
                          initial_state = get(key, envir = initial_cache), ...)
  })
  summary <- cbind(design, do.call(rbind, lapply(runs, function(x) data.frame(
    competition_time = x$competition_time, censored_at = x$censored_at,
    u_max = x$u_max, saturation_index = x$saturation_index,
    peak_width = mean(x$peak_widths), peak_separation = x$peak_separation,
    final_number_of_peaks = x$final_number_of_peaks,
    mean_concentration = x$mean_concentration))))
  list(summary = summary, runs = runs)
}

chiou_wgd_design <- function(volume_ratios = c(1, 2),
                             abundance_exponents = c(0, 1/3, 1),
                             baseline_length = 10, baseline_mass = 20) {
  do.call(rbind, lapply(abundance_exponents, function(exponent) {
    out <- chiou_wgd_adapter(volume_ratios, baseline_length, baseline_mass, exponent)
    out$abundance_exponent <- exponent; out
  }))
}
