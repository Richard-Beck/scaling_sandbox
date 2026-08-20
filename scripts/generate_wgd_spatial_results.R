source("R/vcell_spatial.R")
source("scenarios/wgd_spatial_scale.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

out_dir <- "reports/output/wgd_spatial_scale/raw"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
design <- wgd_spatial_design()
parameters <- wgd_spatial_parameters()
volume_linear_ratio <- 2^(1 / 3)
worker_paths <- file.path(out_dir, sprintf("aspect_ratio_result_%02d.rds", seq_len(nrow(design))))

make_system <- function(length_um, width_um) {
  build_vcell_system(parameters$source_vcml, parameters$simulation,
    physical_extent = c(x = length_um, y = width_um), mesh_size = parameters$mesh)
}

args <- commandArgs(trailingOnly = TRUE)
worker_arg <- grep("^--worker=[1-6]$", args, value = TRUE)
if (length(worker_arg)) {
  worker_index <- as.integer(sub("^--worker=", "", worker_arg[[1]]))
  spec <- design[worker_index, ]
  compiled_kernel <- compile_vcell_kernel()
  if (!compiled_kernel || !exists("reaction_step_cpp", mode = "function")) {
    stop("Rcpp reaction kernel is required for this experiment but did not compile")
  }
  message(sprintf("Worker %d using Rcpp: %s (L=%.2f, W=%.2f um, mesh=80x80)",
                  worker_index, spec$geometry, spec$length_um, spec$width_um))
  system <- make_system(spec$length_um, spec$width_um)
  result <- simulate_vcell(system, end_time = parameters$end_time_s,
    dt = parameters$timestep_s,
    save_times = seq(parameters$start_time_s, parameters$end_time_s,
                     by = parameters$save_interval_s))
  saveRDS(result, worker_paths[[worker_index]], compress = FALSE)
  message("Worker ", worker_index, " completed: ", worker_paths[[worker_index]])
}

if (!length(worker_arg)) {
reuse_workers <- "--reuse" %in% args
workers_to_run <- if (reuse_workers) which(!file.exists(worker_paths)) else seq_len(nrow(design))
if (length(workers_to_run)) {
  compiled_kernel <- compile_vcell_kernel()
  if (!compiled_kernel || !exists("reaction_step_cpp", mode = "function")) {
    stop("Rcpp reaction kernel is required for this experiment but did not compile")
  }
  for (worker_index in workers_to_run) {
    spec <- design[worker_index, ]
    message(sprintf("Running %d/%d: %s", worker_index, nrow(design), spec$geometry))
    system <- make_system(spec$length_um, spec$width_um)
    result <- simulate_vcell(system, end_time = parameters$end_time_s,
      dt = parameters$timestep_s,
      save_times = seq(parameters$start_time_s, parameters$end_time_s,
                       by = parameters$save_interval_s))
    saveRDS(result, worker_paths[[worker_index]], compress = FALSE)
    rm(result, system)
    invisible(gc())
  }
}
if (!all(file.exists(worker_paths))) {
  stop("One or more worker checkpoints were not produced")
}
results <- setNames(lapply(worker_paths, readRDS), design$geometry)
results <- results[design$geometry]
dropped_snapshots <- character()
for (label in names(results)) {
  species_names <- names(results[[label]]$system$model$initial)
  results[[label]]$snapshots <- lapply(results[[label]]$snapshots, function(state) {
    if (is.null(colnames(state))) colnames(state) <- species_names
    state
  })
  finite_snapshot <- vapply(results[[label]]$snapshots,
                            function(state) all(is.finite(state)), logical(1))
  if (any(!finite_snapshot)) {
    dropped_snapshots <- c(dropped_snapshots,
      paste0(label, " at t=", names(results[[label]]$snapshots)[!finite_snapshot]))
    results[[label]]$snapshots <- results[[label]]$snapshots[finite_snapshot]
  }
}
analysis_endpoint <- min(vapply(results, function(result)
  max(as.numeric(names(result$snapshots))), numeric(1)))
if (length(dropped_snapshots)) {
  warning("Excluded non-finite snapshots: ", paste(dropped_snapshots, collapse = "; "))
}
theme_null <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")
geometry_colours <- setNames(c("#56B4E9", "#0072B2", "#003B5C",
                               "#F0A35E", "#D55E00", "#8C2D04"),
                             design$geometry)

relative_x <- function(system) {
  origin <- system$model$geometry$origin[["x"]]
  (system$coords$x - origin) / system$model$geometry$extent[["x"]]
}

front_back <- function(values, u) {
  front <- mean(values[u >= 0.8])
  back <- mean(values[u <= 0.2])
  c(front = front, back = back, ratio = front / back,
    contrast = (front - back) / (front + back))
}

weighted_position <- function(values, u) {
  weights <- pmax(values - min(values), 0)
  if (sum(weights) <= .Machine$double.eps) return(mean(u))
  sum(u * weights) / sum(weights)
}

boundary_mean <- function(values, system) {
  mask <- system$mask
  padded <- matrix(FALSE, nrow(mask) + 2L, ncol(mask) + 2L)
  padded[2:(nrow(mask) + 1L), 2:(ncol(mask) + 1L)] <- mask
  ix <- system$coords$ix + 1L
  iy <- system$coords$iy + 1L
  exposed_x <- (!padded[cbind(ix - 1L, iy)]) +
    (!padded[cbind(ix + 1L, iy)])
  exposed_y <- (!padded[cbind(ix, iy - 1L)]) +
    (!padded[cbind(ix, iy + 1L)])
  exposed_length <- exposed_x * system$dy + exposed_y * system$dx
  weighted.mean(values, exposed_length)
}

boundary_trace <- function(system) {
  mask <- system$mask
  nx <- nrow(mask); ny <- ncol(mask)
  inside <- function(i, j) i >= 1L && i <= nx && j >= 1L && j <= ny && mask[i, j]
  edges <- vector("list", 4L * sum(mask))
  edge_count <- 0L
  add_edge <- function(x0, y0, x1, y1, cell_index) {
    edge_count <<- edge_count + 1L
    edges[[edge_count]] <<- c(x0 = x0, y0 = y0, x1 = x1, y1 = y1,
                              cell_index = cell_index)
  }
  active <- which(mask, arr.ind = TRUE)
  for (row in seq_len(nrow(active))) {
    i <- active[row, 1]; j <- active[row, 2]
    cell <- system$cell_index[i, j]
    # Exposed faces directed counterclockwise around their owning cell.
    if (!inside(i, j - 1L)) add_edge(i - 1L, j - 1L, i, j - 1L, cell)
    if (!inside(i + 1L, j)) add_edge(i, j - 1L, i, j, cell)
    if (!inside(i, j + 1L)) add_edge(i, j, i - 1L, j, cell)
    if (!inside(i - 1L, j)) add_edge(i - 1L, j, i - 1L, j - 1L, cell)
  }
  edges <- as.data.frame(do.call(rbind, edges[seq_len(edge_count)]))
  start_key <- paste(edges$x0, edges$y0, sep = ",")
  end_key <- paste(edges$x1, edges$y1, sep = ",")
  unused <- rep(TRUE, nrow(edges))
  components <- list()
  while (any(unused)) {
    first <- which(unused)[1]
    component <- integer()
    current <- first
    repeat {
      component <- c(component, current)
      unused[current] <- FALSE
      next_edge <- which(unused & start_key == end_key[current])
      if (!length(next_edge)) break
      current <- next_edge[1]
    }
    components[[length(components) + 1L]] <- component
  }
  edge_length <- sqrt(((edges$x1 - edges$x0) * system$dx)^2 +
                        ((edges$y1 - edges$y0) * system$dy)^2)
  component_length <- vapply(components, function(index) sum(edge_length[index]), numeric(1))
  order <- components[[which.max(component_length)]]
  if (end_key[tail(order, 1)] != start_key[order[1]]) {
    stop("Exposed-face traversal did not produce a closed perimeter")
  }
  # Reverse the counterclockwise face order. In clockwise traversal, each
  # face begins at its original endpoint.
  order <- rev(order)
  traverse_x <- edges$x1[order]
  traverse_y <- edges$y1[order]
  rear_candidates <- which(traverse_x == min(traverse_x))
  start <- rear_candidates[which.min(abs(traverse_y[rear_candidates] - ny / 2))]
  order <- c(order[seq.int(start, length(order))],
             if (start > 1L) order[seq_len(start - 1L)])
  lengths <- edge_length[order]
  perimeter_um <- c(0, cumsum(lengths))
  cell_index <- c(as.integer(edges$cell_index[order]),
                  as.integer(edges$cell_index[order[1]]))
  data.frame(perimeter_um = perimeter_um, cell_index = cell_index)
}

count_peaks <- function(values, ix) {
  profile <- aggregate(values, list(ix = ix), mean)
  y <- as.numeric(stats::filter(profile$x, rep(1 / 5, 5), sides = 2))
  y[is.na(y)] <- profile$x[is.na(y)]
  amplitude <- diff(range(y))
  if (amplitude < 1e-8) return(c(peaks = 0, dominance = NA_real_))
  candidates <- which(y[2:(length(y) - 1)] >= y[1:(length(y) - 2)] &
                        y[2:(length(y) - 1)] > y[3:length(y)]) + 1L
  if (y[1] > y[2]) candidates <- c(1L, candidates)
  if (y[length(y)] > y[length(y) - 1]) candidates <- c(candidates, length(y))
  candidates <- candidates[y[candidates] >= min(y) + 0.20 * amplitude]
  if (!length(candidates)) return(c(peaks = 0, dominance = NA_real_))
  ordered <- candidates[order(y[candidates], decreasing = TRUE)]
  accepted <- ordered[1]
  if (length(ordered) > 1) {
    for (candidate in ordered[-1]) {
      nearest <- accepted[which.min(abs(accepted - candidate))]
      between <- seq.int(min(candidate, nearest), max(candidate, nearest))
      prominence <- (y[candidate] - min(y[between])) / amplitude
      if (prominence >= 0.10) accepted <- c(accepted, candidate)
    }
  }
  heights <- sort(y[accepted], decreasing = TRUE)
  dominance <- if (length(heights) >= 2) heights[1] / heights[2] else Inf
  c(peaks = length(accepted), dominance = dominance)
}

metric_rows <- list()
profile_rows <- list()
k <- 1L
for (geometry in names(results)) {
  result <- results[[geometry]]
  system <- result$system
  u <- relative_x(system)
  for (time_label in names(result$snapshots)) {
    time <- as.numeric(time_label)
    state <- result$snapshots[[time_label]]
    cdc <- state[, "cdc42_act"]
    rac <- state[, "Rac_act"]
    rho <- state[, "Rho_act"]
    cdc_fb <- front_back(cdc, u)
    rac_fb <- front_back(rac, u)
    rho_fb <- front_back(rho, u)
    cdc_peaks <- count_peaks(cdc, system$coords$ix)
    rac_peaks <- count_peaks(rac, system$coords$ix)
    metric_rows[[k]] <- data.frame(
      geometry = geometry, time = time,
      cdc42_ratio = cdc_fb[["ratio"]], cdc42_contrast = cdc_fb[["contrast"]],
      rac_ratio = rac_fb[["ratio"]], rac_contrast = rac_fb[["contrast"]],
      rho_front_back_ratio = rho_fb[["ratio"]],
      rho_separation = weighted_position(cdc, u) - weighted_position(rho, u),
      cdc42_boundary_mean = boundary_mean(cdc, system),
      rac_boundary_mean = boundary_mean(rac, system),
      rho_boundary_mean = boundary_mean(rho, system),
      cdc42_peaks = cdc_peaks[["peaks"]], rac_peaks = rac_peaks[["peaks"]],
      cdc42_dominance = cdc_peaks[["dominance"]],
      cdc42_peak_position = u[which.max(cdc)]
    )
    k <- k + 1L

    bins <- pmin(50L, floor(u * 50) + 1L)
    for (species in c("cdc42_act", "Rac_act", "Rho_act")) {
      prof <- aggregate(state[, species], list(bin = bins), mean)
      profile_rows[[length(profile_rows) + 1L]] <- data.frame(
        geometry = geometry, time = time, species = species,
        relative_x = (prof$bin - 0.5) / 50, value = prof$x)
    }
  }
}
metrics <- do.call(rbind, metric_rows)
profiles <- do.call(rbind, profile_rows)
write.csv(metrics, file.path(out_dir, "dynamic_metrics.csv"), row.names = FALSE)

establishment_time <- function(time, contrast) {
  target <- 0.9 * mean(contrast[time >= 180 & time <= 200])
  candidates <- time[contrast >= target & time <= 200]
  if (length(candidates)) min(candidates) else NA_real_
}

summary_rows <- lapply(split(metrics, metrics$geometry), function(x) {
  late <- x$time >= 180 & x$time <= 200
  final <- x[x$time == max(x$time), ]
  spec <- design[design$geometry == x$geometry[1], ]
  data.frame(
    geometry = x$geometry[1], pair = spec$pair, size = spec$size,
    length_um = spec$length_um, width_um = spec$width_um,
    aspect_ratio = spec$length_um / spec$width_um,
    volume_ratio = spec$volume_ratio,
    time_to_90pct_cdc42_polarity_s = establishment_time(x$time, x$cdc42_contrast),
    final_cdc42_front_back = final$cdc42_ratio,
    final_rac_front_back = final$rac_ratio,
    final_rho_front_back = final$rho_front_back_ratio,
    final_rho_separation_cell_lengths = final$rho_separation,
    final_cdc42_boundary_mean_uM = final$cdc42_boundary_mean,
    final_rac_boundary_mean_uM = final$rac_boundary_mean,
    final_rho_boundary_mean_uM = final$rho_boundary_mean,
    final_cdc42_peaks = final$cdc42_peaks,
    max_cdc42_peaks_after_stimulus = max(x$cdc42_peaks[x$time >= 10]),
    cdc42_multi_peak_duration_s = sum(x$cdc42_peaks[x$time >= 10 & x$time <= 200] > 1),
    rac_multi_peak_duration_s = sum(x$rac_peaks[x$time >= 10 & x$time <= 200] > 1),
    max_late_cdc42_peaks = max(x$cdc42_peaks[late]),
    late_peak_position_sd = sd(x$cdc42_peak_position[late]),
    late_contrast_cv = sd(x$cdc42_contrast[late]) / mean(x$cdc42_contrast[late])
  )
})
summary_table <- do.call(rbind, summary_rows)
summary_table$unique_state_stable <-
  summary_table$final_cdc42_peaks == 1 &
  summary_table$max_late_cdc42_peaks == 1 &
  summary_table$late_peak_position_sd < 0.02 &
  summary_table$late_contrast_cv < 0.02
write.csv(summary_table, file.path(out_dir, "summary_metrics.csv"), row.names = FALSE)

conservation_summary <- do.call(rbind, lapply(names(results), function(label) {
  snapshots <- results[[label]]$snapshots
  totals <- do.call(rbind, lapply(snapshots, function(state) c(
    Cdc42 = mean(state[, "cdc42_inact"] + state[, "cdc42_act"]),
    Rac = mean(state[, "Rac_inact"] + state[, "Rac_act"]),
    Rho = mean(state[, "Rho_inact"] + state[, "Rho_act"])
  )))
  reference <- totals[1, ]
  data.frame(geometry = label, species = colnames(totals),
    max_relative_drift = apply(abs(sweep(totals, 2, reference, "-")), 2, max) /
      pmax(abs(reference), 1e-12))
}))
write.csv(conservation_summary, file.path(out_dir, "conservation_summary.csv"),
          row.names = FALSE)

# Experimental design: the same mask and relative stimulus, but physical distances scale.
geometry_data <- do.call(rbind, lapply(names(results), function(label) {
  system <- results[[label]]$system
  data.frame(geometry = label, system$coords,
             relative_x = relative_x(system), dx = system$dx, dy = system$dy)
}))
p_geometry <- ggplot(geometry_data, aes(x, y, fill = relative_x)) +
  geom_tile(aes(width = dx, height = dy)) + coord_fixed() + facet_wrap(~geometry) +
  scale_fill_viridis_c(option = "plasma", name = "Relative\nposition") +
  labs(title = "Matched-volume aspect-ratio design",
       subtitle = "0.5x, 1x, and 2x preserve aspect ratio under implied 3D-volume scaling",
       x = "Physical x (um)", y = "Physical y (um)") + theme_null
ggsave(file.path(out_dir, "01_experimental_design.png"), p_geometry,
       width = 10, height = 4.4, dpi = 200, bg = "white")

contrast_long <- rbind(
  data.frame(metrics[c("geometry", "time")], species = "Cdc42", contrast = metrics$cdc42_contrast),
  data.frame(metrics[c("geometry", "time")], species = "Rac", contrast = metrics$rac_contrast)
)
p_contrast <- ggplot(subset(contrast_long, time <= 200),
                     aes(time, contrast, colour = geometry)) +
  geom_vline(xintercept = 10, linetype = 3, colour = "grey45") +
  geom_line(linewidth = 0.9) + facet_wrap(~species, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = geometry_colours) +
  labs(title = "Front-back polarization dynamics", subtitle = "Dotted line: spatial stimulus becomes uniform",
       x = "Time (s)", y = "(front - back) / (front + back)", colour = NULL) + theme_null

p_separation <- ggplot(subset(metrics, time <= 200),
                       aes(time, rho_separation, colour = geometry)) +
  geom_vline(xintercept = 10, linetype = 3, colour = "grey45") +
  geom_line(linewidth = 0.9) + scale_colour_manual(values = geometry_colours) +
  labs(title = "Cdc42-Rho separation", x = "Time (s)",
       y = "Centroid separation (cell lengths)", colour = NULL) + theme_null

peak_long <- rbind(
  data.frame(metrics[c("geometry", "time")], species = "Cdc42", peaks = metrics$cdc42_peaks),
  data.frame(metrics[c("geometry", "time")], species = "Rac", peaks = metrics$rac_peaks)
)
p_peaks <- ggplot(subset(peak_long, time <= 200),
                  aes(time, peaks, colour = geometry)) +
  geom_vline(xintercept = 10, linetype = 3, colour = "grey45") +
  geom_step(linewidth = 0.9) + facet_wrap(~species, ncol = 1) +
  scale_colour_manual(values = geometry_colours) +
  scale_y_continuous(breaks = 0:4) +
  labs(title = "Distinct longitudinal polarity peaks", x = "Time (s)",
       y = "Peak count", colour = NULL) + theme_null

dynamics_grob <- arrangeGrob(p_contrast, arrangeGrob(p_separation, p_peaks, ncol = 1),
                             ncol = 2, widths = c(1.15, 1),
                             top = "Geometry alone changes polarization dynamics")
ggsave(file.path(out_dir, "02_dynamics_dashboard.png"), dynamics_grob,
       width = 13, height = 8.5, dpi = 200, bg = "white")

profiles$species <- factor(profiles$species,
                           levels = c("cdc42_act", "Rac_act", "Rho_act"),
                           labels = c("Cdc42 active", "Rac active", "Rho active"))
kymograph_plots <- lapply(levels(profiles$species), function(field) {
  ggplot(subset(profiles, species == field & time <= 200),
         aes(relative_x, time, fill = value)) +
    geom_raster() + facet_wrap(~geometry, nrow = 2) +
    scale_fill_viridis_c(option = "magma") +
    labs(title = field, x = "Relative long-axis position (back to front)",
         y = "Time (s)", fill = "uM") +
    theme_null + theme(legend.position = "right")
})
kymograph_grob <- arrangeGrob(grobs = kymograph_plots, ncol = 1,
  top = "Polarization kymographs: species-specific scales shared across geometries")
ggsave(file.path(out_dir, "03_polarization_kymographs.png"), kymograph_grob,
       width = 11, height = 12, dpi = 200, bg = "white")

snapshot_times <- c(10, 20, 50, 100, 150, 200)
snapshot_data <- do.call(rbind, lapply(names(results), function(label) {
  result <- results[[label]]
  data <- snapshot_long(result, snapshot_times,
                        c("cdc42_act", "Rac_act", "Rho_act"))
  data$geometry <- label
  data$dx <- result$system$dx
  data$dy <- result$system$dy
  data
}))
snapshot_data$species <- factor(snapshot_data$species,
  levels = c("cdc42_act", "Rac_act", "Rho_act"),
  labels = c("Cdc42 active", "Rac active", "Rho active"))
snapshot_plots <- lapply(levels(snapshot_data$species), function(field) {
  ggplot(subset(snapshot_data, species == field), aes(x, y, fill = value)) +
    geom_tile(aes(width = dx, height = dy)) + coord_fixed() + facet_grid(geometry ~ time) +
    scale_fill_viridis_c(option = "magma") +
    labs(title = field, x = NULL, y = NULL, fill = "uM") +
    theme_null + theme(legend.position = "right")
})
snapshots_grob <- arrangeGrob(grobs = snapshot_plots, ncol = 1,
                              top = "Spatial fields across matched aspect-ratio pairs")
ggsave(file.path(out_dir, "04_spatial_snapshots.png"), snapshots_grob,
       width = 18, height = 15, dpi = 200, bg = "white")

late_profiles <- subset(profiles, time %in% c(10, 20, 50, 100, 150, 200))
p_profiles <- ggplot(late_profiles, aes(relative_x, value, colour = geometry)) +
  geom_line(linewidth = 0.75) + facet_grid(species ~ time, scales = "free_y") +
  scale_colour_manual(values = geometry_colours) +
  labs(title = "Matched-position longitudinal profiles", x = "Relative long-axis position",
       y = "Mean concentration (uM)", colour = NULL) + theme_null
ggsave(file.path(out_dir, "05_longitudinal_profiles.png"), p_profiles,
       width = 18, height = 8.5, dpi = 200, bg = "white")

boundary_times <- c(10, 20, 50, 100, 150, 200)
boundary_rows <- list()
for (geometry in names(results)) {
  result <- results[[geometry]]
  trace <- boundary_trace(result$system)
  spec <- design[design$geometry == geometry, ]
  for (time in boundary_times) {
    state <- result$snapshots[[as.character(time)]]
    for (species in c("cdc42_act", "Rac_act", "Rho_act")) {
      boundary_rows[[length(boundary_rows) + 1L]] <- data.frame(
        geometry = geometry, pair = spec$pair, size = spec$size,
        time = time, species = species,
        perimeter_um = trace$perimeter_um,
        relative_perimeter = trace$perimeter_um / max(trace$perimeter_um),
        concentration_uM = state[trace$cell_index, species])
    }
  }
}
boundary_profiles <- do.call(rbind, boundary_rows)
write.csv(boundary_profiles, file.path(out_dir, "boundary_concentration_profiles.csv"),
          row.names = FALSE)
boundary_profiles$species <- factor(boundary_profiles$species,
  levels = c("cdc42_act", "Rac_act", "Rho_act"),
  labels = c("Cdc42 active", "Rac active", "Rho active"))
boundary_profiles$size <- factor(boundary_profiles$size,
  levels = c("0.5x", "1x", "2x"))
size_colours <- c(`0.5x` = "#2166AC", `1x` = "#7B3294", `2x` = "#D73027")
plot_boundary_profiles <- function(pair_label, title_label, filename) {
  plot_data <- subset(boundary_profiles, pair == pair_label)
  plot <- ggplot(plot_data,
    aes(relative_perimeter, concentration_uM, colour = size)) +
    geom_line(linewidth = 0.8) +
    facet_grid(species ~ time, scales = "free_y") +
    scale_colour_manual(values = size_colours, drop = FALSE) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
    labs(title = paste("Local concentration along the cell perimeter:", title_label),
         subtitle = "Relative distance starts at the rear tip and proceeds clockwise",
         x = "Relative distance along perimeter", y = "Local concentration (uM)",
         colour = "Relative volume") + theme_null
  ggsave(file.path(out_dir, filename), plot,
         width = 18, height = 9, dpi = 200, bg = "white")
}
plot_boundary_profiles("L:W = 90:90", "aspect ratio 1.00",
                       "07a_boundary_profiles_ar_1.00.png")
plot_boundary_profiles("L:W = 90:70", "aspect ratio 1.29",
                       "07b_boundary_profiles_ar_1.29.png")

size_response <- rbind(
  data.frame(summary_table[c("geometry", "pair", "size", "volume_ratio")],
    readout = "Polarity establishment (s)",
    value = summary_table$time_to_90pct_cdc42_polarity_s),
  data.frame(summary_table[c("geometry", "pair", "size", "volume_ratio")],
    readout = "Final Cdc42 front:back",
    value = summary_table$final_cdc42_front_back),
  data.frame(summary_table[c("geometry", "pair", "size", "volume_ratio")],
    readout = "Final Rac front:back",
    value = summary_table$final_rac_front_back),
  data.frame(summary_table[c("geometry", "pair", "size", "volume_ratio")],
    readout = "Cdc42-Rho separation (cell lengths)",
    value = summary_table$final_rho_separation_cell_lengths),
  data.frame(summary_table[c("geometry", "pair", "size", "volume_ratio")],
    readout = "Post-stimulus maximum peak count",
    value = summary_table$max_cdc42_peaks_after_stimulus),
  data.frame(summary_table[c("geometry", "pair", "size", "volume_ratio")],
    readout = "Cdc42 multi-peak duration (s)",
    value = summary_table$cdc42_multi_peak_duration_s)
)
p_size_response <- ggplot(size_response,
  aes(volume_ratio, value, group = pair, linetype = pair)) +
  geom_line(colour = "grey35", linewidth = 0.7) +
  geom_point(aes(colour = geometry), size = 3) +
  facet_wrap(~readout, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = c(0.5, 1, 2), labels = c("0.5x", "1x", "2x"),
                     expand = expansion(mult = 0.12)) +
  scale_colour_manual(values = geometry_colours) +
  labs(title = "Matched 0.5x/1x/2x comparisons by aspect ratio",
       subtitle = paste0("Final readouts at ", analysis_endpoint, " s"),
       x = "Relative volume", y = NULL, colour = NULL, linetype = "Aspect ratio") + theme_null
ggsave(file.path(out_dir, "06_size_response_summary.png"), p_size_response,
       width = 12, height = 7.5, dpi = 200, bg = "white")

writeLines(c(
  "# WGD matched aspect-ratio simulations",
  "",
  "Six geometries compare matched 0.5x/1x/2x volumes at two teardrop aspect ratios.",
  sprintf("Square series: 2x L=90, W=90 um; 1x L=%.2f, W=%.2f um; 0.5x L=%.2f, W=%.2f um.",
          90 / volume_linear_ratio, 90 / volume_linear_ratio,
          90 / 4^(1 / 3), 90 / 4^(1 / 3)),
  sprintf("Elongated series: 2x L=90, W=70 um; 1x L=%.2f, W=%.2f um; 0.5x L=%.2f, W=%.2f um.",
          90 / volume_linear_ratio, 70 / volume_linear_ratio,
          90 / 4^(1 / 3), 70 / 4^(1 / 3)),
  "All dimensions preserve aspect ratio and scale with the cube root of relative implied 3D volume.",
  "The original teardrop mask is stretched independently to each specified physical length and width on an 80 x 80 mesh.",
  "All initial concentrations, kinetic constants, and physical diffusion coefficients are unchanged.",
  "The reaction system is evaluated with the compiled Rcpp kernel; the workflow stops rather than falling back if compilation fails.",
  "The stimulus is evaluated at baseline-equivalent x, preserving its amplitude and relative position.",
  "The implied third dimension scales by the same linear factor, giving the stated 0.5x, 1x, and 2x volumes and molecule counts at fixed concentration.",
  "Simulations run from 0-200 s at dt=0.02 s, with dynamics saved every second and spatial outputs emphasized at 100, 150, and 200 s.",
  "",
  "Time-to-polarity is the first saved second reaching 90% of the mean Cdc42 contrast at 180-200 s.",
  "Front and back are the terminal 20% of relative cell length.",
  "Rho separation is the active-Cdc42 minus active-Rho longitudinal centroid in cell lengths.",
  "Boundary profiles report local active-species concentration versus relative distance around the perimeter, starting at the rear tip and proceeding clockwise; separate figures show each aspect ratio.",
  "Peaks are local maxima in a five-bin-smoothed longitudinal profile above 20% of its range and with at least 10% prominence.",
  "Late stability is summarized over 180-200 s.",
  "A stable unique state requires one late Cdc42 peak, peak-position SD <0.02 cell lengths, and contrast CV <0.02.",
  "",
  "## Result synopsis",
  paste("Polarity-establishment times (s):",
        paste(summary_table$geometry,
              summary_table$time_to_90pct_cdc42_polarity_s, collapse = "; ")),
  paste("Post-stimulus Cdc42 multi-peak durations (s):",
        paste(summary_table$geometry,
              summary_table$cdc42_multi_peak_duration_s, collapse = "; ")),
  "See the figures and summary_metrics.csv for the paired size and aspect-ratio comparisons."
), file.path(out_dir, "README.md"))

message("Wrote WGD matched aspect-ratio experiment to ", out_dir)
print(summary_table)
}
