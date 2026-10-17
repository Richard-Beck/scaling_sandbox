#!/usr/bin/env Rscript

# Convert an Eroume-Maree solver checkpoint into adapter-neutral raw metrics.
.libPaths(file.path(R.home(), "library"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: export_wgd_eroume_raw.R CHECKPOINT OUTPUT")
checkpoint <- args[[1]]
output <- args[[2]]
result <- readRDS(checkpoint)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

system <- result$system
origin <- system$model$geometry$origin[["x"]]
relative_x <- (system$coords$x - origin) / system$model$geometry$extent[["x"]]
species <- names(system$model$initial)
snapshots <- lapply(result$snapshots, function(state) {
  if (is.null(colnames(state))) colnames(state) <- species
  state
})
times <- as.numeric(names(snapshots))

count_peaks <- function(values) {
  profile <- aggregate(values, list(ix = system$coords$ix), mean)$x
  if (length(profile) < 3L || diff(range(profile)) < 1e-8) return(0L)
  smoothed <- as.numeric(stats::filter(profile, rep(0.2, 5), sides = 2))
  smoothed[is.na(smoothed)] <- profile[is.na(smoothed)]
  cutoff <- min(smoothed) + 0.2 * diff(range(smoothed))
  middle <- which(smoothed[2:(length(smoothed)-1)] >= smoothed[1:(length(smoothed)-2)] &
    smoothed[2:(length(smoothed)-1)] > smoothed[3:length(smoothed)]) + 1L
  candidates <- middle
  if (smoothed[1] > smoothed[2]) candidates <- c(1L, candidates)
  if (tail(smoothed, 1) > smoothed[length(smoothed)-1]) candidates <- c(candidates, length(smoothed))
  as.integer(sum(smoothed[candidates] >= cutoff))
}

rows <- lapply(seq_along(snapshots), function(index) {
  state <- snapshots[[index]]
  active <- state[, "cdc42_act"]
  front <- mean(active[relative_x >= 0.8])
  back <- mean(active[relative_x <= 0.2])
  data.frame(
    time = times[[index]],
    spatial_contrast = (max(active)-min(active))/(max(active)+min(active)+1e-12),
    spatial_cv = stats::sd(active)/mean(active),
    front_back_contrast = (front-back)/(front+back+1e-12),
    n_domains = count_peaks(active),
    native_integral = sum(active) * system$dx * system$dy
  )
})
trajectory_metrics <- do.call(rbind, rows)
final <- trajectory_metrics[nrow(trajectory_metrics), ]
target <- 0.9 * max(trajectory_metrics$spatial_contrast)
reached <- trajectory_metrics$time[trajectory_metrics$spatial_contrast >= target]

pool_names <- list(
  Cdc42 = c("cdc42_inact", "cdc42_act"),
  Rac = c("Rac_inact", "Rac_act"),
  Rho = c("Rho_inact", "Rho_act")
)
pools <- do.call(rbind, lapply(snapshots, function(state) {
  vapply(pool_names, function(columns) mean(rowSums(state[, columns, drop = FALSE])), numeric(1))
}))
reference <- pools[1, ]
conservation_error <- max(abs(sweep(pools, 2, reference, "-") /
  matrix(pmax(abs(reference), 1e-12), nrow(pools), length(reference), byrow = TRUE)))

state_class <- if (final$n_domains == 0 || final$spatial_contrast < 0.05) {
  "homogeneous"
} else if (final$n_domains > 1) "multi_domain" else "polarized"
raw <- data.frame(
  spatial_contrast = final$spatial_contrast,
  spatial_cv = final$spatial_cv,
  n_domains = final$n_domains,
  response_time = if (length(reached)) min(reached) else NA_real_,
  total_output = NA_real_, output_per_volume = NA_real_, efficiency = NA_real_,
  conservation_error = conservation_error, state_class = state_class,
  front_back_contrast = final$front_back_contrast,
  native_integral = final$native_integral
)
write.csv(raw, file.path(output, "raw_metrics.csv"), row.names = FALSE)
write.csv(trajectory_metrics, file.path(output, "trajectory.csv"), row.names = FALSE)
