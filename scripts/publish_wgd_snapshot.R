# Convert completed PDE checkpoints into the small, immutable dataset consumed
# by the public report. This script never runs the PDE solver.
#
# Default input:
#   reports/output/wgd_spatial_scale/raw
# An alternate raw-results directory may be supplied as the first argument.

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Publishing requires the digest package.")
}
source("scenarios/wgd_spatial_scale.R")
source("R/wgd_spatial_analysis.R")

args <- commandArgs(trailingOnly = TRUE)
raw_dir <- if (length(args)) args[[1]] else "reports/output/wgd_spatial_scale/raw"
snapshot_dir <- "reports/assets/wgd_spatial/v1"
dir.create(snapshot_dir, recursive = TRUE, showWarnings = FALSE)

required_raw <- c("summary_metrics.csv", "dynamic_metrics.csv",
                  "conservation_summary.csv", "boundary_concentration_profiles.csv")
missing_raw <- required_raw[!file.exists(file.path(raw_dir, required_raw))]
if (length(missing_raw)) stop("Missing generated results: ", paste(missing_raw, collapse = ", "))

design <- wgd_spatial_design()
parameters <- wgd_spatial_parameters()
selected_times <- parameters$selected_spatial_times_s
selected_species <- c("cdc42_act", "Rac_act", "Rho_act")
spatial_rows <- list()
profile_rows <- list()

aggregate_path <- file.path(raw_dir, "aspect_ratio_results.rds")
worker_paths <- file.path(raw_dir, sprintf("aspect_ratio_result_%02d.rds", seq_len(nrow(design))))
if (!all(file.exists(worker_paths)) && !file.exists(aggregate_path)) {
  stop("Missing aspect-ratio checkpoints in ", raw_dir)
}
aggregate_results <- if (all(file.exists(worker_paths))) NULL else readRDS(aggregate_path)
checkpoint_paths <- if (all(file.exists(worker_paths))) worker_paths else aggregate_path

for (index in seq_len(nrow(design))) {
  label <- design$geometry[[index]]
  result <- if (is.null(aggregate_results)) readRDS(worker_paths[[index]]) else aggregate_results[[label]]
  system <- result$system
  origin_x <- system$model$geometry$origin[["x"]]
  relative_x <- (system$coords$x - origin_x) / system$model$geometry$extent[["x"]]
  bins <- pmin(50L, floor(relative_x * 50) + 1L)

  for (time in selected_times) {
    state <- result$snapshots[[as.character(time)]]
    if (is.null(state)) stop("Checkpoint ", label, " is missing t=", time)
    for (species in selected_species) {
      spatial_rows[[length(spatial_rows) + 1L]] <- data.frame(
        geometry = label, time = time, species = species,
        x = system$coords$x, y = system$coords$y,
        dx = system$dx, dy = system$dy, value = state[, species]
      )
      profile <- aggregate(state[, species], list(bin = bins), mean)
      profile_rows[[length(profile_rows) + 1L]] <- data.frame(
        geometry = label, time = time, species = species,
        relative_x = (profile$bin - 0.5) / 50, value = profile$x
      )
    }
  }
  rm(result)
  invisible(gc())
}

spatial <- do.call(rbind, spatial_rows)
longitudinal <- do.call(rbind, profile_rows)
boundary <- utils::read.csv(file.path(raw_dir, "boundary_concentration_profiles.csv"),
                            stringsAsFactors = FALSE)
# The perimeter trace repeats concentrations on adjacent exposed faces. Retain
# at most 201 evenly spaced observations per curve for report rendering.
boundary_key <- interaction(boundary$geometry, boundary$time, boundary$species, drop = TRUE)
boundary <- do.call(rbind, lapply(split(boundary, boundary_key), function(x) {
  x[unique(round(seq(1, nrow(x), length.out = min(201L, nrow(x))))), , drop = FALSE]
}))
row.names(boundary) <- NULL

utils::write.csv(utils::read.csv(file.path(raw_dir, "summary_metrics.csv")),
                 file.path(snapshot_dir, "summary_metrics.csv"), row.names = FALSE, na = "")
utils::write.csv(utils::read.csv(file.path(raw_dir, "dynamic_metrics.csv")),
                 file.path(snapshot_dir, "dynamic_metrics.csv"), row.names = FALSE, na = "")
utils::write.csv(utils::read.csv(file.path(raw_dir, "conservation_summary.csv")),
                 file.path(snapshot_dir, "conservation_summary.csv"), row.names = FALSE, na = "")
utils::write.csv(boundary, file.path(snapshot_dir, "boundary_profiles.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(longitudinal, file.path(snapshot_dir, "longitudinal_profiles.csv"),
                 row.names = FALSE, na = "")
saveRDS(spatial, file.path(snapshot_dir, "spatial_fields.rds"), compress = "xz")

source_files <- c("R/vcell_vcml.R", "R/vcell_spatial.R", "R/vcell_model.cpp",
                  "R/wgd_spatial_analysis.R", "scenarios/wgd_spatial_scale.R",
                  "scenarios/wgd_spatial_scale.vcml", "scripts/generate_wgd_spatial_results.R")
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
git_commit <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                       error = function(e) NA_character_)
git_status <- tryCatch(system2("git", c("status", "--porcelain", "--", source_files),
                               stdout = TRUE), error = function(e) "unknown")
published_files <- wgd_snapshot_required_files()
manifest <- list(
  schema_version = 1L,
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  source_git_commit = git_commit[[1]],
  source_worktree_dirty = length(git_status) > 0L,
  source_sha256 = setNames(lapply(source_files, sha256), source_files),
  source_vcml_sha256 = sha256(parameters$source_vcml),
  checkpoint_sha256 = setNames(lapply(checkpoint_paths, sha256), basename(checkpoint_paths)),
  parameters = parameters,
  design = design,
  r_version = R.version.string,
  file_sha256 = setNames(lapply(file.path(snapshot_dir, published_files), sha256),
                         published_files)
)
saveRDS(manifest, file.path(snapshot_dir, "manifest.rds"), compress = "xz")

writeLines(c(
  "# WGD spatial report snapshot v1", "",
  "This directory is a committed, plot-ready scientific result. Public report builds",
  "read these files and never invoke the PDE solver.", "",
  paste0("Generated (UTC): ", manifest$generated_utc),
  paste0("Source base commit: `", manifest$source_git_commit, "`"),
  paste0("Source worktree was dirty: `", manifest$source_worktree_dirty, "`"),
  paste0("VCML SHA-256: `", manifest$source_vcml_sha256, "`"), "",
  "Refresh explicitly with:", "",
  "```text", "Rscript scripts/generate_wgd_spatial_results.R",
  "Rscript scripts/publish_wgd_snapshot.R", "```"
), file.path(snapshot_dir, "README.md"))

snapshot <- read_wgd_snapshot(snapshot_dir)
validate_wgd_snapshot_values(snapshot)
message("Published and validated WGD report snapshot at ", snapshot_dir)
