wgd_snapshot_required_files <- function() {
  c(
    "summary_metrics.csv",
    "dynamic_metrics.csv",
    "conservation_summary.csv",
    "boundary_profiles.csv",
    "longitudinal_profiles.csv",
    "spatial_fields.rds"
  )
}

wgd_sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("SHA-256 validation requires the digest package.")
  }
  extension <- tolower(tools::file_ext(path))
  if (extension %in% c("csv", "md", "r", "rmd", "cpp", "vcml")) {
    bytes <- readBin(path, what = "raw", n = file.info(path)$size)
    # Git commonly stores text as LF but checks it out as CRLF on Windows.
    # Hash canonical LF bytes so a scientific snapshot validates identically
    # on Windows, macOS, and Linux without weakening binary-file validation.
    crlf <- bytes == as.raw(13L) & c(bytes[-1L] == as.raw(10L), FALSE)
    if (any(crlf)) bytes <- bytes[!crlf]
    return(digest::digest(bytes, algo = "sha256", serialize = FALSE))
  }
  digest::digest(path, algo = "sha256", file = TRUE)
}

read_wgd_snapshot <- function(path = "reports/assets/wgd_spatial/v1",
                              validate = TRUE) {
  manifest_path <- file.path(path, "manifest.rds")
  if (!file.exists(manifest_path)) {
    stop("Missing WGD report snapshot: ", manifest_path,
         ". Run scripts/publish_wgd_snapshot.R after generating the PDE results.")
  }
  manifest <- readRDS(manifest_path)
  required <- wgd_snapshot_required_files()
  missing <- required[!file.exists(file.path(path, required))]
  if (length(missing)) stop("Incomplete WGD report snapshot: ", paste(missing, collapse = ", "))

  if (validate) {
    if (!requireNamespace("digest", quietly = TRUE)) {
      stop("Snapshot validation requires the digest package.")
    }
    actual <- vapply(required, function(name) {
      wgd_sha256_file(file.path(path, name))
    }, character(1))
    expected <- unlist(manifest$file_sha256[required], use.names = TRUE)
    if (!identical(actual, expected)) stop("WGD report snapshot checksum validation failed.")
  }

  list(
    manifest = manifest,
    summary = utils::read.csv(file.path(path, "summary_metrics.csv"),
                              stringsAsFactors = FALSE),
    dynamics = utils::read.csv(file.path(path, "dynamic_metrics.csv"),
                               stringsAsFactors = FALSE),
    conservation = utils::read.csv(file.path(path, "conservation_summary.csv"),
                                   stringsAsFactors = FALSE),
    boundary = utils::read.csv(file.path(path, "boundary_profiles.csv"),
                               stringsAsFactors = FALSE),
    longitudinal = utils::read.csv(file.path(path, "longitudinal_profiles.csv"),
                                   stringsAsFactors = FALSE),
    spatial = readRDS(file.path(path, "spatial_fields.rds"))
  )
}

validate_wgd_snapshot_values <- function(snapshot) {
  stopifnot(
    identical(snapshot$manifest$schema_version, 1L),
    nrow(snapshot$summary) == 6L,
    all(c("0.5x", "1x", "2x") %in% snapshot$summary$size),
    all(snapshot$conservation$max_relative_drift < 1e-8),
    all(is.finite(snapshot$spatial$value))
  )
  invisible(snapshot)
}

wgd_geometry_colours <- function() {
  setNames(
    c("#56B4E9", "#0072B2", "#003B5C", "#F0A35E", "#D55E00", "#8C2D04"),
    c("AR 1.00 - 0.5x", "AR 1.00 - 1x", "AR 1.00 - 2x",
      "AR 1.29 - 0.5x", "AR 1.29 - 1x", "AR 1.29 - 2x")
  )
}
