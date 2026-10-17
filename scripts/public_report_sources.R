# Public catalogue report scope, including canonical locations after relocation.
# New workstream Rmd files enter the public site through this explicit list.
public_report_sources <- function(root = ".") {
  report_names <- c(
    "buttenschoen_cell_size_rac_polarity", "chiou_competition",
    "holmes_rac_cdc42_size_scaling", "membrane_to_nucleus_signaling",
    "mt_polymer_density_networks", "sedlack_single_cell_oxygen",
    "size_scaling_overview", "wave_pinning_polarity", "wgd_spatial_scaling"
  )
  canonical <- c(
    holmes_rac_cdc42_size_scaling =
      "literature_models/holmes_2012_rac_cdc42_rho/holmes_rac_cdc42_size_scaling.Rmd",
    wave_pinning_polarity =
      "literature_models/mori_2008_wave_pinning/wave_pinning_polarity.Rmd"
  )
  paths <- setNames(file.path(root, "reports", paste0(report_names, ".Rmd")),
                    report_names)
  for (name in names(canonical)) {
    if (!file.exists(paths[[name]])) paths[[name]] <- file.path(root, canonical[[name]])
  }
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Public report sources are missing: ", paste(missing, collapse = ", "))
  }
  paths
}
