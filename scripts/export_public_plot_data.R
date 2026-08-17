# Export the exact data frames used by the R Markdown plots to docs/data.
# Run from the project root:
#   Rscript scripts/export_public_plot_data.R

write_public_csv <- function(data, relative_path) {
  target <- file.path("docs", "data", relative_path)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data, target, row.names = FALSE, na = "")
  message("Wrote ", target)
}

run_report <- function(report_path) {
  exported_script <- tempfile(fileext = ".R")
  knitr::purl(report_path, output = exported_script, quiet = TRUE)
  report_environment <- new.env(parent = globalenv())
  original_directory <- getwd()
  on.exit(setwd(original_directory), add = TRUE)
  setwd(dirname(report_path))
  source(exported_script, local = report_environment, echo = FALSE)
  report_environment
}

size_scaling <- run_report("reports/size_scaling_overview.Rmd")
write_public_csv(size_scaling$external_profiles,
                 "size_scaling_overview/external_profiles.csv")
write_public_csv(size_scaling$external_boundaries,
                 "size_scaling_overview/external_boundaries.csv")
write_public_csv(size_scaling$external_consumption_summary,
                 "size_scaling_overview/external_consumption_summary.csv")
write_public_csv(size_scaling$combined_profiles,
                 "size_scaling_overview/interior_exterior_profiles.csv")
write_public_csv(size_scaling$medium_boundaries,
                 "size_scaling_overview/interior_exterior_boundaries.csv")
write_public_csv(size_scaling$medium_consumption_summary,
                 "size_scaling_overview/interior_exterior_consumption_summary.csv")
write_public_csv(size_scaling$cytoplasmic_profiles,
                 "size_scaling_overview/nuclear_export_concentration_profiles.csv")
write_public_csv(size_scaling$cytoplasmic_sink_profiles,
                 "size_scaling_overview/nuclear_export_sink_profiles.csv")
write_public_csv(size_scaling$cytoplasmic_boundaries,
                 "size_scaling_overview/nuclear_export_boundaries.csv")
write_public_csv(size_scaling$cytoplasmic_consumption_summary,
                 "size_scaling_overview/nuclear_export_consumption_summary.csv")

oxygen <- run_report("reports/sedlack_single_cell_oxygen.Rmd")
write_public_csv(oxygen$kinetics,
                 "sedlack_single_cell_oxygen/oxygen_kinetics.csv")
write_public_csv(oxygen$calibration_summary,
                 "sedlack_single_cell_oxygen/calibration_summary.csv")
write_public_csv(oxygen$validation,
                 "sedlack_single_cell_oxygen/validation.csv")
write_public_csv(oxygen$slice_data,
                 "sedlack_single_cell_oxygen/central_slice_profiles.csv")
write_public_csv(oxygen$size_summary,
                 "sedlack_single_cell_oxygen/size_scaling_summary.csv")

membrane_to_nucleus <- run_report("reports/membrane_to_nucleus_signaling.Rmd")
write_public_csv(membrane_to_nucleus$signal_profiles,
                 "membrane_to_nucleus_signaling/signal_profiles.csv")
write_public_csv(membrane_to_nucleus$profile_boundaries,
                 "membrane_to_nucleus_signaling/profile_boundaries.csv")
write_public_csv(membrane_to_nucleus$signal_fate,
                 "membrane_to_nucleus_signaling/signal_fate.csv")
write_public_csv(membrane_to_nucleus$delivery_scaling,
                 "membrane_to_nucleus_signaling/nuclear_delivery_scaling.csv")
write_public_csv(membrane_to_nucleus$capture_sensitivity,
                 "membrane_to_nucleus_signaling/capture_sensitivity.csv")
write_public_csv(membrane_to_nucleus$regime_heatmap,
                 "membrane_to_nucleus_signaling/parameter_regime_heatmap.csv")
write_public_csv(membrane_to_nucleus$pde_zero_control,
                 "membrane_to_nucleus_signaling/pde_zero_control.csv")
write_public_csv(membrane_to_nucleus$source_scaling_control,
                 "membrane_to_nucleus_signaling/source_scaling_control.csv")
