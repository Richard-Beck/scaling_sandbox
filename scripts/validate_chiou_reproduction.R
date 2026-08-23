# Validate qualitative behaviors of the Chiou port. Use --full for the native
# 500-bin resolution and longer observation window; the default is a fast CI
# diagnostic that is not a source of publication-grade numbers.
source(file.path("R", "chiou_mcas.R"))
source(file.path("scenarios", "chiou_competition.R"))

full <- "--full" %in% commandArgs(trailingOnly = TRUE)
p <- chiou_fig6_parameters()
if (!full) {
  p$n <- 60L
  p$equilibration_time <- 20
  p$competition_limit <- 200
  p$output_interval <- 1
  p$initial_dt <- 1e-3
  p$reaction_error <- 0.2
  p$dt_growth <- 1.5
}
template <- read_chiou_initial_state()

mass_grid <- if (full) c(10, 15, 20, 25) else c(10, 20)
mass_runs <- lapply(mass_grid, function(mass) {
  run <- run_chiou_competition(10, mass, p, template)
  data.frame(total_mass = mass, competition_time = run$competition_time,
             censored_at = run$censored_at, u_max = run$u_max,
             saturation_index = run$saturation_index)
})
mass_summary <- do.call(rbind, mass_runs)
print(mass_summary)

design <- if (full) chiou_published_size_design() else
  chiou_published_size_design(c(10, 15))
size_result <- run_chiou_published_size_experiment(
  design, p, template = template,
  equilibration_time = p$equilibration_time,
  competition_limit = p$competition_limit)
print(size_result$summary)

if (full) {
  dir.create(file.path("reports", "output"), recursive = TRUE, showWarnings = FALSE)
  write.csv(mass_summary, file.path("reports", "output", "chiou_mass_validation.csv"),
            row.names = FALSE)
  write.csv(size_result$summary, file.path("reports", "output", "chiou_fig6_reproduction.csv"),
            row.names = FALSE)
}

finite_mass <- mass_summary[is.finite(mass_summary$competition_time), ]
if (nrow(finite_mass) >= 2L) {
  stopifnot(cor(finite_mass$competition_time, finite_mass$saturation_index,
                method = "spearman") > 0)
}
message(if (full) "Full Chiou validation completed." else "Quick Chiou validation completed.")
