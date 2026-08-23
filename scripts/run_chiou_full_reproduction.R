# Native-resolution Chiou reproduction and validation suite.
# Independent conditions run in parallel and checkpoint under reports/output/.
source(file.path("R", "chiou_mcas.R"))
source(file.path("scenarios", "chiou_competition.R"))

arguments <- commandArgs(trailingOnly = TRUE)
read_option <- function(prefix, default) {
  hit <- arguments[startsWith(arguments, prefix)]
  if (length(hit)) as.numeric(sub(prefix, "", hit[1], fixed = TRUE)) else default
}
workers <- as.integer(read_option("--workers=", 4))
workers <- max(1L, min(workers, parallel::detectCores(logical = FALSE)))

p <- chiou_fig6_parameters()
template <- read_chiou_initial_state()
published <- chiou_published_size_design(parameters = p)
published$experiment <- "published_size"
validation <- data.frame(
  control = "mass/saturation validation", length = 10,
  total_mass = c(10, 15, 20, 25), experiment = "mass_validation",
  stringsAsFactors = FALSE)
requested <- rbind(published, validation)
requested$key <- sprintf("L%05.2f_M%07.3f", requested$length, requested$total_mass)
unique_tasks <- requested[!duplicated(requested$key), c("key", "length", "total_mass")]

output_dir <- file.path("reports", "output", "chiou_full_reproduction")
checkpoint_dir <- file.path(output_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
unique_tasks$checkpoint <- file.path(checkpoint_dir, paste0(unique_tasks$key, ".rds"))

root <- normalizePath(".", winslash = "/")
cluster <- parallel::makeCluster(workers, outfile = "")
on.exit(parallel::stopCluster(cluster), add = TRUE)
parallel::clusterExport(cluster, c("root", "p", "template"), envir = environment())
parallel::clusterEvalQ(cluster, {
  setwd(root)
  source(file.path("R", "chiou_mcas.R"))
  invisible(.chiou_integrate_fast(rep(1, 4), rep(1, 4), 1, 0.01, p))
  NULL
})

run_task <- function(task) {
  checkpoint <- task$checkpoint
  if (file.exists(checkpoint)) {
    saved <- readRDS(checkpoint)
    message("Reusing ", task$key)
    return(saved$summary)
  }
  message("Starting ", task$key)
  started <- Sys.time()
  equilibration_timing <- system.time(initial <- chiou_competition_initial_state(
    task$length, task$total_mass, p, template,
    equilibration_time = p$equilibration_time, engine = "optimized"))
  competition_timing <- system.time(run <- run_chiou_competition(
    task$length, task$total_mass, p, template,
    equilibration_time = p$equilibration_time,
    competition_limit = p$competition_limit,
    engine = "optimized", initial_state = initial))
  summary <- data.frame(
    key = task$key, length = task$length, total_mass = task$total_mass,
    mean_concentration = run$mean_concentration,
    competition_time = run$competition_time, censored_at = run$censored_at,
    u_max = run$u_max, saturation_index = run$saturation_index,
    peak_heights = paste(signif(run$peak_heights, 10), collapse = ";"),
    peak_widths = paste(signif(run$peak_widths, 10), collapse = ";"),
    peak_separation = run$peak_separation,
    final_number_of_peaks = run$final_number_of_peaks,
    equilibration_steps = initial$accepted_steps,
    competition_steps = run$accepted_steps,
    rejected_steps = run$rejected_steps,
    equilibration_elapsed_s = unname(equilibration_timing["elapsed"]),
    competition_elapsed_s = unname(competition_timing["elapsed"]),
    total_elapsed_s = as.numeric(difftime(Sys.time(), started, units = "secs")),
    stringsAsFactors = FALSE)
  trajectory <- data.frame(
    time = run$time, winner_mass = run$winner_mass,
    loser_mass = run$loser_mass,
    winner_fraction = apply(run$half_fractions, 1, max),
    loser_fraction = apply(run$half_fractions, 1, min))
  saveRDS(list(summary = summary, trajectory = trajectory,
               initial_u = run$initial_u, initial_v = run$initial_v,
               final_u = run$final_u, final_v = run$final_v,
               parameters = p), checkpoint)
  message("Finished ", task$key, " in ", round(summary$total_elapsed_s, 1), " s")
  summary
}

task_list <- lapply(seq_len(nrow(unique_tasks)), function(i)
  as.list(unique_tasks[i, , drop = FALSE]))
summaries <- parallel::parLapplyLB(cluster, task_list, run_task)
summary_table <- do.call(rbind, summaries)

requested_summary <- merge(requested, summary_table, by = c("key", "length", "total_mass"),
                           all.x = TRUE, sort = FALSE)
published_summary <- requested_summary[requested_summary$experiment == "published_size", ]
validation_summary <- requested_summary[requested_summary$experiment == "mass_validation", ]
write.csv(summary_table, file.path(output_dir, "unique_run_summary.csv"), row.names = FALSE)
write.csv(published_summary, file.path(output_dir, "fig6_reproduction_summary.csv"), row.names = FALSE)
write.csv(validation_summary, file.path(output_dir, "mass_validation_summary.csv"), row.names = FALSE)

manifest <- list(
  completed_at = Sys.time(), workers = workers,
  parameters = p, requested_runs = nrow(requested), unique_runs = nrow(unique_tasks),
  total_worker_seconds = sum(summary_table$total_elapsed_s),
  outputs = c("unique_run_summary.csv", "fig6_reproduction_summary.csv",
              "mass_validation_summary.csv"))
saveRDS(manifest, file.path(output_dir, "manifest.rds"))
message("Full Chiou suite complete: ", nrow(unique_tasks), " unique runs on ",
        workers, " workers; ", round(manifest$total_worker_seconds, 1),
        " cumulative worker-seconds.")
