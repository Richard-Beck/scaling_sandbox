source(file.path("R", "chiou_mcas.R"))
source(file.path("scenarios", "chiou_competition.R"))

arguments <- commandArgs(trailingOnly = TRUE)
read_option <- function(prefix, default) {
  hit <- arguments[startsWith(arguments, prefix)]
  if (length(hit)) as.numeric(sub(prefix, "", hit[1], fixed = TRUE)) else default
}
n <- as.integer(read_option("--n=", 125))
competition_limit <- read_option("--limit=", 5000)
equilibration_time <- read_option("--equilibration=", 2000)

p <- chiou_fig6_parameters()
p$n <- n
p$competition_limit <- competition_limit
p$equilibration_time <- equilibration_time
template <- read_chiou_initial_state()
invisible(.chiou_integrate_fast(rep(1, 4), rep(1, 4), 1, 0.01, p))

design <- chiou_published_size_design(parameters = p)
initial_cache <- new.env(parent = emptyenv())
rows <- vector("list", nrow(design))
for (i in seq_len(nrow(design))) {
  key <- sprintf("L%.12g_M%.12g", design$length[i], design$total_mass[i])
  equilibration_elapsed <- 0
  if (!exists(key, envir = initial_cache, inherits = FALSE)) {
    timing <- system.time(initial <- chiou_competition_initial_state(
      design$length[i], design$total_mass[i], p, template,
      equilibration_time, engine = "optimized"))
    equilibration_elapsed <- unname(timing["elapsed"])
    assign(key, initial, envir = initial_cache)
  }
  initial <- get(key, envir = initial_cache)
  timing <- system.time(run <- run_chiou_competition(
    design$length[i], design$total_mass[i], p, template,
    equilibration_time, competition_limit, engine = "optimized",
    initial_state = initial))
  rows[[i]] <- data.frame(
    control = design$control[i], length = design$length[i],
    total_mass = design$total_mass[i], n = n,
    equilibration_time = equilibration_time,
    equilibration_elapsed_s = equilibration_elapsed,
    competition_limit = competition_limit,
    simulated_competition_time = tail(run$time, 1),
    competition_time = run$competition_time,
    censored_at = run$censored_at,
    competition_elapsed_s = unname(timing["elapsed"]),
    accepted_steps = run$accepted_steps,
    rejected_steps = run$rejected_steps,
    saturation_index = run$saturation_index
  )
  print(rows[[i]])
}
benchmark <- do.call(rbind, rows)
benchmark$linear_n500_projection_s <- benchmark$competition_elapsed_s * 500 / n

dir.create(file.path("reports", "output"), recursive = TRUE, showWarnings = FALSE)
path <- file.path("reports", "output", sprintf("chiou_runtime_n%d.csv", n))
write.csv(benchmark, path, row.names = FALSE)

cat("\nSequential observed elapsed:", sum(benchmark$equilibration_elapsed_s +
                                           benchmark$competition_elapsed_s), "seconds\n")
cat("Linear-in-N projection of competition portions at N=500:",
    sum(benchmark$linear_n500_projection_s), "seconds\n")
cat("Saved", path, "\n")
