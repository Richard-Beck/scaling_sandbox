# Resolve this relocated bundle for Rscript and source()/agentRrunner calls.
.chiou_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.chiou_script <- if (length(.chiou_source_files)) {
  tail(.chiou_source_files, 1)[[1]]
} else {
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
}
if (is.na(.chiou_script)) stop("Run this file with Rscript or source().")
.chiou_bundle <- normalizePath(file.path(dirname(.chiou_script), ".."), mustWork = TRUE)
setwd(.chiou_bundle)
rm(.chiou_source_files, .chiou_script)

source(file.path("src", "chiou_mcas.R"))
source(file.path("params-or-config", "chiou_competition.R"))

p <- chiou_fig6_parameters(); state <- read_chiou_initial_state()
stopifnot(length(state$u) == 500L, length(state$v) == 500L)
stopifnot(abs(chiou_mass(state$u, state$v, 10) - 20) < 1e-8)
stopifnot(abs(chiou_reaction(2, 3, p) - 10) < 1e-12)
sat <- chiou_saturation_values(p)
stopifnot(abs(sat$u_sat - sqrt(200)) < 1e-12, abs(sat$q_sat - sqrt(0.045)) < 1e-12)
wgd <- chiou_wgd_adapter(2, abundance_exponent = 1)
stopifnot(abs(wgd$length_ratio - 2^(1/3)) < 1e-12)
design <- chiou_published_size_design(c(10, 20))
stopifnot(design$total_mass[design$control == "constant overall concentration"][2] == 40)
stopifnot(design$total_mass[design$control == "constant total M"][2] == 20)
stopifnot(abs(design$total_mass[design$control == "constant protein content in peaks"][2] -
                (20 + 10 * sat$q_sat)) < 1e-12)

# Short conservative numerical smoke test, not a reproduction run.
p$n <- 40L; p$initial_dt <- 1e-3; p$output_interval <- 0.25
u <- rep(1, p$n); v <- rep(1, p$n); u[10] <- 1.1; u[30] <- 0.9
before <- chiou_mass(u, v, 10)
out <- .chiou_integrate(u, v, 10, 1, p, record = TRUE)
after <- chiou_mass(out$u, out$v, 10)
stopifnot(abs(after - before) < 1e-8, all(is.finite(c(out$u, out$v))))
message("Chiou port checks passed.")
