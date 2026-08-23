directory <- file.path("reports", "output", "chiou_full_reproduction")
fig6 <- read.csv(file.path(directory, "fig6_reproduction_summary.csv"),
                 stringsAsFactors = FALSE)
mass <- read.csv(file.path(directory, "mass_validation_summary.csv"),
                 stringsAsFactors = FALSE)

constant_concentration <- fig6[fig6$control == "constant overall concentration", ]
constant_concentration <- constant_concentration[order(constant_concentration$length), ]
finite_times <- constant_concentration$competition_time[
  is.finite(constant_concentration$competition_time)]
stopifnot(all(diff(finite_times) > 0))
stopifnot(is.na(constant_concentration$competition_time[
  constant_concentration$length == 20]))
stopifnot(constant_concentration$censored_at[
  constant_concentration$length == 20] == 1e5)

fixed_mass <- fig6[fig6$control == "constant total M", ]
stopifnot(fixed_mass$competition_time[fixed_mass$length == 20] <
            fixed_mass$competition_time[fixed_mass$length == 10])

fixed_peak <- fig6[fig6$control == "constant protein content in peaks" &
                     fig6$length >= 10, ]
fixed_peak <- fixed_peak[order(fixed_peak$length), ]
stopifnot(all(diff(fixed_peak$competition_time) > 0))
stopifnot(fixed_peak$competition_time[fixed_peak$length == 20] <
            2 * fixed_peak$competition_time[fixed_peak$length == 10])

mass <- mass[order(mass$total_mass), ]
stopifnot(all(diff(mass$competition_time) > 0))
stopifnot(all(diff(mass$saturation_index) > 0))

message("Full Chiou reproduction checks passed.")
