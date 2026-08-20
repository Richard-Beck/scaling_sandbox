reports <- list.files("reports", pattern = "\\.Rmd$", full.names = TRUE)
patterns <- c("simulate_vcell\\s*\\(", "compile_vcell_kernel\\s*\\(",
              "generate_wgd_spatial_results\\.R")
for (report in reports) {
  text <- paste(readLines(report, warn = FALSE), collapse = "\n")
  hits <- patterns[vapply(patterns, grepl, logical(1), x = text, perl = TRUE)]
  if (length(hits)) stop(report, " invokes the WGD simulation path: ", paste(hits, collapse = ", "))
}
message("Report sources do not invoke the WGD PDE simulation path.")
