# Render public report pages for GitHub Pages.
# Run from the project root:
#   Rscript scripts/render_public_reports.R

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Rendering requires the rmarkdown package.")
}

# R installed outside RStudio may not inherit Pandoc's user-local PATH entry.
pandoc_candidates <- c(
  Sys.getenv("RSTUDIO_PANDOC"),
  dirname(Sys.which(if (.Platform$OS.type == "windows") "pandoc.exe" else "pandoc")),
  file.path(Sys.getenv("LOCALAPPDATA"), "Pandoc")
)
pandoc_candidates <- unique(pandoc_candidates[nzchar(pandoc_candidates)])
pandoc_binary <- if (.Platform$OS.type == "windows") "pandoc.exe" else "pandoc"
pandoc_directory <- pandoc_candidates[file.exists(file.path(pandoc_candidates, pandoc_binary))][1]
if (is.na(pandoc_directory) || !nzchar(pandoc_directory)) {
  stop("Pandoc was not found. Set RSTUDIO_PANDOC to its installation directory.")
}
Sys.setenv(RSTUDIO_PANDOC = pandoc_directory)
if (!rmarkdown::pandoc_available()) stop("rmarkdown could not use Pandoc.")

output_directory <- file.path("docs", "reports")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
if (!file.exists("reports/output/mt_transport_results.rds")) {
  source("scripts/generate_mt_transport_results.R")
}
reports <- sort(list.files("reports", pattern = "\\.Rmd$", full.names = TRUE))
if (!length(reports)) stop("No R Markdown reports were found in reports/.")
for (report in reports) {
  output <- rmarkdown::render(report, output_dir = output_directory, quiet = TRUE)
  message("Rendered ", output)
}
