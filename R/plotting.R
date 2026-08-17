# Small ggplot2 helpers shared across reports.

plot_radial_profile <- function(grid, concentration, main = "Radial profile",
                                xlab = "Radius (um)", ylab = "Concentration") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plotting helpers require the ggplot2 package.")
  }
  ggplot2::ggplot(
    data.frame(radius_um = grid$centers_um, concentration = concentration),
    ggplot2::aes(x = radius_um, y = concentration)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(title = main, x = xlab, y = ylab) +
    ggplot2::theme_minimal()
}

plot_scaling_curve <- function(volume_multiplier, value, main = "Scaling curve",
                               ylab = "Relative value") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plotting helpers require the ggplot2 package.")
  }
  ggplot2::ggplot(
    data.frame(volume_multiplier = volume_multiplier, value = value),
    ggplot2::aes(x = volume_multiplier, y = value)
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point() +
    ggplot2::labs(title = main, x = "Cell-volume multiplier", y = ylab) +
    ggplot2::theme_minimal()
}
