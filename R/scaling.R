# Scaling relationships shared across scenarios.

#' Convert a volume multiplier to a multiplier for a quantity with exponent.
#'
#' @param volume_multiplier Cell-volume change relative to baseline.
#' @param exponent Scaling exponent; 1 for volume, 2/3 for area, 1/3 for length.
#' @return A numeric multiplier.
scale_factor <- function(volume_multiplier, exponent) {
  stopifnot(all(volume_multiplier > 0), length(exponent) == 1L)
  volume_multiplier ^ exponent
}

#' Geometry for a spherical cell scaled from a baseline radius.
spherical_geometry <- function(radius_um, volume_multiplier = 1) {
  stopifnot(radius_um > 0, all(volume_multiplier > 0))
  radius <- radius_um * scale_factor(volume_multiplier, 1 / 3)
  list(
    radius_um = radius,
    volume_um3 = 4 / 3 * pi * radius ^ 3,
    area_um2 = 4 * pi * radius ^ 2
  )
}

#' Relative source loading against downstream processing capacity.
relative_loading <- function(volume_multiplier, source_exponent, capacity_exponent) {
  scale_factor(volume_multiplier, source_exponent - capacity_exponent)
}

#' Michaelis-Menten substrate concentration in a well-mixed steady state.
#'
#' Returns Inf when source exceeds the maximum processing capacity.
mm_steady_state <- function(source_rate, vmax, km) {
  stopifnot(source_rate >= 0, vmax > 0, km > 0)
  if (source_rate >= vmax) return(Inf)
  km * source_rate / (vmax - source_rate)
}
