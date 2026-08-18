# Shariff-style 3D microtubule network scenario.

mt_polymer_density_scenario <- function() {
  list(
    name = "3D Shariff-style MT networks at matched polymer density",
    geometry_3d = mt_geometry(c(18, 12, 3.5), c(5.2, 4.0, 2.1)),
    parameters = list(
      n_base = 175L,
      step_um = 0.2,
      mean_length_um = 35,
      sd_length_um = 15,
      cos_alpha = 0.9,
      centrosome_diameter_um = 0.4,
      max_direction_attempts = 250L,
      pruning_iterations = 10L,
      pruning_removal_fraction = 0.5,
      pruning_kernel_um = 1.5,
      pruning_grid_n = c(24L, 16L, 8L)
    )
  )
}
