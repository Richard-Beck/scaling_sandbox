wgd_spatial_design <- function() {
  sizes <- c(0.5, 1, 2)
  scale_from_2x <- sizes^(1 / 3) / 2^(1 / 3)
  data.frame(
    geometry = c(
      "AR 1.00 - 0.5x", "AR 1.00 - 1x", "AR 1.00 - 2x",
      "AR 1.29 - 0.5x", "AR 1.29 - 1x", "AR 1.29 - 2x"
    ),
    pair = c(rep("L:W = 90:90", 3), rep("L:W = 90:70", 3)),
    size = rep(c("0.5x", "1x", "2x"), 2),
    length_um = 90 * rep(scale_from_2x, 2),
    width_um = c(90 * scale_from_2x, 70 * scale_from_2x),
    volume_ratio = rep(sizes, 2),
    stringsAsFactors = FALSE
  )
}

wgd_spatial_parameters <- function() {
  list(
    simulation = "simulation2",
    mesh = c(x = 80L, y = 80L),
    start_time_s = 0,
    end_time_s = 200,
    timestep_s = 0.02,
    save_interval_s = 1,
    stimulus_end_s = 10,
    selected_spatial_times_s = c(10, 20, 50, 100, 150, 200),
    source_vcml = "scenarios/wgd_spatial_scale.vcml"
  )
}
