# Generate cached network, junction, and transport results for the MT report.
# Run from the project root.

source("R/mt_network.R")
source("R/mt_transport.R")
source("scenarios/mt_polymer_density.R")

scenario <- mt_polymer_density_scenario()
p <- scenario$parameters
generator_names <- c("step_um", "mean_length_um", "sd_length_um", "cos_alpha",
                     "centrosome_diameter_um", "max_direction_attempts")
generator_parameters <- p[generator_names]

base <- do.call(generate_matched_mt_networks,
  c(list(base_geometry = scenario$geometry_3d, volume_multiplier = 2,
         n_base = p$n_base, seed = 9101), generator_parameters))

prune_pair <- function(score_mode, seed) {
  list(
    `1x volume` = prune_mt_network(base[["1x volume"]],
      iterations = p$pruning_iterations, score_mode = score_mode,
      removal_fraction = p$pruning_removal_fraction,
      density_kernel_um = p$pruning_kernel_um, grid_n = p$pruning_grid_n,
      refill_to_target = TRUE, seed = seed),
    `2x volume` = prune_mt_network(base[["2x volume"]],
      iterations = p$pruning_iterations, score_mode = score_mode,
      removal_fraction = p$pruning_removal_fraction,
      density_kernel_um = p$pruning_kernel_um, grid_n = p$pruning_grid_n,
      refill_to_target = TRUE, seed = seed + 1L)
  )
}

li_pair <- list(
  `1x volume` = do.call(generate_mt_network,
    c(list(geometry = scenario$geometry_3d,
           target_density = base[["1x volume"]]$target_density,
           boundary_mode = "li", seed = 9401), generator_parameters)),
  `2x volume` = do.call(generate_mt_network,
    c(list(geometry = scale_mt_geometry(scenario$geometry_3d, 2),
           target_density = base[["1x volume"]]$target_density,
           boundary_mode = "li", seed = 9402), generator_parameters))
)

network_sets <- list(
  `Unpruned Shariff` = base,
  `Integrated-density pruning` = prune_pair("integrated", 9201),
  `Distal-tip pruning` = prune_pair("distal", 9301),
  `Li rebound/regeneration` = li_pair
)

junctions <- list()
transport <- list()
junction_sensitivity <- data.frame()
transport_summary <- data.frame()
network_summary <- data.frame()
membrane_bins <- data.frame()
membrane_metrics <- data.frame()
membrane_summaries <- list()

job_id <- 0L
for (variant in names(network_sets)) {
  junctions[[variant]] <- list()
  transport[[variant]] <- list()
  membrane_summaries[[variant]] <- list()
  for (size in names(network_sets[[variant]])) {
    job_id <- job_id + 1L
    network <- network_sets[[variant]][[size]]
    message("Detecting junctions: ", variant, ", ", size)
    max_detection <- detect_mt_junctions(network,
      encounter_radius_um = 0.5, search_radius_um = 0.5,
      merge_radius_um = 0.5, parallel_angle_deg = 10,
      bundle_min_length_um = 1, exclude_mtoc_radius_um = 0.75)
    active_detection <- subset_mt_junctions(max_detection, 0.25, 0.5)
    junctions[[variant]][[size]] <- active_detection
    for (radius in c(0.1, 0.25, 0.5)) {
      sensitivity <- subset_mt_junctions(max_detection, radius, 0.5)
      junction_sensitivity <- rbind(junction_sensitivity, data.frame(
        variant = variant, size = size, encounter_radius_um = radius,
        junction_count = nrow(sensitivity$junctions),
        high_multiplicity_fraction = if (nrow(sensitivity$junctions))
          mean(sensitivity$junctions$n_microtubules > 5) else 0,
        mean_multiplicity = if (nrow(sensitivity$junctions))
          mean(sensitivity$junctions$n_microtubules) else NA_real_))
    }
    message("Simulating cargo: ", variant, ", ", size)
    cargo <- simulate_mt_transport(network, active_detection,
      n_cargo = 250, source_radius_um = 0.75, capture_radius_um = 0.4,
      hop_radius_um = 0.5, speed_um_s = 1, deadline_s = 120,
      cortex_distance_um = 0.3, hop_search_time_s = 1,
      pause_probability = 0.5, max_hops = 30, seed = 10000 + job_id)
    transport[[variant]][[size]] <- cargo
    membrane <- summarize_membrane_distribution(cargo, network$geometry)
    membrane_summaries[[variant]][[size]] <- membrane
    membrane_bins <- rbind(membrane_bins,
      cbind(variant = variant, size = size, membrane$bins))
    membrane_metrics <- rbind(membrane_metrics,
      cbind(variant = variant, size = size, membrane$metrics))
    transport_summary <- rbind(transport_summary,
      cbind(variant = variant, size = size, summarize_mt_transport(cargo)))
    network_summary <- rbind(network_summary,
      cbind(variant = variant, size = size, summarize_mt_network(network),
        junction_count = nrow(active_detection$junctions),
        mean_junction_multiplicity = if (nrow(active_detection$junctions))
          mean(active_detection$junctions$n_microtubules) else NA_real_,
        high_multiplicity_fraction = if (nrow(active_detection$junctions))
          mean(active_detection$junctions$n_microtubules > 5) else 0,
        bundle_count = nrow(active_detection$bundles)))
  }
}

membrane_size_divergence <- do.call(rbind, lapply(names(network_sets), function(variant) {
  cbind(variant = variant, compare_membrane_distributions(
    membrane_summaries[[variant]][["1x volume"]],
    membrane_summaries[[variant]][["2x volume"]]))
}))

results <- list(scenario = scenario, network_sets = network_sets,
                junctions = junctions, transport = transport,
                network_summary = network_summary,
                transport_summary = transport_summary,
                junction_sensitivity = junction_sensitivity,
                membrane_bins = membrane_bins,
                membrane_metrics = membrane_metrics,
                membrane_size_divergence = membrane_size_divergence,
                transport_parameters = list(n_cargo = 250,
                  encounter_radius_um = 0.25, source_radius_um = 0.75,
                  capture_radius_um = 0.4, hop_radius_um = 0.5,
                  speed_um_s = 1, deadline_s = 120,
                  cortex_distance_um = 0.3, pause_probability = 0.5))

dir.create("reports/output", recursive = TRUE, showWarnings = FALSE)
saveRDS(results, "reports/output/mt_transport_results.rds")
message("Saved reports/output/mt_transport_results.rds")
