# Refresh cargo trajectories and membrane-distribution summaries while reusing
# the expensive cached networks and junction detections. Run from project root.

source("R/mt_network.R")
source("R/mt_transport.R")

cache_path <- "reports/output/mt_transport_results.rds"
results <- readRDS(cache_path)
variants <- names(results$network_sets)
sizes <- c("1x volume", "2x volume")
transport_summary <- data.frame()
membrane_bins <- data.frame()
membrane_metrics <- data.frame()
membrane_summaries <- list()

job_id <- 0L
for (variant in variants) {
  membrane_summaries[[variant]] <- list()
  for (size in sizes) {
    job_id <- job_id + 1L
    message("Refreshing cargo: ", variant, ", ", size)
    network <- results$network_sets[[variant]][[size]]
    cargo <- simulate_mt_transport(network, results$junctions[[variant]][[size]],
      n_cargo = 250, source_radius_um = 0.75, capture_radius_um = 0.4,
      hop_radius_um = 0.5, speed_um_s = 1, deadline_s = 120,
      cortex_distance_um = 0.3, hop_search_time_s = 1,
      pause_probability = 0.5, max_hops = 30, seed = 10000 + job_id)
    results$transport[[variant]][[size]] <- cargo
    transport_summary <- rbind(transport_summary,
      cbind(variant = variant, size = size, summarize_mt_transport(cargo)))
    membrane <- summarize_membrane_distribution(cargo, network$geometry)
    membrane_summaries[[variant]][[size]] <- membrane
    membrane_bins <- rbind(membrane_bins,
      cbind(variant = variant, size = size, membrane$bins))
    membrane_metrics <- rbind(membrane_metrics,
      cbind(variant = variant, size = size, membrane$metrics))
  }
}

results$transport_summary <- transport_summary
results$membrane_bins <- membrane_bins
results$membrane_metrics <- membrane_metrics
results$membrane_size_divergence <- do.call(rbind, lapply(variants, function(variant) {
  cbind(variant = variant, compare_membrane_distributions(
    membrane_summaries[[variant]][["1x volume"]],
    membrane_summaries[[variant]][["2x volume"]]))
}))
saveRDS(results, cache_path)
message("Updated ", cache_path)
