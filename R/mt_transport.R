# Transport over static 3D microtubule networks.

.segment_closest_points <- function(p1, q1, p2, q2) {
  u <- q1 - p1
  v <- q2 - p2
  w <- p1 - p2
  a <- sum(u * u)
  b <- sum(u * v)
  c <- sum(v * v)
  d <- sum(u * w)
  e <- sum(v * w)
  denominator <- a * c - b * b
  small <- 1e-12
  s_n <- denominator
  s_d <- denominator
  t_n <- denominator
  t_d <- denominator
  if (denominator < small) {
    s_n <- 0
    s_d <- 1
    t_n <- e
    t_d <- c
  } else {
    s_n <- b * e - c * d
    t_n <- a * e - b * d
    if (s_n < 0) {
      s_n <- 0
      t_n <- e
      t_d <- c
    } else if (s_n > s_d) {
      s_n <- s_d
      t_n <- e + b
      t_d <- c
    }
  }
  if (t_n < 0) {
    t_n <- 0
    if (-d < 0) s_n <- 0 else if (-d > a) s_n <- s_d else {
      s_n <- -d
      s_d <- a
    }
  } else if (t_n > t_d) {
    t_n <- t_d
    if (-d + b < 0) s_n <- 0 else if (-d + b > a) s_n <- s_d else {
      s_n <- -d + b
      s_d <- a
    }
  }
  s <- if (abs(s_n) < small) 0 else s_n / s_d
  t <- if (abs(t_n) < small) 0 else t_n / t_d
  point_1 <- p1 + s * u
  point_2 <- p2 + t * v
  list(distance = sqrt(sum((point_1 - point_2) ^ 2)), point_1 = point_1,
       point_2 = point_2, s = s, t = t)
}

.transport_segments <- function(network) {
  segments <- mt_segments(network)
  segments$row_id <- seq_len(nrow(segments))
  segments$arc_start_um <- ave(segments$length_um, segments$filament,
    FUN = function(lengths) c(0, head(cumsum(lengths), -1L)))
  segments
}

.cluster_pair_encounters <- function(pair_encounters, merge_radius_um) {
  if (!nrow(pair_encounters)) {
    return(list(junctions = data.frame(), memberships = data.frame()))
  }
  clusters <- list()
  assignment <- integer(nrow(pair_encounters))
  for (i in seq_len(nrow(pair_encounters))) {
    point <- as.numeric(pair_encounters[i, c("x", "y", "z")])
    if (!length(clusters)) {
      clusters[[1]] <- list(center = point, rows = i)
      assignment[i] <- 1L
      next
    }
    distances <- vapply(clusters, function(cluster) {
      sqrt(sum((point - cluster$center) ^ 2))
    }, numeric(1))
    nearest <- which.min(distances)
    if (distances[nearest] <= merge_radius_um) {
      clusters[[nearest]]$rows <- c(clusters[[nearest]]$rows, i)
      clusters[[nearest]]$center <- colMeans(pair_encounters[
        clusters[[nearest]]$rows, c("x", "y", "z"), drop = FALSE])
      assignment[i] <- nearest
    } else {
      clusters[[length(clusters) + 1L]] <- list(center = point, rows = i)
      assignment[i] <- length(clusters)
    }
  }

  junction_rows <- list()
  membership_rows <- list()
  for (junction_id in seq_along(clusters)) {
    rows <- clusters[[junction_id]]$rows
    encounters <- pair_encounters[rows, , drop = FALSE]
    mt_ids <- sort(unique(c(encounters$filament_1, encounters$filament_2)))
    junction_rows[[junction_id]] <- data.frame(junction_id = junction_id,
      x = mean(encounters$x), y = mean(encounters$y), z = mean(encounters$z),
      n_microtubules = length(mt_ids),
      minimum_3d_separation_um = min(encounters$distance_um),
      mean_crossing_angle_deg = mean(encounters$angle_deg),
      max_crossing_angle_deg = max(encounters$angle_deg),
      mt_ids = I(list(mt_ids)))
    memberships <- rbind(
      data.frame(filament = encounters$filament_1, arc_um = encounters$arc_1_um),
      data.frame(filament = encounters$filament_2, arc_um = encounters$arc_2_um))
    memberships <- stats::aggregate(arc_um ~ filament, memberships, mean)
    memberships$junction_id <- junction_id
    membership_rows[[junction_id]] <- memberships[, c("junction_id", "filament", "arc_um")]
  }
  list(junctions = do.call(rbind, junction_rows),
       memberships = do.call(rbind, membership_rows))
}

#' Detect transport-relevant 3D MT encounters and merge them into junctions.
#'
#' A fast midpoint-radius search supplies candidate segment pairs; exact 3D
#' segment distances are then calculated. Sustained near-parallel encounters
#' are recorded as bundles and excluded from the junction table.
detect_mt_junctions <- function(network, encounter_radius_um = 0.25,
                                search_radius_um = 0.5,
                                merge_radius_um = 0.5,
                                parallel_angle_deg = 10,
                                bundle_min_length_um = 1,
                                exclude_mtoc_radius_um = 0.75,
                                max_neighbors = 80L) {
  if (!requireNamespace("RANN", quietly = TRUE)) stop("RANN is required")
  stopifnot(encounter_radius_um > 0, search_radius_um >= encounter_radius_um,
            merge_radius_um > 0, parallel_angle_deg >= 0)
  segments <- .transport_segments(network)
  midpoints <- cbind((segments$x + segments$xend) / 2,
                     (segments$y + segments$yend) / 2,
                     (segments$z + segments$zend) / 2)
  expanded_radius <- search_radius_um + max(segments$length_um)
  neighbors <- RANN::nn2(midpoints, midpoints,
    k = min(max_neighbors, nrow(midpoints)), searchtype = "radius",
    radius = expanded_radius)$nn.idx
  candidate_key_list <- vector("list", nrow(neighbors))
  for (i in seq_len(nrow(neighbors))) {
    js <- neighbors[i, ]
    js <- js[js > i]
    if (length(js)) {
      js <- js[segments$filament[js] != segments$filament[i]]
      if (length(js)) candidate_key_list[[i]] <- paste(i, js, sep = ":")
    }
  }
  candidate_keys <- unique(unlist(candidate_key_list, use.names = FALSE))
  if (!length(candidate_keys)) {
    return(list(junctions = data.frame(), memberships = data.frame(),
                bundles = data.frame(), pair_encounters = data.frame()))
  }
  pairs <- do.call(rbind, strsplit(candidate_keys, ":", fixed = TRUE))
  pairs <- matrix(as.integer(pairs), ncol = 2L)
  pair_labels <- paste(pmin(segments$filament[pairs[, 1]], segments$filament[pairs[, 2]]),
                       pmax(segments$filament[pairs[, 1]], segments$filament[pairs[, 2]]),
                       sep = ":")
  midpoint_distances <- sqrt(rowSums((midpoints[pairs[, 1], , drop = FALSE] -
                                       midpoints[pairs[, 2], , drop = FALSE]) ^ 2))
  pair_groups_raw <- split(seq_len(nrow(pairs)), pair_labels)
  candidate_extent <- vapply(pair_groups_raw, function(rows) {
    length(unique(pairs[rows, 1])) * mean(segments$length_um[pairs[rows, 1]])
  }, numeric(1))
  # Exact segment geometry is evaluated only for the two closest midpoint
  # candidates per MT pair. This retains the pairwise minimum calculation while
  # avoiding an exact calculation for every nearby 0.2-um segment combination.
  shortlisted <- unlist(lapply(pair_groups_raw, function(rows) {
    head(rows[order(midpoint_distances[rows])], 2L)
  }), use.names = FALSE)
  pairs <- pairs[shortlisted, , drop = FALSE]
  close_rows <- list()
  used <- 0L
  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k, 1]
    j <- pairs[k, 2]
    closest <- .segment_closest_points(
      c(segments$x[i], segments$y[i], segments$z[i]),
      c(segments$xend[i], segments$yend[i], segments$zend[i]),
      c(segments$x[j], segments$y[j], segments$z[j]),
      c(segments$xend[j], segments$yend[j], segments$zend[j]))
    if (closest$distance <= search_radius_um) {
      tangent_1 <- c(segments$xend[i] - segments$x[i],
                     segments$yend[i] - segments$y[i],
                     segments$zend[i] - segments$z[i]) / segments$length_um[i]
      tangent_2 <- c(segments$xend[j] - segments$x[j],
                     segments$yend[j] - segments$y[j],
                     segments$zend[j] - segments$z[j]) / segments$length_um[j]
      used <- used + 1L
      close_rows[[used]] <- data.frame(segment_1 = i, segment_2 = j,
        filament_1 = segments$filament[i], filament_2 = segments$filament[j],
        distance_um = closest$distance,
        x = mean(c(closest$point_1[1], closest$point_2[1])),
        y = mean(c(closest$point_1[2], closest$point_2[2])),
        z = mean(c(closest$point_1[3], closest$point_2[3])),
        angle_deg = acos(pmin(1, abs(sum(tangent_1 * tangent_2)))) * 180 / pi,
        arc_1_um = segments$arc_start_um[i] + closest$s * segments$length_um[i],
        arc_2_um = segments$arc_start_um[j] + closest$t * segments$length_um[j])
    }
  }
  close <- if (length(close_rows)) do.call(rbind, close_rows) else data.frame()
  if (!nrow(close)) {
    return(list(junctions = data.frame(), memberships = data.frame(),
                bundles = data.frame(), pair_encounters = data.frame()))
  }
  close$pair <- paste(pmin(close$filament_1, close$filament_2),
                      pmax(close$filament_1, close$filament_2), sep = ":")
  pair_groups <- split(seq_len(nrow(close)), close$pair)
  pair_rows <- list()
  bundle_rows <- list()
  for (pair_name in names(pair_groups)) {
    rows <- pair_groups[[pair_name]]
    eligible <- rows[close$distance_um[rows] <= encounter_radius_um]
    if (!length(eligible)) next
    best <- eligible[which.min(close$distance_um[eligible])]
    extent <- candidate_extent[[pair_name]]
    is_bundle <- close$angle_deg[best] < parallel_angle_deg &&
      extent >= bundle_min_length_um
    if (is_bundle) {
      bundle_rows[[length(bundle_rows) + 1L]] <- data.frame(
        filament_1 = close$filament_1[best], filament_2 = close$filament_2[best],
        minimum_3d_separation_um = close$distance_um[best],
        angle_deg = close$angle_deg[best], close_extent_um = extent)
    } else {
      pair_rows[[length(pair_rows) + 1L]] <- close[best, setdiff(names(close), "pair")]
    }
  }
  pair_encounters <- if (length(pair_rows)) do.call(rbind, pair_rows) else data.frame()
  if (nrow(pair_encounters) && exclude_mtoc_radius_um > 0) {
    encounter_points <- as.matrix(pair_encounters[, c("x", "y", "z")])
    mtoc_distances <- sqrt(rowSums((encounter_points - matrix(network$geometry$mtoc_um,
      nrow(pair_encounters), 3L, byrow = TRUE)) ^ 2))
    pair_encounters <- pair_encounters[mtoc_distances > exclude_mtoc_radius_um,
                                       , drop = FALSE]
  }
  clustered <- .cluster_pair_encounters(pair_encounters, merge_radius_um)
  list(junctions = clustered$junctions, memberships = clustered$memberships,
       bundles = if (length(bundle_rows)) do.call(rbind, bundle_rows) else data.frame(),
       pair_encounters = pair_encounters,
       parameters = list(encounter_radius_um = encounter_radius_um,
         search_radius_um = search_radius_um, merge_radius_um = merge_radius_um,
         parallel_angle_deg = parallel_angle_deg,
         bundle_min_length_um = bundle_min_length_um,
         exclude_mtoc_radius_um = exclude_mtoc_radius_um))
}

#' Rebuild junctions at a smaller encounter radius from a max-radius detection.
subset_mt_junctions <- function(detection, encounter_radius_um,
                                merge_radius_um = 0.5) {
  pairs <- detection$pair_encounters
  if (nrow(pairs)) pairs <- pairs[pairs$distance_um <= encounter_radius_um, , drop = FALSE]
  clustered <- .cluster_pair_encounters(pairs, merge_radius_um)
  list(junctions = clustered$junctions, memberships = clustered$memberships,
       bundles = detection$bundles, pair_encounters = pairs,
       parameters = utils::modifyList(detection$parameters,
                                      list(encounter_radius_um = encounter_radius_um,
                                           merge_radius_um = merge_radius_um)))
}

.nearest_track <- function(point, segments, exclude_filament = integer(),
                           minimum_forward_um = 0.1) {
  keep <- !(segments$filament %in% exclude_filament)
  keep <- keep & segments$length_um >= minimum_forward_um
  candidates <- segments[keep, , drop = FALSE]
  if (!nrow(candidates)) return(NULL)
  starts <- as.matrix(candidates[, c("x", "y", "z")])
  ends <- as.matrix(candidates[, c("xend", "yend", "zend")])
  vectors <- ends - starts
  fractions <- rowSums((matrix(point, nrow(candidates), 3, byrow = TRUE) - starts) * vectors) /
    rowSums(vectors ^ 2)
  fractions <- pmax(0, pmin(1, fractions))
  closest <- starts + fractions * vectors
  distances <- sqrt(rowSums((closest - matrix(point, nrow(candidates), 3,
                                               byrow = TRUE)) ^ 2))
  best <- which.min(distances)
  list(filament = candidates$filament[best],
       arc_um = candidates$arc_start_um[best] + fractions[best] * candidates$length_um[best],
       point = closest[best, ], distance_um = distances[best])
}

.sample_cargo_source <- function(geometry, radius_um) {
  for (attempt in seq_len(10000L)) {
    point <- .sample_point_in_sphere(geometry$mtoc_um, radius_um)
    if (.inside_cytoplasm(point, geometry)) return(point)
  }
  stop("Could not sample cargo source region")
}

.cell_boundary_distance <- function(point, geometry) {
  level <- .ellipsoid_level(point, geometry$cell_axes_um)
  if (level <= 0) return(min(geometry$cell_axes_um))
  sqrt(sum((point / sqrt(level) - point) ^ 2))
}

.filament_contact_arcs <- function(mt, geometry, cortex_distance_um) {
  points <- mt$points
  arcs <- c(0, cumsum(sqrt(rowSums((points[-1, , drop = FALSE] -
                                     points[-nrow(points), , drop = FALSE]) ^ 2))))
  distances <- apply(points, 1, .cell_boundary_distance, geometry = geometry)
  arcs[distances <= cortex_distance_um]
}

.draw_junction_delay <- function(multiplicity) {
  if (multiplicity <= 2) return(stats::rlnorm(1, log(1.5), 0.25))
  if (multiplicity <= 5) return(stats::runif(1, 3, 10))
  stats::rlnorm(1, log(20), 0.8)
}

#' Simulate simple plus-end-directed cargo transport on a static MT network.
simulate_mt_transport <- function(network, junction_detection, n_cargo = 250L,
                                  source_radius_um = 0.75,
                                  capture_radius_um = 0.4,
                                  hop_radius_um = 0.5,
                                  speed_um_s = 1,
                                  deadline_s = 120,
                                  cortex_distance_um = 0.3,
                                  hop_search_time_s = 1,
                                  pause_probability = 0.5,
                                  max_hops = 30L, seed = 1L) {
  stopifnot(n_cargo > 0, source_radius_um > 0, capture_radius_um > 0,
            hop_radius_um > 0, speed_um_s > 0, deadline_s > 0,
            pause_probability >= 0, pause_probability <= 1)
  set.seed(seed)
  segments <- .transport_segments(network)
  junctions <- junction_detection$junctions
  memberships <- junction_detection$memberships
  contact_arcs <- lapply(network$filaments, .filament_contact_arcs,
                         geometry = network$geometry,
                         cortex_distance_um = cortex_distance_um)
  results <- vector("list", n_cargo)

  for (cargo_id in seq_len(n_cargo)) {
    source <- .sample_cargo_source(network$geometry, source_radius_um)
    attachment <- .nearest_track(source, segments)
    elapsed <- 0
    path_length <- 0
    hops <- 0L
    junction_count <- 0L
    max_junction_multiplicity <- 0L
    visited_junctions <- integer()
    success <- FALSE
    final_point <- source

    if (!is.null(attachment) && attachment$distance_um <= capture_radius_um) {
      filament <- attachment$filament
      arc <- attachment$arc_um
      path_length <- attachment$distance_um
      final_point <- attachment$point
      repeat {
        mt <- network$filaments[[filament]]
        length_um <- mt$realized_length_um
        contacts <- contact_arcs[[filament]]
        contacts <- contacts[contacts >= arc - 1e-9]
        target_arc <- if (length(contacts)) min(contacts) else length_um

        if (nrow(memberships)) {
          route_junctions <- memberships[memberships$filament == filament &
            memberships$arc_um > arc & memberships$arc_um <= target_arc,
            , drop = FALSE]
          route_junctions <- route_junctions[order(route_junctions$arc_um), , drop = FALSE]
          for (row in seq_len(nrow(route_junctions))) {
            junction_id <- route_junctions$junction_id[row]
            if (!(junction_id %in% visited_junctions)) {
              visited_junctions <- c(visited_junctions, junction_id)
              junction_count <- junction_count + 1L
              multiplicity <- junctions$n_microtubules[
                match(junction_id, junctions$junction_id)]
              max_junction_multiplicity <- max(max_junction_multiplicity, multiplicity)
              if (stats::runif(1) < pause_probability) {
                elapsed <- elapsed + .draw_junction_delay(multiplicity)
              }
            }
          }
        }
        travelled <- max(0, target_arc - arc)
        elapsed <- elapsed + travelled / speed_um_s
        path_length <- path_length + travelled
        points <- mt$points
        point_arcs <- c(0, cumsum(sqrt(rowSums((points[-1, , drop = FALSE] -
                                                 points[-nrow(points), , drop = FALSE]) ^ 2))))
        final_point <- points[which.min(abs(point_arcs - target_arc)), ]
        if (elapsed > deadline_s) break
        if (length(contacts)) {
          success <- TRUE
          break
        }
        if (hops >= max_hops) break
        next_attachment <- .nearest_track(final_point, segments,
                                          exclude_filament = filament)
        if (is.null(next_attachment) || next_attachment$distance_um > hop_radius_um) break
        hops <- hops + 1L
        elapsed <- elapsed + hop_search_time_s
        path_length <- path_length + next_attachment$distance_um
        filament <- next_attachment$filament
        arc <- next_attachment$arc_um
        final_point <- next_attachment$point
        if (elapsed > deadline_s) break
      }
    }
    straight <- sqrt(sum((final_point - source) ^ 2))
    results[[cargo_id]] <- data.frame(cargo_id = cargo_id, success = success,
      delivery_time_s = if (success) elapsed else NA_real_, hops = hops,
      path_length_um = path_length, straight_distance_um = straight,
      route_inefficiency = if (success && straight > 1e-9) path_length / straight else NA_real_,
      junctions_encountered = junction_count,
      max_junction_multiplicity = max_junction_multiplicity,
      source_x = source[1], source_y = source[2], source_z = source[3],
      final_x = final_point[1], final_y = final_point[2], final_z = final_point[3])
  }
  output <- do.call(rbind, results)
  row.names(output) <- NULL
  output
}

summarize_mt_transport <- function(results) {
  delivered <- results[results$success, , drop = FALSE]
  data.frame(n_cargo = nrow(results), success_probability = mean(results$success),
    median_delivery_time_s = if (nrow(delivered)) stats::median(delivered$delivery_time_s) else NA,
    mean_delivery_time_s = if (nrow(delivered)) mean(delivered$delivery_time_s) else NA,
    median_hops = stats::median(results$hops), mean_hops = mean(results$hops),
    median_route_inefficiency = if (nrow(delivered))
      stats::median(delivered$route_inefficiency) else NA,
    mean_route_inefficiency = if (nrow(delivered)) mean(delivered$route_inefficiency) else NA)
}

.ellipsoid_surface_jacobian <- function(azimuth, equal_area_latitude, axes) {
  latitude <- asin(equal_area_latitude)
  cos_latitude <- cos(latitude)
  sin_latitude <- sin(latitude)
  r_azimuth <- c(-axes[1] * cos_latitude * sin(azimuth),
                  axes[2] * cos_latitude * cos(azimuth), 0)
  r_latitude <- c(-axes[1] * sin_latitude * cos(azimuth),
                  -axes[2] * sin_latitude * sin(azimuth),
                  axes[3] * cos_latitude)
  cross <- c(r_azimuth[2] * r_latitude[3] - r_azimuth[3] * r_latitude[2],
             r_azimuth[3] * r_latitude[1] - r_azimuth[1] * r_latitude[3],
             r_azimuth[1] * r_latitude[2] - r_azimuth[2] * r_latitude[1])
  sqrt(sum(cross ^ 2)) / max(cos_latitude, 1e-12)
}

.jensen_shannon_bits <- function(p, q) {
  p <- p / sum(p)
  q <- q / sum(q)
  midpoint <- (p + q) / 2
  kl <- function(a, b) sum(ifelse(a > 0, a * log2(a / b), 0))
  (kl(p, midpoint) + kl(q, midpoint)) / 2
}

#' Bin successful deliveries on an ellipsoidal membrane.
#'
#' Longitude and sin(latitude) form equal-solid-angle bins. Expected bin mass
#' is corrected by the ellipsoid surface Jacobian, so enrichment is relative to
#' uniform membrane area rather than uniform angle.
summarize_membrane_distribution <- function(cargo_results, geometry,
                                            n_azimuth = 12L,
                                            n_latitude = 4L,
                                            pseudocount = 0.5) {
  delivered <- cargo_results[cargo_results$success, , drop = FALSE]
  if (!nrow(delivered)) stop("No successful deliveries to summarize")
  points <- as.matrix(delivered[, c("final_x", "final_y", "final_z")])
  scaled <- sweep(points, 2, geometry$cell_axes_um, "/")
  directions <- scaled / sqrt(rowSums(scaled ^ 2))
  azimuth <- atan2(directions[, 2], directions[, 1])
  equal_area_latitude <- directions[, 3]
  azimuth_breaks <- seq(-pi, pi, length.out = n_azimuth + 1L)
  latitude_breaks <- seq(-1, 1, length.out = n_latitude + 1L)
  azimuth_bin <- pmax(1L, pmin(n_azimuth,
    findInterval(azimuth, azimuth_breaks, all.inside = TRUE)))
  latitude_bin <- pmax(1L, pmin(n_latitude,
    findInterval(equal_area_latitude, latitude_breaks, all.inside = TRUE)))
  linear_bin <- azimuth_bin + (latitude_bin - 1L) * n_azimuth
  counts <- tabulate(linear_bin, nbins = n_azimuth * n_latitude)

  grid <- expand.grid(azimuth_bin = seq_len(n_azimuth),
                      latitude_bin = seq_len(n_latitude))
  grid$azimuth_mid <- (azimuth_breaks[grid$azimuth_bin] +
                         azimuth_breaks[grid$azimuth_bin + 1L]) / 2
  grid$equal_area_latitude_mid <- (latitude_breaks[grid$latitude_bin] +
                                     latitude_breaks[grid$latitude_bin + 1L]) / 2
  grid$surface_weight <- mapply(.ellipsoid_surface_jacobian,
    grid$azimuth_mid, grid$equal_area_latitude_mid,
    MoreArgs = list(axes = geometry$cell_axes_um))
  grid$expected_probability <- grid$surface_weight / sum(grid$surface_weight)
  grid$count <- counts
  grid$observed_probability <- (counts + pseudocount) /
    (sum(counts) + pseudocount * length(counts))
  grid$enrichment <- grid$observed_probability / grid$expected_probability
  grid$log2_enrichment <- log2(grid$enrichment)

  observed <- grid$observed_probability
  expected <- grid$expected_probability
  metrics <- data.frame(n_delivered = nrow(delivered),
    total_variation_from_uniform = 0.5 * sum(abs(observed - expected)),
    js_from_uniform_bits = .jensen_shannon_bits(observed, expected),
    polarization_resultant = sqrt(sum(colMeans(directions) ^ 2)))
  list(bins = grid, metrics = metrics, directions = directions)
}

#' Compare two binned membrane-delivery distributions.
compare_membrane_distributions <- function(first, second) {
  stopifnot(nrow(first$bins) == nrow(second$bins))
  data.frame(js_between_sizes_bits = .jensen_shannon_bits(
    first$bins$observed_probability, second$bins$observed_probability),
    total_variation_between_sizes = 0.5 * sum(abs(
      first$bins$observed_probability - second$bins$observed_probability)))
}
