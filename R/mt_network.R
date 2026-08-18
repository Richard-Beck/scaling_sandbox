# Shariff-style generative model for 3D microtubule distributions.

.unit_vector <- function(x) {
  norm <- sqrt(sum(x ^ 2))
  if (!is.finite(norm) || norm < .Machine$double.eps) return(rep(0, length(x)))
  x / norm
}

#' Define a 3D ellipsoidal cell and excluded ellipsoidal nucleus.
mt_geometry <- function(cell_axes_um, nucleus_axes_um,
                        nucleus_center_um = c(0, 0, 0), mtoc_um = NULL) {
  stopifnot(length(cell_axes_um) == 3L, length(nucleus_axes_um) == 3L,
            length(nucleus_center_um) == 3L, all(cell_axes_um > 0),
            all(nucleus_axes_um > 0), all(nucleus_axes_um < cell_axes_um))
  if (is.null(mtoc_um)) {
    mtoc_um <- nucleus_center_um + c(nucleus_axes_um[1] + 0.35, 0, 0)
  }
  stopifnot(length(mtoc_um) == 3L)
  structure(list(dimension = 3L, cell_axes_um = as.numeric(cell_axes_um),
                 nucleus_axes_um = as.numeric(nucleus_axes_um),
                 nucleus_center_um = as.numeric(nucleus_center_um),
                 mtoc_um = as.numeric(mtoc_um)), class = "mt_geometry")
}

#' Isotropically enlarge a 3D geometry by a volume multiplier.
scale_mt_geometry <- function(geometry, volume_multiplier) {
  stopifnot(inherits(geometry, "mt_geometry"), volume_multiplier > 0)
  factor <- volume_multiplier ^ (1 / 3)
  mt_geometry(geometry$cell_axes_um * factor,
              geometry$nucleus_axes_um * factor,
              geometry$nucleus_center_um * factor,
              geometry$mtoc_um * factor)
}

#' Cytoplasmic volume of an MT geometry.
mt_cytoplasm_measure <- function(geometry) {
  4 * pi / 3 * (prod(geometry$cell_axes_um) - prod(geometry$nucleus_axes_um))
}

.ellipsoid_level <- function(point, axes, center = c(0, 0, 0)) {
  sum(((point - center) / axes) ^ 2)
}

.inside_cytoplasm <- function(point, geometry, tolerance = 1e-10) {
  .ellipsoid_level(point, geometry$cell_axes_um) <= 1 + tolerance &&
    .ellipsoid_level(point, geometry$nucleus_axes_um,
                     geometry$nucleus_center_um) >= 1 - tolerance
}

.random_unit_vector_3d <- function() .unit_vector(stats::rnorm(3L))

.sample_point_in_sphere <- function(center, radius_um) {
  center + .random_unit_vector_3d() * radius_um * stats::runif(1) ^ (1 / 3)
}

.sample_centrosomal_origin <- function(geometry, diameter_um) {
  for (attempt in seq_len(10000L)) {
    candidate <- .sample_point_in_sphere(geometry$mtoc_um, diameter_um / 2)
    if (.inside_cytoplasm(candidate, geometry)) return(candidate)
  }
  stop("Could not sample a cytosolic point in the centrosomal sphere")
}

.sample_cone_direction <- function(tangent, cos_alpha) {
  tangent <- .unit_vector(tangent)
  # Uniform solid-angle sampling within the spherical cap.
  axial <- stats::runif(1, cos_alpha, 1)
  azimuth <- stats::runif(1, 0, 2 * pi)
  reference <- diag(3)[, which.min(abs(tangent))]
  basis_1 <- .unit_vector(c(tangent[2] * reference[3] - tangent[3] * reference[2],
                            tangent[3] * reference[1] - tangent[1] * reference[3],
                            tangent[1] * reference[2] - tangent[2] * reference[1]))
  basis_2 <- c(tangent[2] * basis_1[3] - tangent[3] * basis_1[2],
               tangent[3] * basis_1[1] - tangent[1] * basis_1[3],
               tangent[1] * basis_1[2] - tangent[2] * basis_1[1])
  radial <- sqrt(max(0, 1 - axial ^ 2))
  .unit_vector(axial * tangent + radial *
                 (cos(azimuth) * basis_1 + sin(azimuth) * basis_2))
}

.sample_positive_normal <- function(mean_um, sd_um) {
  if (sd_um == 0) return(mean_um)
  repeat {
    value <- stats::rnorm(1, mean_um, sd_um)
    if (value > 0) return(value)
  }
}

.polyline_length <- function(points) {
  if (nrow(points) < 2L) return(0)
  sum(sqrt(rowSums((points[-1, , drop = FALSE] -
                      points[-nrow(points), , drop = FALSE]) ^ 2)))
}

.trim_polyline <- function(points, target_length_um) {
  if (target_length_um <= 0 || nrow(points) < 2L) return(points[1, , drop = FALSE])
  lengths <- sqrt(rowSums((points[-1, , drop = FALSE] -
                             points[-nrow(points), , drop = FALSE]) ^ 2))
  cumulative <- cumsum(lengths)
  if (target_length_um >= tail(cumulative, 1)) return(points)
  segment <- which(cumulative >= target_length_um)[1]
  before <- if (segment == 1L) 0 else cumulative[segment - 1L]
  fraction <- (target_length_um - before) / lengths[segment]
  endpoint <- points[segment, ] + fraction * (points[segment + 1L, ] - points[segment, ])
  rbind(points[seq_len(segment), , drop = FALSE], endpoint)
}

.grow_shariff_mt <- function(geometry, requested_length_um, step_um, cos_alpha,
                             centrosome_diameter_um, max_direction_attempts) {
  origin <- .sample_centrosomal_origin(geometry, centrosome_diameter_um)
  points <- matrix(origin, nrow = 1L)
  direction <- NULL
  realized <- 0
  terminated <- FALSE

  while (realized < requested_length_um - 1e-12) {
    step_length <- min(step_um, requested_length_um - realized)
    accepted <- FALSE
    for (attempt in seq_len(max_direction_attempts)) {
      candidate_direction <- if (is.null(direction)) {
        .random_unit_vector_3d()
      } else {
        .sample_cone_direction(direction, cos_alpha)
      }
      candidate <- points[nrow(points), ] + step_length * candidate_direction
      if (.inside_cytoplasm(candidate, geometry)) {
        direction <- candidate_direction
        points <- rbind(points, candidate)
        realized <- realized + step_length
        accepted <- TRUE
        break
      }
    }
    if (!accepted) {
      terminated <- TRUE
      break
    }
  }
  list(points = points, requested_length_um = requested_length_um,
       realized_length_um = .polyline_length(points),
       boundary_terminated = terminated)
}

.violates_global_turn_rule <- function(directions, window_steps,
                                       global_angle_deg, max_large_pairs) {
  if (nrow(directions) < 2L) return(FALSE)
  recent <- tail(directions, window_steps)
  dots <- pmax(-1, pmin(1, tcrossprod(recent)))
  angles <- acos(dots) * 180 / pi
  sum(angles[upper.tri(angles)] > global_angle_deg) > max_large_pairs
}

.grow_li_once <- function(geometry, requested_length_um, step_um, cos_alpha,
                          centrosome_diameter_um, max_direction_attempts,
                          local_turn_deg, global_turn_deg,
                          global_window_steps, max_global_large_pairs) {
  origin <- .sample_centrosomal_origin(geometry, centrosome_diameter_um)
  points <- matrix(origin, nrow = 1L)
  directions <- matrix(numeric(), nrow = 0L, ncol = 3L)
  direction <- NULL
  realized <- 0

  while (realized < requested_length_um - 1e-12) {
    step_length <- min(step_um, requested_length_um - realized)
    accepted_direction <- NULL

    # Ordinary Shariff-style persistent step.
    for (attempt in seq_len(max_direction_attempts)) {
      candidate_direction <- if (is.null(direction)) .random_unit_vector_3d() else
        .sample_cone_direction(direction, cos_alpha)
      candidate <- points[nrow(points), ] + step_length * candidate_direction
      if (.inside_cytoplasm(candidate, geometry)) {
        accepted_direction <- candidate_direction
        break
      }
    }

    # If blocked, relax collinearity but retain Li et al.'s local turn limit.
    if (is.null(accepted_direction) && !is.null(direction)) {
      local_cos <- cos(local_turn_deg * pi / 180)
      for (attempt in seq_len(max_direction_attempts)) {
        candidate_direction <- .sample_cone_direction(direction, local_cos)
        candidate <- points[nrow(points), ] + step_length * candidate_direction
        if (.inside_cytoplasm(candidate, geometry)) {
          accepted_direction <- candidate_direction
          break
        }
      }
    }
    if (is.null(accepted_direction)) return(NULL)

    proposed_directions <- rbind(directions, accepted_direction)
    if (.violates_global_turn_rule(proposed_directions, global_window_steps,
                                   global_turn_deg, max_global_large_pairs)) {
      return(NULL)
    }
    direction <- accepted_direction
    directions <- proposed_directions
    points <- rbind(points, points[nrow(points), ] + step_length * direction)
    realized <- realized + step_length
  }
  list(points = points, requested_length_um = requested_length_um,
       realized_length_um = .polyline_length(points),
       boundary_terminated = FALSE)
}

.grow_li_mt <- function(geometry, requested_length_um, step_um, cos_alpha,
                        centrosome_diameter_um, max_direction_attempts,
                        trajectory_retries, local_turn_deg, global_turn_deg,
                        global_window_steps, max_global_large_pairs) {
  for (attempt in seq_len(trajectory_retries)) {
    mt <- .grow_li_once(geometry, requested_length_um, step_um, cos_alpha,
      centrosome_diameter_um, max_direction_attempts, local_turn_deg,
      global_turn_deg, global_window_steps, max_global_large_pairs)
    if (!is.null(mt)) return(mt)
  }
  stop("Li-style generator could not realize an MT of ",
       signif(requested_length_um, 4), " um after ", trajectory_retries,
       " complete trajectory attempts")
}

#' Generate a 3D Shariff-style persistent MT network.
#'
#' Specify either `n_filaments` or `target_density`. In density-target mode the
#' last path is trimmed so independently generated cell sizes match exactly.
generate_mt_network <- function(geometry, n_filaments = NULL,
                                target_density = NULL, step_um = 0.2,
                                mean_length_um = 25, sd_length_um = 15,
                                cos_alpha = 0.9,
                                centrosome_diameter_um = 0.4,
                                max_direction_attempts = 250L,
                                boundary_mode = c("shariff", "li"),
                                trajectory_retries = 100L,
                                local_turn_deg = 63.9,
                                global_turn_deg = 120,
                                global_window_steps = 30L,
                                max_global_large_pairs = 3L,
                                max_filaments = 10000L, seed = NULL) {
  boundary_mode <- match.arg(boundary_mode)
  stopifnot(inherits(geometry, "mt_geometry"), xor(is.null(n_filaments),
                                                    is.null(target_density)),
            step_um > 0, mean_length_um > 0, sd_length_um >= 0,
            cos_alpha >= -1, cos_alpha <= 1,
            centrosome_diameter_um > 0, max_direction_attempts > 0)
  if (!is.null(seed)) set.seed(seed)
  target_length <- if (is.null(target_density)) Inf else
    target_density * mt_cytoplasm_measure(geometry)
  filaments <- list()
  total_length <- 0

  repeat {
    if (!is.null(n_filaments) && length(filaments) >= n_filaments) break
    if (is.null(n_filaments) && total_length >= target_length - 1e-10) break
    if (length(filaments) >= max_filaments) stop("max_filaments reached")
    requested <- .sample_positive_normal(mean_length_um, sd_length_um)
    mt <- if (boundary_mode == "shariff") {
      .grow_shariff_mt(geometry, requested, step_um, cos_alpha,
                       centrosome_diameter_um, max_direction_attempts)
    } else {
      .grow_li_mt(geometry, requested, step_um, cos_alpha,
        centrosome_diameter_um, max_direction_attempts, trajectory_retries,
        local_turn_deg, global_turn_deg, global_window_steps,
        max_global_large_pairs)
    }
    if (mt$realized_length_um <= 1e-10) next
    if (is.null(n_filaments)) {
      remaining <- target_length - total_length
      if (mt$realized_length_um > remaining) {
        mt$points <- .trim_polyline(mt$points, remaining)
        mt$realized_length_um <- .polyline_length(mt$points)
        mt$requested_length_um <- min(mt$requested_length_um, remaining)
      }
    }
    filaments[[length(filaments) + 1L]] <- mt
    total_length <- total_length + mt$realized_length_um
  }

  achieved_density <- total_length / mt_cytoplasm_measure(geometry)
  structure(list(geometry = geometry, filaments = filaments,
                 target_density = if (is.null(target_density)) achieved_density else target_density,
                 total_length_um = total_length,
                 parameters = list(n_filaments = n_filaments, step_um = step_um,
                   mean_length_um = mean_length_um, sd_length_um = sd_length_um,
                   cos_alpha = cos_alpha,
                   centrosome_diameter_um = centrosome_diameter_um,
                   max_direction_attempts = max_direction_attempts,
                   boundary_mode = boundary_mode,
                   trajectory_retries = trajectory_retries,
                   local_turn_deg = local_turn_deg,
                   global_turn_deg = global_turn_deg,
                   global_window_steps = global_window_steps,
                   max_global_large_pairs = max_global_large_pairs,
                   seed = seed)),
            class = "mt_network")
}

#' Generate 1x and 2x networks with matched realized polymer density.
generate_matched_mt_networks <- function(base_geometry, volume_multiplier = 2,
                                         n_base = 175L, seed = 1L, ...) {
  baseline <- generate_mt_network(base_geometry, n_filaments = n_base,
                                  seed = seed, ...)
  enlarged <- generate_mt_network(scale_mt_geometry(base_geometry, volume_multiplier),
                                  target_density = baseline$target_density,
                                  seed = seed + 1L, ...)
  list(`1x volume` = baseline, `2x volume` = enlarged)
}

#' Convert paths to a segment table for projections and spatial analysis.
mt_segments <- function(network) {
  pieces <- lapply(seq_along(network$filaments), function(i) {
    mt <- network$filaments[[i]]
    points <- mt$points
    if (nrow(points) < 2L) return(NULL)
    starts <- points[-nrow(points), , drop = FALSE]
    ends <- points[-1, , drop = FALSE]
    data.frame(filament = i, segment = seq_len(nrow(points) - 1L),
      x = starts[, 1], y = starts[, 2], z = starts[, 3],
      xend = ends[, 1], yend = ends[, 2], zend = ends[, 3],
      z_mid = (starts[, 3] + ends[, 3]) / 2,
      length_um = sqrt(rowSums((ends - starts) ^ 2)))
  })
  do.call(rbind, Filter(Negate(is.null), pieces))
}

.smooth_1d <- function(values, kernel) {
  radius <- (length(kernel) - 1L) %/% 2L
  padded <- c(rep(0, radius), values, rep(0, radius))
  vapply(seq_along(values), function(i) {
    sum(padded[i + seq_along(kernel) - 1L] * kernel)
  }, numeric(1))
}

.smooth_density_grid <- function(grid, spacing_um, kernel_um) {
  kernels <- lapply(spacing_um, function(spacing) {
    radius <- max(1L, ceiling(3 * kernel_um / spacing))
    offsets <- seq.int(-radius, radius)
    weights <- exp(-(offsets * spacing) ^ 2 / (2 * kernel_um ^ 2))
    weights / sum(weights)
  })
  out <- grid
  for (j in seq_len(dim(grid)[2])) for (k in seq_len(dim(grid)[3])) {
    out[, j, k] <- .smooth_1d(grid[, j, k], kernels[[1]])
  }
  grid <- out
  for (i in seq_len(dim(grid)[1])) for (k in seq_len(dim(grid)[3])) {
    out[i, , k] <- .smooth_1d(grid[i, , k], kernels[[2]])
  }
  grid <- out
  for (i in seq_len(dim(grid)[1])) for (j in seq_len(dim(grid)[2])) {
    out[i, j, ] <- .smooth_1d(grid[i, j, ], kernels[[3]])
  }
  out
}

.mt_density_field <- function(network, grid_n, kernel_um) {
  stopifnot(length(grid_n) == 3L, all(grid_n >= 3), kernel_um > 0)
  axes <- network$geometry$cell_axes_um
  breaks <- lapply(seq_len(3L), function(i) seq(-axes[i], axes[i],
                                                length.out = grid_n[i] + 1L))
  spacing <- vapply(breaks, function(x) x[2] - x[1], numeric(1))
  segments <- mt_segments(network)
  midpoints <- cbind((segments$x + segments$xend) / 2,
                     (segments$y + segments$yend) / 2,
                     (segments$z + segments$zend) / 2)
  indices <- sapply(seq_len(3L), function(i) {
    pmax(1L, pmin(grid_n[i], findInterval(midpoints[, i], breaks[[i]],
                                          all.inside = TRUE)))
  })
  linear <- indices[, 1] + (indices[, 2] - 1L) * grid_n[1] +
    (indices[, 3] - 1L) * grid_n[1] * grid_n[2]
  totals <- tapply(segments$length_um, linear, sum)
  density <- array(0, dim = grid_n)
  density[as.integer(names(totals))] <- as.numeric(totals) / prod(spacing)
  list(values = .smooth_density_grid(density, spacing, kernel_um),
       breaks = breaks, grid_n = grid_n)
}

.sample_density_field <- function(field, points) {
  indices <- sapply(seq_len(3L), function(i) {
    pmax(1L, pmin(field$grid_n[i], findInterval(points[, i], field$breaks[[i]],
                                                all.inside = TRUE)))
  })
  if (is.null(dim(indices))) indices <- matrix(indices, nrow = 1L)
  field$values[cbind(indices[, 1], indices[, 2], indices[, 3])]
}

.score_mt_density <- function(network, field, score_mode) {
  vapply(network$filaments, function(mt) {
    points <- mt$points
    if (score_mode == "distal") {
      return(.sample_density_field(field, points[nrow(points), , drop = FALSE]))
    }
    starts <- points[-nrow(points), , drop = FALSE]
    ends <- points[-1, , drop = FALSE]
    lengths <- sqrt(rowSums((ends - starts) ^ 2))
    densities <- .sample_density_field(field, (starts + ends) / 2)
    stats::weighted.mean(densities, lengths)
  }, numeric(1))
}

.removal_probabilities <- function(scores, removal_fraction, score_power) {
  if (!length(scores)) return(numeric())
  if (all(scores <= 0) || removal_fraction <= 0) return(rep(removal_fraction, length(scores)))
  weights <- (scores / mean(scores)) ^ score_power
  objective <- function(scale) mean(pmin(1, scale * weights)) - removal_fraction
  upper <- 1
  while (objective(upper) < 0) upper <- upper * 2
  scale <- stats::uniroot(objective, c(0, upper))$root
  pmin(1, scale * weights)
}

.refill_mt_network <- function(network, seed) {
  measure <- mt_cytoplasm_measure(network$geometry)
  target_length <- network$target_density * measure
  current_length <- sum(vapply(network$filaments, `[[`, numeric(1),
                               "realized_length_um"))
  remaining <- target_length - current_length
  if (remaining <= 1e-10) return(network)
  args <- network$parameters
  args$n_filaments <- NULL
  args$seed <- seed
  refill <- do.call(generate_mt_network, c(list(geometry = network$geometry,
    target_density = remaining / measure), args))
  network$filaments <- c(network$filaments, refill$filaments)
  network$total_length_um <- current_length + refill$total_length_um
  network
}

#' Iteratively prune MTs from locally dense regions.
#'
#' `score_mode = "integrated"` uses the length-weighted mean density along each
#' path; `"distal"` samples only its plus end. By default removed polymer is
#' regenerated after every pass to retain the original density target.
prune_mt_network <- function(network, iterations = 10L,
                             score_mode = c("integrated", "distal"),
                             removal_fraction = 0.5,
                             density_kernel_um = 1.5,
                             grid_n = c(24L, 16L, 8L), score_power = 1,
                             refill_to_target = TRUE, seed = 1L) {
  stopifnot(inherits(network, "mt_network"), iterations >= 0,
            removal_fraction >= 0, removal_fraction <= 0.5,
            density_kernel_um > 0, score_power > 0)
  score_mode <- match.arg(score_mode)
  set.seed(seed)
  result <- unserialize(serialize(network, NULL))
  history <- data.frame()
  for (iteration in seq_len(iterations)) {
    before_count <- length(result$filaments)
    field <- .mt_density_field(result, as.integer(grid_n), density_kernel_um)
    scores <- .score_mt_density(result, field, score_mode)
    weights <- if (all(scores <= 0)) rep(1, length(scores)) else
      (scores + .Machine$double.eps) ^ score_power
    remove_n <- min(length(scores) - 1L, floor(removal_fraction * length(scores)))
    remove <- rep(FALSE, length(scores))
    if (remove_n > 0) {
      remove[sample(seq_along(scores), size = remove_n,
                    replace = FALSE, prob = weights)] <- TRUE
    }
    result$filaments <- result$filaments[!remove]
    result$total_length_um <- sum(vapply(result$filaments, `[[`, numeric(1),
                                         "realized_length_um"))
    density_after_prune <- result$total_length_um / mt_cytoplasm_measure(result$geometry)
    if (refill_to_target) result <- .refill_mt_network(result, seed + iteration)
    history <- rbind(history, data.frame(iteration = iteration,
      score_mode = score_mode, before_count = before_count,
      removed_count = sum(remove), realized_removal_fraction = mean(remove),
      mean_score = mean(scores), density_after_prune = density_after_prune,
      after_refill_count = length(result$filaments),
      final_density = result$total_length_um / mt_cytoplasm_measure(result$geometry)))
  }
  result$pruning <- list(parameters = list(iterations = iterations,
    score_mode = score_mode, removal_fraction = removal_fraction,
    density_kernel_um = density_kernel_um, grid_n = grid_n,
    score_power = score_power, refill_to_target = refill_to_target, seed = seed),
    history = history)
  result
}

mt_path_points <- function(network, max_filaments = 100L, seed = 1L) {
  ids <- seq_along(network$filaments)
  if (length(ids) > max_filaments) {
    set.seed(seed)
    ids <- sort(sample(ids, max_filaments))
  }
  do.call(rbind, lapply(ids, function(i) {
    points <- network$filaments[[i]]$points
    data.frame(filament = i, point = seq_len(nrow(points)),
               x = points[, 1], y = points[, 2], z = points[, 3])
  }))
}

#' Summarize calibration and path completion.
summarize_mt_network <- function(network) {
  lengths <- vapply(network$filaments, `[[`, numeric(1), "realized_length_um")
  requested <- vapply(network$filaments, `[[`, numeric(1), "requested_length_um")
  terminated <- vapply(network$filaments, `[[`, logical(1), "boundary_terminated")
  measure <- mt_cytoplasm_measure(network$geometry)
  data.frame(dimension = "3D", cytoplasm_volume_um3 = measure,
    filament_count = length(lengths), total_length_um = sum(lengths),
    achieved_density = sum(lengths) / measure,
    mean_realized_length_um = mean(lengths),
    median_realized_length_um = stats::median(lengths),
    mean_requested_length_um = mean(requested),
    boundary_terminated_fraction = mean(terminated))
}

mt_boundary_data <- function(geometry, projection = "xy", n = 240L) {
  indices <- match(strsplit(projection, "")[[1]], c("x", "y", "z"))
  if (anyNA(indices) || length(indices) != 2L) stop("projection must be xy, xz, or yz")
  theta <- seq(0, 2 * pi, length.out = n)
  make_boundary <- function(axes, center, boundary) {
    data.frame(horizontal = center[indices[1]] + axes[indices[1]] * cos(theta),
      vertical = center[indices[2]] + axes[indices[2]] * sin(theta), boundary = boundary)
  }
  rbind(make_boundary(geometry$cell_axes_um, c(0, 0, 0), "Cell"),
        make_boundary(geometry$nucleus_axes_um, geometry$nucleus_center_um, "Nucleus"))
}

#' Plot an orthogonal projection of a 3D MT network.
plot_mt_network <- function(network, projection = "xy", title = NULL,
                            max_filaments = 100L, seed = 1L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
  axes <- strsplit(projection, "")[[1]]
  paths <- mt_path_points(network, max_filaments, seed)
  boundary <- mt_boundary_data(network$geometry, projection)
  plot <- ggplot2::ggplot() +
    ggplot2::geom_path(data = boundary,
      ggplot2::aes(horizontal, vertical, group = boundary, linetype = boundary),
      colour = "grey30", linewidth = 0.5) +
    ggplot2::geom_path(data = paths,
      ggplot2::aes(x = .data[[axes[1]]], y = .data[[axes[2]]], group = filament),
      colour = "#126782", linewidth = 0.28, alpha = 0.65) +
    ggplot2::coord_equal() +
    ggplot2::scale_linetype_manual(values = c(Cell = "solid", Nucleus = "dashed")) +
    ggplot2::labs(title = title, x = paste0(axes[1], " (um)"),
      y = paste0(axes[2], " (um)"), linetype = "Boundary") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
  plot
}
