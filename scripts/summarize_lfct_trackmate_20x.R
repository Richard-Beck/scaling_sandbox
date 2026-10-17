#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(xml2)
  library(ggplot2)
  library(scales)
})

input_dir <- "data/lfct_trackmate_20x"
output_dir <- "reports/lfct_20x_track_summary_assets"
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
shard_dir <- file.path(output_dir, "shards")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(shard_dir, recursive = TRUE, showWarnings = FALSE)

cell_labels <- c(
  mcf10a = "MCF10A",
  hek293 = "HEK293",
  mdamb = "MDA-MB-231",
  u87 = "U87"
)

num_attr <- function(nodes, name) as.numeric(xml_attr(nodes, name))

roi_bbox <- function(text, x, y) {
  vals <- suppressWarnings(as.numeric(strsplit(trimws(text), "[[:space:]]+")[[1]]))
  if (length(vals) < 4L || length(vals) %% 2L != 0L || anyNA(vals)) {
    return(c(xmin = x, xmax = x, ymin = y, ymax = y))
  }
  xx <- vals[seq.int(1L, length(vals), 2L)] + x
  yy <- vals[seq.int(2L, length(vals), 2L)] + y
  c(xmin = min(xx), xmax = max(xx), ymin = min(yy), ymax = max(yy))
}

add_neighbor_metrics <- function(spots) {
  spots[, c("nn_center", "nn_box_gap") := {
    n <- .N
    if (n < 2L) {
      list(rep(NA_real_, n), rep(NA_real_, n))
    } else {
      center_dist <- sqrt(outer(x, x, "-")^2 + outer(y, y, "-")^2)
      diag(center_dist) <- Inf
      center_nn <- apply(center_dist, 1L, min)
      box_nn <- numeric(n)
      for (i in seq_len(n)) {
        gx <- pmax(0, xmin - xmax[i], xmin[i] - xmax)
        gy <- pmax(0, ymin - ymax[i], ymin[i] - ymax)
        gap <- sqrt(gx^2 + gy^2)
        gap[i] <- Inf
        box_nn[i] <- min(gap)
      }
      list(center_nn, box_nn)
    }
  }, by = frame]
  spots[, isolated_0 := nn_box_gap > 0]
  spots[, isolated_5 := nn_box_gap > 5]
  spots[, isolated_10 := nn_box_gap > 10]
  spots[]
}

rolling_mean <- function(x, n = 21L) {
  if (length(x) < n) return(rep(mean(x, na.rm = TRUE), length(x)))
  frollmean(x, n = n, align = "center", na.rm = TRUE)
}

all_spots <- list()
all_edges <- list()
all_tracks <- list()
all_corr <- list()
all_msd <- list()
field_meta <- list()

xml_files <- sort(list.files(input_dir, pattern = "[.]xml$", recursive = TRUE, full.names = TRUE))
if (!length(xml_files)) stop("No TrackMate XML files found in ", input_dir)

aggregate_mode <- identical(Sys.getenv("LFCT_AGGREGATE", unset = "0"), "1")
worker_id_text <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
worker_id <- if (nzchar(worker_id_text)) as.integer(worker_id_text) else NA_integer_

if (aggregate_mode) {
  shard_files <- sort(list.files(shard_dir, pattern = "[.]rds$", full.names = TRUE))
  if (length(shard_files) != length(xml_files)) {
    stop("Expected ", length(xml_files), " completed shards, found ", length(shard_files))
  }
  shards <- lapply(shard_files, readRDS)
  all_spots <- lapply(shards, `[[`, "spots")
  all_edges <- lapply(shards, `[[`, "edges")
  all_tracks <- lapply(shards, `[[`, "tracks")
  all_corr <- lapply(shards, `[[`, "corr")
  all_msd <- lapply(shards, `[[`, "msd")
  field_meta <- lapply(shards, `[[`, "field_meta")
} else {
if (!is.na(worker_id)) {
  if (worker_id < 1L || worker_id > length(xml_files)) stop("Invalid array task index: ", worker_id)
  xml_files <- xml_files[worker_id]
}

for (file_index in seq_along(xml_files)) {
  path <- xml_files[file_index]
  cell_key <- basename(dirname(path))
  cell_line <- unname(cell_labels[cell_key])
  field <- tools::file_path_sans_ext(basename(path))
  message(sprintf("[%d/%d] %s %s", file_index, length(xml_files), cell_line, field))

  doc <- read_xml(path, options = "HUGE")
  image_data <- xml_find_first(doc, ".//Settings/ImageData")
  spot_nodes <- xml_find_all(doc, ".//Model/AllSpots/SpotsInFrame/Spot")
  spots <- data.table(
    id = num_attr(spot_nodes, "ID"),
    frame = num_attr(spot_nodes, "FRAME"),
    x = num_attr(spot_nodes, "POSITION_X"),
    y = num_attr(spot_nodes, "POSITION_Y"),
    area = num_attr(spot_nodes, "AREA"),
    radius = num_attr(spot_nodes, "RADIUS"),
    circularity = num_attr(spot_nodes, "CIRCULARITY"),
    solidity = num_attr(spot_nodes, "SOLIDITY")
  )
  bbox <- t(vapply(
    seq_along(spot_nodes),
    function(i) roi_bbox(xml_text(spot_nodes[[i]]), spots$x[i], spots$y[i]),
    numeric(4L)
  ))
  spots[, `:=`(
    xmin = bbox[, "xmin"], xmax = bbox[, "xmax"],
    ymin = bbox[, "ymin"], ymax = bbox[, "ymax"]
  )]
  add_neighbor_metrics(spots)
  spots[, `:=`(cell_line = cell_line, field = field)]

  track_nodes <- xml_find_all(doc, ".//Model/AllTracks/Track")
  tracks <- rbindlist(lapply(track_nodes, function(track_node) {
    data.table(
      track_id = as.numeric(xml_attr(track_node, "TRACK_ID")),
      duration = as.numeric(xml_attr(track_node, "TRACK_DURATION")),
      start = as.numeric(xml_attr(track_node, "TRACK_START")),
      stop = as.numeric(xml_attr(track_node, "TRACK_STOP")),
      n_spots = as.numeric(xml_attr(track_node, "NUMBER_SPOTS")),
      n_splits = as.numeric(xml_attr(track_node, "NUMBER_SPLITS")),
      total_distance = as.numeric(xml_attr(track_node, "TOTAL_DISTANCE_TRAVELED")),
      net_displacement = as.numeric(xml_attr(track_node, "TRACK_DISPLACEMENT")),
      confinement_ratio = as.numeric(xml_attr(track_node, "CONFINEMENT_RATIO")),
      mean_speed = as.numeric(xml_attr(track_node, "TRACK_MEAN_SPEED"))
    )
  }))
  tracks[, `:=`(cell_line = cell_line, field = field)]

  edge_list <- lapply(track_nodes, function(track_node) {
    edge_nodes <- xml_find_all(track_node, "./Edge")
    if (!length(edge_nodes)) return(NULL)
    data.table(
      track_id = as.numeric(xml_attr(track_node, "TRACK_ID")),
      source = num_attr(edge_nodes, "SPOT_SOURCE_ID"),
      target = num_attr(edge_nodes, "SPOT_TARGET_ID"),
      xml_speed = num_attr(edge_nodes, "SPEED")
    )
  })
  edges <- rbindlist(edge_list)
  lookup <- spots[, .(id, frame, x, y, isolated_5)]
  edges <- merge(edges, lookup, by.x = "source", by.y = "id", all.x = TRUE)
  setnames(edges, c("frame", "x", "y", "isolated_5"), c("source_frame", "source_x", "source_y", "source_isolated_5"))
  edges <- merge(edges, lookup, by.x = "target", by.y = "id", all.x = TRUE)
  setnames(edges, c("frame", "x", "y", "isolated_5"), c("target_frame", "target_x", "target_y", "target_isolated_5"))
  edges[, `:=`(
    dx = target_x - source_x,
    dy = target_y - source_y,
    delta_frame = target_frame - source_frame
  )]
  edges[, speed := sqrt(dx^2 + dy^2) / delta_frame]
  edges[, isolated_step := source_isolated_5 & target_isolated_5]
  edges[, `:=`(cell_line = cell_line, field = field)]

  # Build unique continuations. Branches and merges are not crossed when calculating
  # direction autocorrelation or MSD, so division geometry does not masquerade as motion.
  edge_core <- edges[, .(source, target, source_frame, target_frame, dx, dy)]
  out_count <- edge_core[, .N, by = source]
  in_count <- edge_core[, .N, by = target]
  unique_out <- edge_core[out_count[N == 1L], on = "source"]
  next_edge <- setNames(seq_len(nrow(unique_out)), as.character(unique_out$source))
  indegree <- setNames(in_count$N, as.character(in_count$target))
  pos_x <- setNames(spots$x, as.character(spots$id))
  pos_y <- setNames(spots$y, as.character(spots$id))
  pos_frame <- setNames(spots$frame, as.character(spots$id))
  corr_rows <- vector("list", nrow(edge_core) * 3L)
  msd_rows <- vector("list", nrow(edge_core) * 5L)
  corr_n <- 0L
  msd_n <- 0L
  corr_lags <- c(1L, 3L, 5L)
  msd_lags <- c(1L, 3L, 5L, 10L, 20L)
  for (i in seq_len(nrow(edge_core))) {
    start_id <- edge_core$source[i]
    current_id <- start_id
    v0x <- edge_core$dx[i]
    v0y <- edge_core$dy[i]
    v0norm <- sqrt(v0x^2 + v0y^2)
    for (step in seq_len(max(msd_lags) + 1L)) {
      edge_idx <- unname(next_edge[as.character(current_id)])
      if (is.na(edge_idx)) break
      next_id <- unique_out$target[edge_idx]
      next_indegree <- unname(indegree[as.character(next_id)])
      if (!is.na(next_indegree) && next_indegree != 1L) break
      if (step %in% msd_lags) {
        msd_n <- msd_n + 1L
        ddx <- pos_x[[as.character(next_id)]] - pos_x[[as.character(start_id)]]
        ddy <- pos_y[[as.character(next_id)]] - pos_y[[as.character(start_id)]]
        msd_rows[[msd_n]] <- data.table(frame = pos_frame[[as.character(start_id)]], lag = step, squared_displacement = ddx^2 + ddy^2)
      }
      # The edge at source frame t+lag is compared with the edge at t.
      if ((step - 1L) %in% corr_lags && step > 1L && v0norm > 0) {
        vnorm <- sqrt(unique_out$dx[edge_idx]^2 + unique_out$dy[edge_idx]^2)
        if (vnorm > 0) {
          corr_n <- corr_n + 1L
          corr_rows[[corr_n]] <- data.table(
            frame = pos_frame[[as.character(start_id)]],
            lag = step - 1L,
            directional_correlation = (v0x * unique_out$dx[edge_idx] + v0y * unique_out$dy[edge_idx]) / (v0norm * vnorm)
          )
        }
      }
      current_id <- next_id
    }
  }
  corr <- rbindlist(corr_rows[seq_len(corr_n)])
  msd <- rbindlist(msd_rows[seq_len(msd_n)])
  if (nrow(corr)) corr[, `:=`(cell_line = cell_line, field = field)]
  if (nrow(msd)) msd[, `:=`(cell_line = cell_line, field = field)]

  field_meta[[length(field_meta) + 1L]] <- data.table(
    cell_line = cell_line,
    field = field,
    n_frames = as.numeric(xml_attr(image_data, "nframes")),
    n_detections = nrow(spots),
    n_tracks = nrow(tracks),
    n_divisions = sum(tracks$n_splits, na.rm = TRUE)
  )
  all_spots[[length(all_spots) + 1L]] <- spots[, .(cell_line, field, id, frame, x, y, area, radius, circularity, solidity, nn_center, nn_box_gap, isolated_0, isolated_5, isolated_10)]
  all_edges[[length(all_edges) + 1L]] <- edges[, .(cell_line, field, track_id, source, target, source_frame, target_frame, delta_frame, speed, isolated_step)]
  all_tracks[[length(all_tracks) + 1L]] <- tracks
  all_corr[[length(all_corr) + 1L]] <- corr
  all_msd[[length(all_msd) + 1L]] <- msd
  rm(doc, spot_nodes, track_nodes, spots, tracks, edges, corr, msd)
  gc(FALSE)
}

if (!is.na(worker_id)) {
  shard_path <- file.path(shard_dir, sprintf("shard_%03d.rds", worker_id))
  saveRDS(list(
    spots = all_spots[[1L]], edges = all_edges[[1L]], tracks = all_tracks[[1L]],
    corr = all_corr[[1L]], msd = all_msd[[1L]], field_meta = field_meta[[1L]]
  ), shard_path, compress = "gzip")
  message("Wrote shard ", shard_path)
  quit(save = "no", status = 0L)
}
}

spots <- rbindlist(all_spots)
edges <- rbindlist(all_edges)
tracks <- rbindlist(all_tracks)
corr <- rbindlist(all_corr)
msd <- rbindlist(all_msd)
fields <- rbindlist(field_meta)

field_summary <- merge(
  fields,
  tracks[, .(
    duration_q10 = quantile(duration, 0.10, na.rm = TRUE),
    duration_q25 = quantile(duration, 0.25, na.rm = TRUE),
    duration_median = median(duration, na.rm = TRUE),
    duration_mean = mean(duration, na.rm = TRUE),
    duration_q75 = quantile(duration, 0.75, na.rm = TRUE),
    duration_q90 = quantile(duration, 0.90, na.rm = TRUE),
    track_speed_median = median(mean_speed, na.rm = TRUE),
    confinement_median = median(confinement_ratio, na.rm = TRUE)
  ), by = .(cell_line, field)],
  by = c("cell_line", "field")
)
field_summary <- merge(
  field_summary,
  edges[delta_frame == 1 & is.finite(speed), .(
    speed_q10 = quantile(speed, 0.10),
    speed_q25 = quantile(speed, 0.25),
    speed_median = median(speed),
    speed_mean = mean(speed),
    speed_q75 = quantile(speed, 0.75),
    speed_q90 = quantile(speed, 0.90)
  ), by = .(cell_line, field)],
  by = c("cell_line", "field")
)
field_summary <- merge(
  field_summary,
  spots[, .(
    nn_center_median = median(nn_center, na.rm = TRUE),
    nn_box_gap_median = median(nn_box_gap, na.rm = TRUE),
    isolated_5_fraction = mean(isolated_5, na.rm = TRUE)
  ), by = .(cell_line, field)],
  by = c("cell_line", "field")
)

cell_summary <- field_summary[, .(
  n_fields = .N,
  n_frames_min = min(n_frames),
  n_frames_max = max(n_frames),
  n_detections = sum(n_detections),
  n_tracks = sum(n_tracks),
  n_divisions = sum(n_divisions),
  median_field_duration = median(duration_median),
  median_field_speed = median(track_speed_median),
  median_field_nn = median(nn_center_median),
  median_isolated_5_fraction = median(isolated_5_fraction)
), by = cell_line]

speed_summary <- edges[delta_frame == 1 & is.finite(speed), .(
  speed_q10 = quantile(speed, 0.10), speed_q25 = quantile(speed, 0.25),
  speed_median = median(speed), speed_mean = mean(speed),
  speed_q75 = quantile(speed, 0.75), speed_q90 = quantile(speed, 0.90)
), by = cell_line]
duration_summary <- tracks[is.finite(duration), .(
  duration_q10 = quantile(duration, 0.10), duration_q25 = quantile(duration, 0.25),
  duration_median = median(duration), duration_mean = mean(duration),
  duration_q75 = quantile(duration, 0.75), duration_q90 = quantile(duration, 0.90)
), by = cell_line]
neighbor_summary <- spots[is.finite(nn_center), .(
  nn_q10 = quantile(nn_center, 0.10), nn_q25 = quantile(nn_center, 0.25),
  nn_median = median(nn_center), nn_mean = mean(nn_center),
  nn_q75 = quantile(nn_center, 0.75), nn_q90 = quantile(nn_center, 0.90),
  isolated_5_fraction = mean(isolated_5, na.rm = TRUE)
), by = cell_line]
corr_summary <- dcast(
  corr[, .(directional_correlation = mean(directional_correlation, na.rm = TRUE)), by = .(cell_line, lag)],
  cell_line ~ lag, value.var = "directional_correlation"
)
setnames(corr_summary, intersect(c("1", "3", "5"), names(corr_summary)), paste0("direction_corr_lag", intersect(c("1", "3", "5"), names(corr_summary))))
msd_cell <- msd[, .(msd = mean(squared_displacement)), by = .(cell_line, lag)]
msd_wide <- dcast(msd_cell[lag %in% c(1, 5)], cell_line ~ lag, value.var = "msd")
msd_wide[, msd_alpha_1_5 := log(`5` / `1`) / log(5)]
msd_wide[, c("1", "5") := NULL]
cell_metric_summary <- Reduce(
  function(x, y) merge(x, y, by = "cell_line"),
  list(speed_summary, duration_summary, neighbor_summary, corr_summary, msd_wide)
)

fwrite(field_summary, file.path(table_dir, "field_summary.csv"))
fwrite(cell_summary, file.path(table_dir, "cell_line_summary.csv"))
fwrite(cell_metric_summary, file.path(table_dir, "cell_line_metric_summary.csv"))

speed_time <- edges[delta_frame == 1 & is.finite(speed), .(
  q10 = quantile(speed, 0.10), median = median(speed), q90 = quantile(speed, 0.90),
  isolated_median = if (any(isolated_step, na.rm = TRUE)) median(speed[isolated_step %in% TRUE]) else NA_real_,
  n = .N
), by = .(cell_line, field, frame = source_frame)]
setorder(speed_time, cell_line, field, frame)
speed_time[, `:=`(
  median_roll = rolling_mean(median),
  isolated_roll = rolling_mean(isolated_median)
), by = .(cell_line, field)]
speed_cell_time <- edges[delta_frame == 1 & is.finite(speed), .(
  q10 = quantile(speed, 0.10), median = median(speed), q90 = quantile(speed, 0.90)
), by = .(cell_line, frame = source_frame)]
setorder(speed_cell_time, cell_line, frame)
speed_cell_time[, `:=`(
  q10_roll = rolling_mean(q10), median_roll = rolling_mean(median), q90_roll = rolling_mean(q90)
), by = cell_line]

neighbor_time <- spots[is.finite(nn_center), .(
  center_median = median(nn_center),
  box_gap_median = median(nn_box_gap),
  isolated_0 = mean(isolated_0),
  isolated_5 = mean(isolated_5),
  isolated_10 = mean(isolated_10)
), by = .(cell_line, field, frame)]
setorder(neighbor_time, cell_line, field, frame)
neighbor_time[, `:=`(
  center_roll = rolling_mean(center_median),
  box_gap_roll = rolling_mean(box_gap_median),
  isolated_0_roll = rolling_mean(isolated_0),
  isolated_5_roll = rolling_mean(isolated_5),
  isolated_10_roll = rolling_mean(isolated_10)
), by = .(cell_line, field)]

corr_time <- corr[, .(correlation = mean(directional_correlation, na.rm = TRUE), n = .N), by = .(cell_line, field, frame, lag)]
setorder(corr_time, cell_line, field, lag, frame)
corr_time[, correlation_roll := rolling_mean(correlation), by = .(cell_line, field, lag)]

msd_curve <- msd[, .(msd = mean(squared_displacement), n = .N), by = .(cell_line, field, lag)]
msd_time <- msd[lag %in% c(1, 5), .(msd = mean(squared_displacement), n = .N), by = .(cell_line, field, frame, lag)]
msd_time <- dcast(msd_time, cell_line + field + frame ~ lag, value.var = "msd")
setnames(msd_time, c("1", "5"), c("msd1", "msd5"))
setorder(msd_time, cell_line, field, frame)
msd_time[, `:=`(msd1_roll = rolling_mean(msd1), msd5_roll = rolling_mean(msd5)), by = .(cell_line, field)]
msd_time[, alpha_1_5 := log(msd5_roll / msd1_roll) / log(5)]

fwrite(speed_time, file.path(table_dir, "speed_over_time.csv"))
fwrite(neighbor_time, file.path(table_dir, "neighbor_over_time.csv"))
fwrite(corr_time, file.path(table_dir, "directional_autocorrelation_over_time.csv"))
fwrite(msd_curve, file.path(table_dir, "msd_by_lag.csv"))
fwrite(msd_time, file.path(table_dir, "msd_persistence_over_time.csv"))

theme_report <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", strip.text = element_text(face = "bold"))
field_palette <- hue_pal()(max(fields[, .N, by = cell_line]$N))

save_plot <- function(name, plot, width = 12, height = 7) {
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = 170, bg = "white")
}

p_speed <- ggplot(speed_time, aes(frame, median_roll, group = field, color = field)) +
  geom_ribbon(data = speed_cell_time, aes(x = frame, ymin = q10_roll, ymax = q90_roll), inherit.aes = FALSE, fill = "grey70", alpha = 0.35) +
  geom_line(alpha = 0.65, linewidth = 0.55) +
  geom_line(data = speed_cell_time, aes(x = frame, y = median_roll), inherit.aes = FALSE, color = "#202020", linewidth = 0.8) +
  facet_wrap(~cell_line, scales = "free_x", ncol = 2) +
  labs(x = "Frame", y = "Speed (pixels/frame)", color = "Field", title = "Cell speed over time", subtitle = "Field medians, pooled median (black), and pooled 10th–90th percentile band; 21-frame rolling means") +
  theme_report
save_plot("01_speed_over_time", p_speed)

p_speed_iso <- ggplot(speed_time, aes(frame, group = field, color = field)) +
  geom_line(aes(y = median_roll), alpha = 0.25, linewidth = 0.5) +
  geom_line(aes(y = isolated_roll), alpha = 0.85, linewidth = 0.7) +
  facet_wrap(~cell_line, scales = "free", ncol = 2) +
  labs(x = "Frame", y = "Median speed (pixels/frame)", color = "Field", title = "Speed among conservatively isolated cells", subtitle = "Dark lines: both step endpoints have >5-pixel ROI-box separation; pale lines: all cells") +
  theme_report
save_plot("02_isolated_speed_over_time", p_speed_iso)

p_nn <- ggplot(neighbor_time, aes(frame, center_roll, group = field, color = field)) +
  geom_line(alpha = 0.7, linewidth = 0.6) +
  facet_wrap(~cell_line, scales = "free", ncol = 2) +
  labs(x = "Frame", y = "Median nearest-centroid distance (pixels)", color = "Field", title = "Nearest-neighbor distance over time", subtitle = "21-frame centered rolling means") +
  theme_report
save_plot("03_nearest_neighbor_over_time", p_nn)

iso_long <- melt(
  neighbor_time,
  id.vars = c("cell_line", "field", "frame"),
  measure.vars = c("isolated_0_roll", "isolated_5_roll", "isolated_10_roll"),
  variable.name = "threshold", value.name = "fraction"
)
iso_long[, threshold := factor(threshold, levels = c("isolated_0_roll", "isolated_5_roll", "isolated_10_roll"), labels = c(">0 px", ">5 px", ">10 px"))]
p_iso <- ggplot(iso_long, aes(frame, fraction, group = interaction(field, threshold), color = threshold)) +
  geom_line(alpha = 0.55, linewidth = 0.55) +
  facet_wrap(~cell_line, scales = "free_x", ncol = 2) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(x = "Frame", y = "Fraction of detections", color = "Minimum separation", title = "Conservatively isolated cells over time", subtitle = "Separation is the lower bound from exact ROI-derived bounding boxes; one line per field and threshold") +
  theme_report
save_plot("04_isolation_fraction_over_time", p_iso)

p_duration <- ggplot(tracks[is.finite(duration)], aes(cell_line, duration, fill = cell_line)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.65) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white") +
  geom_point(data = field_summary, aes(cell_line, duration_mean), inherit.aes = FALSE, position = position_jitter(width = 0.09), size = 2, color = "#202020") +
  coord_cartesian(ylim = quantile(tracks$duration, c(0, 0.99), na.rm = TRUE)) +
  labs(x = NULL, y = "Track duration (frames)", title = "Track-duration distributions", subtitle = "Black points are field means; violin and box plots show individual tracks") +
  guides(fill = "none") + theme_report
save_plot("05_track_duration", p_duration, height = 6)

p_corr <- ggplot(corr_time, aes(frame, correlation_roll, group = interaction(field, lag), color = factor(lag))) +
  geom_hline(yintercept = 0, color = "grey65", linewidth = 0.4) +
  geom_line(alpha = 0.55, linewidth = 0.55) +
  facet_wrap(~cell_line, scales = "free_x", ncol = 2) +
  coord_cartesian(ylim = c(-0.3, 1)) +
  labs(x = "Starting frame", y = "Velocity-direction autocorrelation", color = "Lag (frames)", title = "Directional persistence over time", subtitle = "Cosine similarity of motion vectors at lags 1, 3, and 5; 21-frame rolling means; one line per field") +
  theme_report
save_plot("06_directional_persistence_over_time", p_corr)

p_msd <- ggplot(msd_curve, aes(lag, msd, group = field, color = field)) +
  geom_line(alpha = 0.7, linewidth = 0.65) + geom_point(size = 1.5) +
  facet_wrap(~cell_line, scales = "free_y", ncol = 2) +
  scale_x_log10(breaks = c(1, 3, 5, 10, 20)) + scale_y_log10() +
  labs(x = "Lag (frames, log scale)", y = expression(MSD~(pixels^2)), color = "Field", title = "Mean squared displacement", subtitle = "Linear trajectory segments only; divisions and merges are not crossed") +
  theme_report
save_plot("07_msd_curves", p_msd)

p_alpha <- ggplot(msd_time[is.finite(alpha_1_5)], aes(frame, alpha_1_5, group = field, color = field)) +
  geom_hline(yintercept = 1, color = "grey60", linetype = "dashed", linewidth = 0.45) +
  geom_hline(yintercept = 2, color = "grey35", linetype = "dotted", linewidth = 0.45) +
  geom_line(alpha = 0.65, linewidth = 0.55) +
  facet_wrap(~cell_line, scales = "free_x", ncol = 2) +
  coord_cartesian(ylim = c(0, 2.5)) +
  labs(x = "Starting frame", y = expression(alpha[1:5]), color = "Field", title = "MSD-based persistence over time", subtitle = expression(alpha[1:5] == log(MSD[5]/MSD[1])/log(5)~"(1 = diffusive; 2 = ballistic)")) +
  theme_report
save_plot("08_msd_persistence_over_time", p_alpha)

message("Wrote summaries to ", output_dir)
