#!/usr/bin/env Rscript

# Demo workflow for running sagui on a cutout cube and selecting a GC candidate.
# Edit the configuration block below to tune the segmentation and the GC heuristic.

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) &&
      dir.exists("/Users/rd23aag/Documents/GitHub/sagui")) {
    pkgload::load_all("/Users/rd23aag/Documents/GitHub/sagui", quiet = TRUE)
  } else {
    library(sagui)
  }

  library(FITSio)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
input_path <- "CO_84_131.fits"
output_dir <- "sagui_outputs"

seg_ncomp <- 20
starlet_J <- 5
starlet_scales <- 2:5

# GC candidate heuristic:
# favor bright, central, and relatively compact regions.
# No roundness prior is imposed: the final GC mask keeps the segmented morphology.
gc_search_radius_px <- 10
peak_weight <- 0.45
flux_weight <- 0.25
central_weight <- 0.20
compact_weight <- 0.10

# Optional second pass inside the selected GC mask.
# Keep this NULL unless you explicitly want to split the GC internally.
gc_resegment_ncomp <- NULL

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
scale01 <- function(x) {
  if (length(x) == 0L) {
    return(x)
  }
  rng <- range(x, na.rm = TRUE, finite = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

plot_logical_mask <- function(mask) {
  df <- data.frame(
    x = rep(seq_len(nrow(mask)), ncol(mask)),
    y = rep(seq_len(ncol(mask)), each = nrow(mask)),
    value = as.vector(mask)
  )

  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    coord_fixed() +
    scale_fill_manual(values = c("FALSE" = "black", "TRUE" = "#FFDE38")) +
    theme_void() +
    theme(legend.position = "none")
}

plot_numeric_map <- function(mat, title = NULL) {
  df <- data.frame(
    x = rep(seq_len(nrow(mat)), ncol(mat)),
    y = rep(seq_len(ncol(mat)), each = nrow(mat)),
    value = as.vector(mat)
  )

  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    coord_fixed() +
    scale_fill_viridis_c(na.value = "black") +
    theme_void() +
    labs(title = title)
}

crop_cube_to_mask <- function(cube, mask, pad = 1L) {
  idx <- which(mask, arr.ind = TRUE)
  if (!nrow(idx)) {
    stop("Mask is empty; cannot crop cube.")
  }

  x1 <- max(1L, min(idx[, 1]) - pad)
  x2 <- min(dim(cube)[1], max(idx[, 1]) + pad)
  y1 <- max(1L, min(idx[, 2]) - pad)
  y2 <- min(dim(cube)[2], max(idx[, 2]) + pad)

  cube_crop <- cube[x1:x2, y1:y2, , drop = FALSE]
  mask_crop <- mask[x1:x2, y1:y2, drop = FALSE]
  cube_crop[!array(rep(as.vector(mask_crop), times = dim(cube_crop)[3]), dim = dim(cube_crop))] <- NA_real_

  list(
    cube = cube_crop,
    mask = mask_crop,
    bbox = c(xmin = x1, xmax = x2, ymin = y1, ymax = y2)
  )
}

summarize_regions <- function(cluster_map, image) {
  labels <- sort(unique(na.omit(as.vector(cluster_map))))
  if (!length(labels)) {
    stop("No positive labels found in cluster_map.")
  }

  cx0 <- (nrow(cluster_map) + 1) / 2
  cy0 <- (ncol(cluster_map) + 1) / 2

  stats <- lapply(labels, function(label_id) {
    idx <- which(cluster_map == label_id, arr.ind = TRUE)
    vals <- image[cluster_map == label_id]

    x_center <- mean(idx[, 1])
    y_center <- mean(idx[, 2])
    distance_to_center <- sqrt((x_center - cx0)^2 + (y_center - cy0)^2)

    data.frame(
      region = label_id,
      n_pix = nrow(idx),
      x_center = x_center,
      y_center = y_center,
      distance_to_center = distance_to_center,
      peak_collapsed = max(vals, na.rm = TRUE),
      sum_collapsed = sum(vals, na.rm = TRUE),
      mean_collapsed = mean(vals, na.rm = TRUE)
    )
  })

  do.call(rbind, stats)
}

choose_gc_candidate <- function(region_stats, search_radius_px) {
  region_stats$peak_score <- scale01(region_stats$peak_collapsed)
  region_stats$flux_score <- scale01(region_stats$sum_collapsed)
  region_stats$central_score <- 1 - scale01(region_stats$distance_to_center)
  region_stats$compact_score <- 1 - scale01(region_stats$n_pix)

  region_stats$gc_score <- peak_weight * region_stats$peak_score +
    flux_weight * region_stats$flux_score +
    central_weight * region_stats$central_score +
    compact_weight * region_stats$compact_score

  local <- region_stats[region_stats$distance_to_center <= search_radius_px, , drop = FALSE]
  if (nrow(local) == 0L) {
    local <- region_stats
  }

  local[which.max(local$gc_score), , drop = FALSE]
}

# -----------------------------------------------------------------------------
# Run sagui
# -----------------------------------------------------------------------------
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

x <- FITSio::readFITS(input_path)

seg <- segment_regions(
  input = x,
  Ncomp = seg_ncomp,
  use_starlet_mask = TRUE,
  starlet_J = starlet_J,
  starlet_scales = starlet_scales,
  hclust_method = "ward.D2"
)

sed <- extract_region_sed(
  cube = x,
  labels = seg$cluster_map,
  band_values = if (!is.null(x$axDat)) FITSio::axVec(3, x$axDat) else NULL
)

region_stats <- summarize_regions(seg$cluster_map, seg$collapsed)
gc_candidate <- choose_gc_candidate(region_stats, gc_search_radius_px)
gc_region <- gc_candidate$region[[1]]

# This is the final point-source-compatible segmentation:
# a single segmented region selected as the best GC candidate.
gc_mask <- seg$cluster_map == gc_region
gc_final_map <- matrix(NA_integer_, nrow(seg$cluster_map), ncol(seg$cluster_map))
gc_final_map[gc_mask] <- 1L
gc_score_map <- matrix(NA_real_, nrow(seg$cluster_map), ncol(seg$cluster_map))
gc_score_map[] <- 0
gc_score_map[gc_mask] <- gc_candidate$gc_score[[1]]

gc_cube <- mask_cube(x$imDat, gc_mask, mode = "na")
gc_crop <- crop_cube_to_mask(x$imDat, gc_mask, pad = 1L)
gc_crop_map <- matrix(NA_integer_, nrow(gc_crop$mask), ncol(gc_crop$mask))
gc_crop_map[gc_crop$mask] <- 1L
gc_input <- x
gc_input$imDat <- gc_crop$cube

gc_sed <- extract_region_sed(
  cube = gc_crop$cube,
  labels = gc_crop_map,
  band_values = if (!is.null(x$axDat)) FITSio::axVec(3, x$axDat) else NULL
)

gc_seg <- NULL
if (!is.null(gc_resegment_ncomp) && sum(gc_mask, na.rm = TRUE) > gc_resegment_ncomp) {
  gc_seg <- segment_regions(
    input = gc_input,
    Ncomp = gc_resegment_ncomp,
    use_starlet_mask = FALSE,
    hclust_method = "ward.D2"
  )
}

# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------
FITSio::writeFITSim(seg$cluster_map,
  file = file.path(output_dir, sprintf("CO_84_131_cluster_map_N%d.fits", seg_ncomp))
)
FITSio::writeFITSim(seg$collapsed,
  file = file.path(output_dir, "CO_84_131_collapsed_white_light.fits")
)
FITSio::writeFITSim(seg$mask * 1,
  file = file.path(output_dir, "CO_84_131_starlet_mask.fits")
)
FITSio::writeFITSim(gc_mask * 1,
  file = file.path(output_dir, "CO_84_131_gc_candidate_mask.fits")
)
FITSio::writeFITSim(gc_final_map,
  file = file.path(output_dir, "CO_84_131_gc_final_map.fits")
)
FITSio::writeFITSim(gc_crop$cube,
  file = file.path(output_dir, "CO_84_131_gc_cropped_cube.fits")
)
FITSio::writeFITSim(gc_crop$mask * 1,
  file = file.path(output_dir, "CO_84_131_gc_cropped_mask.fits")
)
if (!is.null(gc_seg)) {
  FITSio::writeFITSim(gc_seg$cluster_map,
    file = file.path(output_dir, sprintf("CO_84_131_gc_subcluster_map_N%d.fits", gc_resegment_ncomp))
  )
}

write.csv(region_stats,
  file = file.path(output_dir, "CO_84_131_region_summary.csv"),
  row.names = FALSE
)
write.csv(gc_candidate,
  file = file.path(output_dir, "CO_84_131_gc_candidate.csv"),
  row.names = FALSE
)
write.csv(as.data.frame(as.list(gc_crop$bbox)),
  file = file.path(output_dir, "CO_84_131_gc_crop_bbox.csv"),
  row.names = FALSE
)
write.csv(gc_sed$flux_long,
  file = file.path(output_dir, "CO_84_131_gc_sed_flux_long.csv"),
  row.names = FALSE
)
write.csv(gc_sed$flux_wide,
  file = file.path(output_dir, "CO_84_131_gc_sed_flux_wide.csv"),
  row.names = FALSE
)
write.csv(sed$flux_long,
  file = file.path(output_dir, "CO_84_131_region_sed_flux_long.csv"),
  row.names = FALSE
)
write.csv(sed$flux_wide,
  file = file.path(output_dir, "CO_84_131_region_sed_flux_wide.csv"),
  row.names = FALSE
)

ggsave(
  filename = file.path(output_dir, sprintf("CO_84_131_region_map_N%d.png", seg_ncomp)),
  plot = plot_region_map(seg, palette = "magma", border_color = "black", border_linewidth = 0.7, background_color = "transparent"),
  width = 7, height = 5, dpi = 220, bg = "white"
)
ggsave(
  filename = file.path(output_dir, "CO_84_131_starlet_mask.png"),
  plot = plot_logical_mask(seg$mask),
  width = 6, height = 4, dpi = 220, bg = "white"
)
ggsave(
  filename = file.path(output_dir, "CO_84_131_gc_candidate_mask.png"),
  plot = plot_logical_mask(gc_mask),
  width = 6, height = 4, dpi = 220, bg = "white"
)
ggsave(
  filename = file.path(output_dir, "CO_84_131_gc_final_map.png"),
  plot = plot_region_map(
    list(cluster_map = gc_final_map),
    palette = c("#FFDE38", "#FFDE38"),
    border_color = "black",
    border_linewidth = 0.8,
    background_color = "transparent"
  ),
  width = 6, height = 4, dpi = 220, bg = "white"
)
ggsave(
  filename = file.path(output_dir, "CO_84_131_gc_cropped_mask.png"),
  plot = plot_logical_mask(gc_crop$mask),
  width = 4, height = 4, dpi = 220, bg = "white"
)
if (!is.null(gc_seg)) {
  ggsave(
    filename = file.path(output_dir, sprintf("CO_84_131_gc_subcluster_map_N%d.png", gc_resegment_ncomp)),
    plot = plot_region_map(gc_seg, palette = "magma", border_color = "black", border_linewidth = 0.7, background_color = "transparent"),
    width = 6, height = 4, dpi = 220, bg = "white"
  )
}
ggsave(
  filename = file.path(output_dir, "CO_84_131_collapsed_white_light.png"),
  plot = plot_numeric_map(seg$collapsed, title = "Collapsed white-light image"),
  width = 6, height = 4, dpi = 220, bg = "white"
)

cat("Finished sagui cutout demo.\n")
cat("Input:", normalizePath(input_path), "\n")
cat("Output directory:", normalizePath(output_dir), "\n")
cat("Detected regions:", length(unique(na.omit(as.vector(seg$cluster_map)))), "\n")
cat("GC candidate region:", gc_region, "\n")
cat("GC candidate score:", round(gc_candidate$gc_score[[1]], 4), "\n")
cat("GC pixels kept:", sum(gc_mask, na.rm = TRUE), "\n")
cat("GC crop bbox:", paste(names(gc_crop$bbox), gc_crop$bbox, collapse = ", "), "\n")
if (!is.null(gc_seg)) {
  cat("GC sub-segmentation regions:", length(unique(na.omit(as.vector(gc_seg$cluster_map)))), "\n")
}
