#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(guara)
  library(ragg)
})

repo_dir <- normalizePath(file.path(getwd()), mustWork = TRUE)
crp_dir <- Sys.getenv(
  "CRP8_SEGMENTATION_DIR",
  unset = "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
)

web_out_dir <- file.path(repo_dir, "images", "examples", "target_snr", "sagui10")
paper_out_dir <- Sys.getenv(
  "SAGUI10_TARGET_SNR_PAPER_DIR",
  unset = file.path(crp_dir, "results/figures/paper_repro/target_snr_selection/sagui10")
)
dir.create(web_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_out_dir, recursive = TRUE, showWarnings = FALSE)

cube_path <- Sys.getenv(
  "SAGUI10_PSF_CUBE",
  unset = file.path(crp_dir, "data/PSF_matched/datacube_sagui-10_psfmatched.fits")
)
if (!file.exists(cube_path)) {
  stop("Sagui-10 PSF-matched cube not found: ", cube_path)
}

source_sagui <- function() {
  files <- list.files(file.path(repo_dir, "R"), full.names = TRUE)
  files <- files[endsWith(files, ".R")]
  files <- setdiff(files, file.path(repo_dir, c("R/RcppExports.R", "R/zzz.R")))
  for (f in files) {
    source(f)
  }
}

palette_van_gogh_div <- function(n = 256) {
  stops <- c(
    "#0A1026",
    "#182A6B",
    "#2F5F9A",
    "#4FA7A6",
    "#D8E1D3",
    "#F1C76A",
    "#A96A2A"
  )
  grDevices::colorRampPalette(stops, space = "Lab")(n)
}

make_clean_map_theme <- function(base_size = 13) {
  theme_void(base_size = base_size) +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

compute_region_snr <- function(cube_fits, labels, bands) {
  sed <- extract_region_sed(
    cube = cube_fits,
    labels = labels,
    band_values = bands,
    error_fallback = "flux_over_sqrt_n"
  )
  .compute_region_snr_from_sed(
    sed$flux_long,
    snr_stat = "integrated"
  )
}

source_sagui()

set.seed(42)
cube_fits <- FITSio::readFITS(cube_path)
cube <- cube_fits$imDat
bands <- c(
  "F090W", "F115W", "F150W", "F182M", "F200W", "F210M",
  "F277W", "F335M", "F356W", "F410M", "F444W"
)

collapsed <- guara::collapse_white_light(cube, kclip = 2)
starlet_dec <- guara::starlet_mask(collapsed, J = 5)
starlet_rec <- guara::starlet_reconstruct(
  starlet_dec,
  keep_scales = 2:5,
  include_coarse = FALSE,
  denoise_k = 2,
  mode = "soft"
)
support_mask <- is.finite(starlet_rec) & starlet_rec > 0

masked_fits <- cube_fits
masked_fits$imDat <- guara::mask_cube(cube, support_mask, mode = "na")

candidate_n <- c(40L, 20L, 15L, 10L, 8L, 5L, 4L, 3L, 2L)
candidate_results <- vector("list", length(candidate_n))
names(candidate_results) <- paste0("n", candidate_n)

for (i in seq_along(candidate_n)) {
  ncomp <- candidate_n[i]
  message("Running Sagui-10 candidate Ncomp=", ncomp)
  seg <- segment_regions(
    input = masked_fits,
    Ncomp = ncomp,
    cluster_pretransform = "none",
    use_starlet_mask = FALSE,
    collapse_fn = function(x) guara::collapse_white_light(x, kclip = 2),
    mask_mode = "na"
  )

  snr <- compute_region_snr(cube_fits, seg$cluster_map, bands)
  candidate_results[[i]] <- list(
    Ncomp = ncomp,
    segmentation = seg,
    region_snr = snr,
    snr_min = min(snr, na.rm = TRUE),
    snr_median = stats::median(snr, na.rm = TRUE),
    snr_q16 = unname(stats::quantile(snr, 0.16, na.rm = TRUE)),
    snr_q84 = unname(stats::quantile(snr, 0.84, na.rm = TRUE))
  )
}

candidate_table <- do.call(
  rbind,
  lapply(candidate_results, function(x) {
    data.frame(
      Ncomp = x$Ncomp,
      min_snr = x$snr_min,
      median_snr = x$snr_median,
      q16_snr = x$snr_q16,
      q84_snr = x$snr_q84
    )
  })
)
candidate_table <- candidate_table[order(candidate_table$Ncomp, decreasing = TRUE), ]

targets <- c(5, 10, 15, 20)
selection_table <- do.call(
  rbind,
  lapply(targets, function(target) {
    eligible <- candidate_table[candidate_table$min_snr >= target, , drop = FALSE]
    if (!nrow(eligible)) {
      selected <- candidate_table[which.min(candidate_table$Ncomp), , drop = FALSE]
      status <- "not_reached"
    } else {
      selected <- eligible[which.max(eligible$Ncomp), , drop = FALSE]
      status <- "reached"
    }

    data.frame(
      target_min_snr = target,
      selected_Ncomp = selected$Ncomp,
      achieved_min_snr = selected$min_snr,
      achieved_median_snr = selected$median_snr,
      status = status
    )
  })
)

write.csv(
  candidate_table,
  file.path(web_out_dir, "sagui10_target_snr_candidates.csv"),
  row.names = FALSE
)
write.csv(
  selection_table,
  file.path(web_out_dir, "sagui10_target_snr_selection.csv"),
  row.names = FALSE
)

panel_plot <- function(target) {
  selected_n <- selection_table$selected_Ncomp[selection_table$target_min_snr == target]
  item <- candidate_results[[paste0("n", selected_n)]]
  map <- item$segmentation$cluster_map
  n_regions <- length(unique(as.integer(map[is.finite(map) & map > 0])))

  plot_region_map(
    list(cluster_map = map),
    palette = palette_van_gogh_div(max(3L, n_regions)),
    border_color = "black",
    border_linewidth = 0.85,
    background_color = "white"
  ) +
    make_clean_map_theme(base_size = 12)
}

panels <- lapply(targets, panel_plot)

for (i in seq_along(targets)) {
  file_name <- sprintf("sagui10_target_snr_%02d.png", targets[i])
  for (target_path in c(file.path(web_out_dir, file_name), file.path(paper_out_dir, file_name))) {
    ragg::agg_png(target_path, width = 3.0, height = 3.0, units = "in", res = 400, background = "white")
    print(panels[[i]])
    grDevices::dev.off()
  }
}

write.csv(
  candidate_table,
  file.path(paper_out_dir, "sagui10_target_snr_candidates.csv"),
  row.names = FALSE
)
write.csv(
  selection_table,
  file.path(paper_out_dir, "sagui10_target_snr_selection.csv"),
  row.names = FALSE
)

message("Wrote clean website images to: ", web_out_dir)
message("Wrote clean paper images to: ", paper_out_dir)
message("Selection table: ", file.path(paper_out_dir, "sagui10_target_snr_selection.csv"))
