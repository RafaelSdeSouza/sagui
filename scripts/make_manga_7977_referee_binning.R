#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(capivara)
  library(ggplot2)
  library(patchwork)
  library(reshape2)
})

args <- commandArgs(trailingOnly = TRUE)

arg_or <- function(i, default) {
  if (length(args) < i) {
    return(default)
  }
  value <- args[[i]]
  if (is.null(value) || length(value) == 0 || is.na(value) || !nzchar(value)) {
    return(default)
  }
  value
}

fits_path <- arg_or(1, "/private/tmp/manga-7977-12701-LOGCUBE.fits")
if (!file.exists(fits_path)) {
  fits_path <- "/Users/rd23aag/Downloads/manga-7977-12701-LOGCUBE.fits"
}

out_dir <- arg_or(2, file.path(getwd(), "outputs", "referee_manga_7977"))
ncomp_arg <- arg_or(3, "20")
ncomp <- suppressWarnings(as.integer(ncomp_arg))
if (!is.finite(ncomp) || ncomp < 1L) {
  stop(
    "Ncomp must be a positive integer. Received: ", sQuote(ncomp_arg), "\n",
    "Usage: Rscript scripts/make_manga_7977_referee_binning.R ",
    "[fits_path] [out_dir] [Ncomp]",
    call. = FALSE
  )
}

if (!file.exists(fits_path)) {
  stop(
    "Input FITS not found. Copy manga-7977-12701-LOGCUBE.fits to /private/tmp ",
    "or pass its path as the first argument.",
    call. = FALSE
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

van_gogh <- function(n = 256) {
  grDevices::colorRampPalette(
    c("#213E60", "#345E9D", "#69A7B0", "#F4F2EF", "#F2C75C", "#E68C3A"),
    space = "Lab"
  )(n)
}

to_df <- function(mat, name) {
  df <- reshape2::melt(mat, varnames = c("Row", "Col"), value.name = "value")
  df$panel <- name
  df
}

asinh01 <- function(x) {
  y <- asinh(x)
  rng <- range(y[is.finite(y)], na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) <= 0) {
    return(rep(NA_real_, length(x)))
  }
  (y - rng[1]) / diff(rng)
}

style_void <- function() {
  theme_void() +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(hjust = 0.5, face = "bold", colour = "black", size = 13),
      plot.margin = margin(3, 3, 3, 3)
    )
}

message("Reading ", fits_path)
cube <- FITSio::readFITS(fits_path)

dim0 <- dim(cube$imDat)
if (length(dim0) != 3L) {
  stop("Expected a 3D LOGCUBE, got dimensions: ", paste(dim0, collapse = " x "), call. = FALSE)
}

# FITS readers may place the spectral axis first. Capivara expects [x, y, lambda].
wave_dim <- which.max(dim0)
if (wave_dim != 3L) {
  cube$imDat <- aperm(cube$imDat, c(setdiff(seq_len(3), wave_dim), wave_dim))
  message("Permuted cube from ", paste(dim0, collapse = " x "), " to ", paste(dim(cube$imDat), collapse = " x "))
}

message("Running Capivara sparse-Ward segmentation with Ncomp = ", ncomp)
res <- segment_large(
  cube,
  Ncomp = 40,
  use_starlet_mask = TRUE,
  starlet_J = 5,
  starlet_scales = 2:5,
  include_coarse = FALSE,
  denoise_k = 0,
  positive_only = TRUE,
  mask_mode = "na",
  knn_k = 50,
  auto_k = TRUE,
  feature_scale = "robust_col",
  spatial_weight = 0.06,
  valid_mode = "signal",
  verbose = TRUE
)

prefix <- file.path(out_dir, sprintf("manga_7977_12701_n%02d", ncomp))

bin_plot <- plot_cluster(res, palette = van_gogh(ncomp)) +
  style_void()

ggsave(
  filename = paste0(prefix, "_binning_scheme.png"),
  plot = bin_plot,
  width = 4.2,
  height = 4.2,
  dpi = 450,
  bg = "white"
)

ggsave(
  filename = paste0(prefix, "_binning_scheme.pdf"),
  plot = bin_plot,
  width = 4.2,
  height = 4.2,
  device = cairo_pdf,
  bg = "white"
)

star <- res$starlet_info
white <- star$collapsed
white[white <= 0] <- NA_real_
white_df <- to_df(white, "white")
white_df$value <- asinh01(white_df$value)

mask_df <- to_df(star$mask * 1, "mask")

white_plot <- ggplot(white_df, aes(Row, Col, fill = value)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_gradientn(colours = van_gogh(256), na.value = "white") +
  labs(title = "White-light image") +
  style_void() +
  theme(legend.position = "none")

mask_plot <- ggplot(mask_df, aes(Row, Col, fill = factor(value))) +
  geom_raster() +
  coord_fixed() +
  scale_fill_manual(values = c("0" = "white", "1" = "#213E60"), na.value = "white") +
  labs(title = "Starlet support") +
  style_void() +
  theme(legend.position = "none")

bin_plot_titled <- bin_plot + labs(title = sprintf("Capivara bins, N = %d", res$Ncomp))

panel <- white_plot + mask_plot + bin_plot_titled +
  plot_layout(nrow = 1, widths = c(1, 1, 1)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(
  filename = paste0(prefix, "_referee_panel.png"),
  plot = panel,
  width = 9.2,
  height = 3.15,
  dpi = 450,
  bg = "white"
)

ggsave(
  filename = paste0(prefix, "_referee_panel.pdf"),
  plot = panel,
  width = 9.2,
  height = 3.15,
  device = cairo_pdf,
  bg = "white"
)

saveRDS(res, paste0(prefix, "_result.rds"))

write.csv(
  data.frame(
    plateifu = "7977-12701",
    Ncomp_requested = ncomp,
    Ncomp_returned = res$Ncomp,
    valid_spaxels = res$backend_info$valid_pixels,
    mask_spaxels = sum(star$mask, na.rm = TRUE),
    backend = res$backend,
    knn_k = res$backend_info$knn_k,
    feature_scale = res$backend_info$feature_scale,
    spatial_weight = res$backend_info$spatial_weight
  ),
  paste0(prefix, "_metadata.csv"),
  row.names = FALSE
)

message("Wrote:")
message("  ", paste0(prefix, "_binning_scheme.png"))
message("  ", paste0(prefix, "_binning_scheme.pdf"))
message("  ", paste0(prefix, "_referee_panel.png"))
message("  ", paste0(prefix, "_referee_panel.pdf"))
