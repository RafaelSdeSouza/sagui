#!/usr/bin/env Rscript

# Minimal clump-aware SAGUI demo using the JellyFish J11 cube.
#
# The package function is generic; this script only supplies one local example.
# Segmentation is run on the PSF-matched cube. If an original cube is supplied,
# flux-preserving regional photometry is extracted from that original cube.

suppressPackageStartupMessages({
  if (identical(tolower(Sys.getenv("SAGUI_USE_LOCAL", "true")), "true") &&
      file.exists("DESCRIPTION") &&
      requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(sagui)
  }
  library(FITSio)
  library(ggplot2)
})

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = default)
  if (!nzchar(value)) default else value
}

env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(value)) default else value
}

env_num <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(value)) default else value
}

env_bool <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  value %in% c("1", "true", "t", "yes", "y")
}

read_optional_clump_catalog <- function(path, dims) {
  if (!file.exists(path)) return(NULL)
  x <- utils::read.csv(path)
  if (!all(c("clump_id", "x", "y") %in% names(x))) return(NULL)
  add_one <- any(x$x == 0 | x$y == 0, na.rm = TRUE)
  out <- data.frame(
    clump_id = x$clump_id,
    row = as.integer(x$x) + if (add_one) 1L else 0L,
    col = as.integer(x$y) + if (add_one) 1L else 0L
  )
  ok <- out$row >= 1L & out$row <= dims[1] & out$col >= 1L & out$col <= dims[2]
  out[ok, , drop = FALSE]
}

header_value <- function(header, key, default = NA_character_) {
  hit <- grep(paste0("^", key, "\\s*="), header, value = TRUE)
  if (!length(hit)) return(default)
  sub("^.*=\\s*'?([^'/]+).*$", "\\1", hit[[1]])
}

filters_from_header <- function(header, nb) {
  values <- vapply(seq_len(nb), function(i) {
    trimws(header_value(header, paste0("FILTER", i), NA_character_))
  }, character(1))
  missing <- !nzchar(values) | is.na(values)
  values[missing] <- paste0("band", seq_len(nb))[missing]
  values
}

project_dir <- env_chr("JELLYFISH_PROJECT_DIR", "/Users/rd23aag/Documents/GitHub/JellyFish_Andressa")
psf_cube_path <- env_chr("SAGUI_CLUMP_PSF_CUBE", file.path(project_dir, "J11", "datacube_filtered_psfmatched.fits"))
original_cube_path <- env_chr("SAGUI_CLUMP_ORIGINAL_CUBE", file.path(project_dir, "J11", "datacube_filtered.fits"))
clump_catalog_path <- env_chr("SAGUI_CLUMP_CATALOG", file.path(project_dir, "J11", "a2744_clumps_pixels_J11.csv"))
out_dir <- env_chr("SAGUI_CLUMP_DEMO_OUT", file.path(project_dir, "J11", "sagui_clump_mode_demo"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading PSF-matched cube: ", psf_cube_path)
cube <- FITSio::readFITS(psf_cube_path)
original_cube <- if (file.exists(original_cube_path)) FITSio::readFITS(original_cube_path) else cube

clump_seeds <- NULL
if (env_bool("SAGUI_CLUMP_USE_CATALOG", FALSE)) {
  clump_seeds <- read_optional_clump_catalog(clump_catalog_path, dim(cube$imDat)[1:2])
  message("Using catalogue clump seeds: ", if (is.null(clump_seeds)) 0L else length(unique(clump_seeds$clump_id)))
}

seg <- sagui::segment_clumps(
  input = cube,
  Ncomp_body = env_int("SAGUI_CLUMP_NCOMP_BODY", 80L),
  clump_seeds = clump_seeds,
  n_clumps = if (is.null(clump_seeds)) env_int("SAGUI_CLUMP_N_CLUMPS", 32L) else NULL,
  max_clumps = env_int("SAGUI_CLUMP_MAX_CLUMPS", 50L),
  score_quantile = env_num("SAGUI_CLUMP_SCORE_Q", 0.995),
  min_seed_distance = env_num("SAGUI_CLUMP_MIN_SEED_DISTANCE", 3),
  clump_levels = env_int("SAGUI_CLUMP_LEVELS", 3L),
  grow_radius = env_int("SAGUI_CLUMP_GROW_RADIUS", 3L),
  max_radius = env_int("SAGUI_CLUMP_MAX_RADIUS", 7L),
  min_pixels = env_int("SAGUI_CLUMP_MIN_PIXELS", 9L),
  footprint_quantile = env_num("SAGUI_CLUMP_FOOTPRINT_Q", 0.10),
  contrast_weight = env_num("SAGUI_CLUMP_CONTRAST_WEIGHT", 0.5),
  sed_weight = env_num("SAGUI_CLUMP_SED_WEIGHT", 0.5),
  small_sigma = env_num("SAGUI_CLUMP_SMALL_SIGMA", 0.7),
  large_sigma = env_num("SAGUI_CLUMP_LARGE_SIGMA", 4),
  feature_scale = env_chr("SAGUI_CLUMP_FEATURE_SCALE", "robust_col"),
  spatial_weight = env_num("SAGUI_CLUMP_SPATIAL_WEIGHT", 0.15),
  knn_k = env_int("SAGUI_CLUMP_KNN", 50L),
  verbose = TRUE
)

saveRDS(seg, file.path(out_dir, "j11_segment_clumps_result.rds"))
utils::write.csv(seg$clump_seeds, file.path(out_dir, "j11_detected_clump_seeds.csv"), row.names = FALSE)
utils::write.csv(seg$clump_footprints$footprint_summary, file.path(out_dir, "j11_clump_footprints.csv"), row.names = FALSE)
utils::write.csv(seg$clump_profiles$summary, file.path(out_dir, "j11_clump_profile_regions.csv"), row.names = FALSE)

FITSio::writeFITSim(ifelse(is.finite(seg$cluster_map), seg$cluster_map, 0L), file.path(out_dir, "j11_segment_clumps_map.fits"), type = "single")
FITSio::writeFITSim(ifelse(is.finite(seg$clump_profiles$profile_labels), seg$clump_profiles$profile_labels, 0L), file.path(out_dir, "j11_clump_profiles_map.fits"), type = "single")
FITSio::writeFITSim(ifelse(is.finite(seg$clump_score$score), seg$clump_score$score, 0), file.path(out_dir, "j11_clump_score.fits"), type = "single")

ggplot2::ggsave(file.path(out_dir, "j11_clump_score.png"), sagui::plot_clump_diagnostics(seg, mode = "score"), width = 7.2, height = 4.2, dpi = 260, bg = "black")
ggplot2::ggsave(file.path(out_dir, "j11_clump_overlay.png"), sagui::plot_clump_diagnostics(seg, mode = "overlay"), width = 7.2, height = 4.2, dpi = 260, bg = "black")
ggplot2::ggsave(file.path(out_dir, "j11_clump_profiles.png"), sagui::plot_clump_diagnostics(seg, mode = "profiles"), width = 7.2, height = 4.2, dpi = 260, bg = "white")
ggplot2::ggsave(file.path(out_dir, "j11_segment_clumps_map.png"), sagui::plot_clump_diagnostics(seg, mode = "merged"), width = 7.2, height = 4.2, dpi = 260, bg = "white")

filters <- filters_from_header(original_cube$header, dim(original_cube$imDat)[3])
phot <- sagui::RegionPhotometry(
  cube = original_cube,
  labels = seg$cluster_map,
  band_values = filters,
  error_fallback = "mad_sky"
)
utils::write.csv(as.data.frame(phot$flux_long), file.path(out_dir, "j11_segment_clumps_region_seds_long.csv"), row.names = FALSE)

message("Wrote demo products to: ", out_dir)
