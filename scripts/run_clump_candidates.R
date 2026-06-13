#!/usr/bin/env Rscript

# General SAGUI clump-candidate tutorial workflow.
#
# Minimal use:
#
#   SAGUI_PSF_CUBE=/path/to/psf_matched_cube.fits \
#   SAGUI_RAW_CUBE=/path/to/original_flux_cube.fits \
#   SAGUI_OUT_DIR=/path/to/output \
#   SAGUI_CLUMP_RELAX=1 \
#   Rscript scripts/run_clump_candidates.R
#
# The segmentation is measured on SAGUI_PSF_CUBE. Flux-preserving regional
# photometry is measured on SAGUI_RAW_CUBE, falling back to SAGUI_PSF_CUBE if
# no raw cube is supplied.

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(sagui)
  }
  library(FITSio)
  library(ggplot2)
})

env_chr <- function(name, default = "") {
  value <- Sys.getenv(name, unset = default)
  if (!nzchar(value)) default else value
}

env_num <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, NA_character_)))
  if (length(value) == 1L && is.finite(value)) value else default
}

env_int <- function(name, default) {
  as.integer(round(env_num(name, default)))
}

env_int_null <- function(name) {
  value <- suppressWarnings(as.integer(round(as.numeric(Sys.getenv(name, NA_character_)))))
  if (length(value) == 1L && is.finite(value)) value else NULL
}

env_logical <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(name, NA_character_))
  if (is.na(value) || !nzchar(value)) return(default)
  value %in% c("true", "t", "1", "yes", "y")
}

parse_ints <- function(value, default = NULL) {
  if (is.null(value) || !nzchar(value)) return(default)
  out <- suppressWarnings(as.integer(strsplit(value, "[,; ]+")[[1]]))
  out[is.finite(out)]
}

safe_stem <- function(path) {
  x <- basename(path)
  sub("\\.[Ff][Ii][Tt][Ss](\\.gz)?$", "", x)
}

read_header_value <- function(header, key, default = NA_character_) {
  hit <- grep(paste0("^", key, "\\s*="), header, value = TRUE)
  if (!length(hit)) return(default)
  trimws(sub("^.*=\\s*'?([^'/]+).*$", "\\1", hit[[1]]))
}

filters_from_header <- function(header, nb) {
  filters <- vapply(seq_len(nb), function(i) {
    read_header_value(header, paste0("FILTER", i), NA_character_)
  }, character(1))
  missing <- !nzchar(filters) | is.na(filters)
  filters[missing] <- paste0("band", seq_len(nb))[missing]
  filters
}

read_seed_table <- function(path, dims) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  seeds <- utils::read.csv(path)
  if (all(c("clump_id", "row", "col") %in% names(seeds))) {
    out <- data.frame(
      clump_id = as.integer(seeds$clump_id),
      row = as.integer(seeds$row),
      col = as.integer(seeds$col)
    )
  } else if (all(c("clump_id", "x", "y") %in% names(seeds))) {
    add_one <- any(seeds$x == 0 | seeds$y == 0, na.rm = TRUE)
    out <- data.frame(
      clump_id = as.integer(seeds$clump_id) + if (add_one) 1L else 0L,
      row = as.integer(seeds$x) + if (add_one) 1L else 0L,
      col = as.integer(seeds$y) + if (add_one) 1L else 0L
    )
  } else {
    stop("Seed table must contain either clump_id,row,col or clump_id,x,y.")
  }

  ok <- out$row >= 1L & out$row <= dims[[1]] &
    out$col >= 1L & out$col <= dims[[2]] &
    is.finite(out$clump_id)
  out[ok, , drop = FALSE]
}

normalize01 <- function(x, q = c(0.005, 0.998)) {
  values <- x[is.finite(x)]
  if (!length(values)) return(x * 0)
  lim <- stats::quantile(values, q, na.rm = TRUE, names = FALSE)
  y <- pmin(lim[[2]], pmax(lim[[1]], x))
  dim(y) <- dim(x)
  y <- y - min(y[is.finite(y)], na.rm = TRUE)
  den <- max(y[is.finite(y)], na.rm = TRUE)
  if (!is.finite(den) || den <= 0) den <- 1
  y / den
}

asinh_stretch <- function(x, a = 10, q = c(0.005, 0.998)) {
  y <- normalize01(x, q = q)
  asinh(a * y) / asinh(a)
}

collapse_bands <- function(cube, bands) {
  bands <- bands[bands >= 1L & bands <= dim(cube)[3]]
  if (!length(bands)) bands <- seq_len(dim(cube)[3])
  Reduce(`+`, lapply(bands, function(j) cube[, , j])) / length(bands)
}

hex_from_rgb <- function(r, g, b) {
  grDevices::rgb(
    pmax(0, pmin(1, as.vector(r))),
    pmax(0, pmin(1, as.vector(g))),
    pmax(0, pmin(1, as.vector(b)))
  )
}

default_rgb_bands <- function(nb) {
  thirds <- cut(seq_len(nb), breaks = 3, labels = FALSE)
  list(
    blue = which(thirds == 1L),
    green = which(thirds == 2L),
    red = which(thirds == 3L)
  )
}

make_rgb_background <- function(cube, red_bands, green_bands, blue_bands) {
  red <- collapse_bands(cube, red_bands)
  green <- collapse_bands(cube, green_bands)
  blue <- collapse_bands(cube, blue_bands)

  r <- pmin(1, 1.10 * asinh_stretch(red, a = 10))
  g <- pmin(1, 1.05 * asinh_stretch(green, a = 10))
  b <- pmin(1, 1.30 * asinh_stretch(blue, a = 10))

  grid <- expand.grid(Row = seq_len(dim(cube)[1]), Col = seq_len(dim(cube)[2]))
  grid$fill_col <- hex_from_rgb(r, g, b)
  grid
}

accepted_clump_quality <- function(seg, clump_relax) {
  sagui::clump_evidence_table(
    clump_labels = seg$clump_footprints$clump_labels,
    score_info = seg$clump_score,
    clump_relax = clump_relax
  )
}

profile_labels_for_clumps <- function(seg, clump_ids) {
  labels <- seg$clump_profiles$profile_labels
  summary <- seg$clump_profiles$summary
  keep_regions <- summary$profile_region[summary$clump_id %in% clump_ids]
  labels[!(labels %in% keep_regions)] <- NA_integer_
  labels
}

support_from_cube <- function(cube, bands, mode = c("adaptive", "finite", "positive")) {
  mode <- match.arg(mode)
  if (identical(mode, "adaptive")) return(NULL)

  subcube <- cube[, , bands, drop = FALSE]
  finite_count <- apply(is.finite(subcube), c(1, 2), sum)
  min_frac <- env_num("SAGUI_MIN_FINITE_FRAC", 0.5)
  support <- finite_count >= ceiling(length(bands) * min_frac)

  if (identical(mode, "positive")) {
    image <- collapse_bands(cube, bands)
    support <- support & is.finite(image) & image > 0
  }

  support
}

label_components <- function(mask, connectivity = 8L) {
  mask <- is.finite(mask) & mask
  labels <- matrix(NA_integer_, nrow = nrow(mask), ncol = ncol(mask))
  offsets <- if (connectivity == 4L) {
    rbind(c(-1L, 0L), c(1L, 0L), c(0L, -1L), c(0L, 1L))
  } else {
    as.matrix(expand.grid(dr = -1L:1L, dc = -1L:1L))[ -5L, , drop = FALSE]
  }

  current <- 0L
  candidates <- which(mask & is.na(labels), arr.ind = TRUE)
  for (start_idx in seq_len(nrow(candidates))) {
    start <- candidates[start_idx, ]
    if (!is.na(labels[start[1], start[2]])) next
    current <- current + 1L
    queue_r <- integer(length(mask))
    queue_c <- integer(length(mask))
    head <- 1L
    tail <- 1L
    queue_r[tail] <- start[1]
    queue_c[tail] <- start[2]
    labels[start[1], start[2]] <- current

    while (head <= tail) {
      r <- queue_r[head]
      c <- queue_c[head]
      head <- head + 1L
      for (k in seq_len(nrow(offsets))) {
        rr <- r + offsets[k, 1]
        cc <- c + offsets[k, 2]
        if (rr < 1L || rr > nrow(mask) || cc < 1L || cc > ncol(mask)) next
        if (!mask[rr, cc] || !is.na(labels[rr, cc])) next
        tail <- tail + 1L
        queue_r[tail] <- rr
        queue_c[tail] <- cc
        labels[rr, cc] <- current
      }
    }
  }
  labels
}

central_body_mask <- function(cube, bands, support = NULL, k = 1.0, q = 0.90) {
  image <- collapse_bands(cube, bands)
  if (is.null(support)) support <- is.finite(image)
  support <- is.finite(support) & support & is.finite(image)
  values <- image[support]
  if (!length(values)) return(matrix(FALSE, nrow = nrow(image), ncol = ncol(image)))

  bg <- stats::median(values, na.rm = TRUE)
  sig <- stats::mad(values, center = bg, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(sig) || sig <= 0) sig <- stats::sd(values, na.rm = TRUE)
  if (!is.finite(sig) || sig <= 0) sig <- 1
  qcut <- stats::quantile(values, probs = min(max(q, 0), 1), na.rm = TRUE, names = FALSE)
  bright <- support & image >= max(bg + k * sig, qcut, na.rm = TRUE)
  if (!any(bright, na.rm = TRUE)) return(matrix(FALSE, nrow = nrow(image), ncol = ncol(image)))

  labels <- label_components(bright, connectivity = 8L)
  peak <- which(ifelse(support, image, NA_real_) == max(image[support], na.rm = TRUE), arr.ind = TRUE)
  peak <- peak[1, , drop = FALSE]
  peak_label <- labels[peak[1], peak[2]]
  if (!is.finite(peak_label) || is.na(peak_label)) {
    tab <- sort(table(labels[is.finite(labels)]), decreasing = TRUE)
    if (!length(tab)) return(matrix(FALSE, nrow = nrow(image), ncol = ncol(image)))
    peak_label <- as.integer(names(tab)[1])
  }
  labels == peak_label
}

plot_rgb_clump_candidates <- function(seg,
                                      cube,
                                      quality,
  clump_relax,
  red_bands,
  green_bands,
  blue_bands,
  show_labels = TRUE,
  label_size = 3.4) {
  bg <- make_rgb_background(cube, red_bands, green_bands, blue_bands)
  accepted_ids <- quality$clump_id[quality$accepted %in% TRUE]

  labels <- seg$clump_footprints$clump_labels
  labels[!(labels %in% accepted_ids)] <- NA_integer_

  fg <- expand.grid(Row = seq_len(nrow(labels)), Col = seq_len(ncol(labels)))
  fg$clump_id <- as.vector(labels)
  fg <- fg[is.finite(fg$clump_id) & fg$clump_id > 0, , drop = FALSE]
  fg <- merge(
    fg,
    quality[, c("clump_id", "quality")],
    by = "clump_id",
    all.x = TRUE,
    sort = FALSE
  )
  fg$fill_col <- ifelse(fg$quality == "strong", "#E68C3A", "#66A7B0")

  centroids <- data.frame()
  if (nrow(fg)) {
    centroids <- do.call(rbind, lapply(sort(unique(fg$clump_id)), function(id) {
      idx <- fg$clump_id == id
      data.frame(
        clump_id = id,
        Row = mean(fg$Row[idx]),
        Col = mean(fg$Col[idx]),
        stringsAsFactors = FALSE
      )
    }))
  }

  p <- ggplot() +
    geom_raster(data = bg, aes(Row, Col, fill = fill_col), alpha = 1) +
    geom_tile(data = fg, aes(Row, Col, fill = fill_col), alpha = 0.72, width = 1, height = 1) +
    scale_fill_identity(guide = "none") +
    coord_fixed(expand = FALSE) +
    theme_void(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "black", colour = NA),
      panel.background = element_rect(fill = "black", colour = NA)
    )

  if (isTRUE(show_labels) && nrow(centroids)) {
    p <- p +
      geom_text(
        data = centroids,
        aes(Row, Col, label = clump_id),
        colour = "black",
        size = label_size + 0.5,
        fontface = "bold"
      ) +
      geom_text(
        data = centroids,
        aes(Row, Col, label = clump_id),
        colour = "white",
        size = label_size,
        fontface = "bold"
      )
  }

  attr(p, "accepted_ids") <- accepted_ids
  attr(p, "clump_relax") <- clump_relax
  p
}

write_fits_label <- function(labels, path) {
  if (!requireNamespace("FITSio", quietly = TRUE)) return(invisible(FALSE))
  FITSio::writeFITSim(ifelse(is.finite(labels), labels, 0L), path, type = "single")
  invisible(TRUE)
}

run_optional_bagpipes <- function(photometry,
                                  filters,
                                  profile_summary,
                                  accepted_ids,
                                  out_dir,
                                  redshift,
                                  unit) {
  if (!requireNamespace("saguiSED", quietly = TRUE)) {
    stop("SAGUI_RUN_BAGPIPES=true requires the saguiSED package.")
  }
  if (!is.finite(redshift)) {
    stop("SAGUI_RUN_BAGPIPES=true requires SAGUI_REDSHIFT.")
  }

  filter_set_kind <- env_chr("SAGUI_FILTER_SET", "jwst_nircam")
  if (filter_set_kind == "jwst_nircam") {
    filter_set <- saguiSED::jwst_nircam_filter_set(filters)
  } else {
    stop(
      "Only SAGUI_FILTER_SET=jwst_nircam is automatic here. ",
      "For other surveys, create a filter_set and run saguiSED on the saved flux table."
    )
  }

  sed_table <- saguiSED::as_sagui_sed_table(
    photometry,
    filter_set = filter_set,
    unit = unit,
    redshift = redshift
  )

  regions <- sort(unique(profile_summary$profile_region[profile_summary$clump_id %in% accepted_ids]))
  if (!length(regions)) {
    stop("No accepted clump-profile regions available for Bagpipes.")
  }

  fit <- saguiSED::fit_region_seds(
    sed_table,
    backend = "bagpipes",
    model = saguiSED::bagpipes_model(systematic_frac = env_num("SAGUI_BAGPIPES_SYSFRAC", 0.06), metallicity = "free"),
    redshift = redshift,
    regions = regions,
    out_dir = file.path(out_dir, "bagpipes_clump_profiles"),
    python = Sys.getenv("SAGUI_PYTHON", unset = Sys.which("python3")),
    overwrite = env_logical("SAGUI_OVERWRITE_BAGPIPES", TRUE)
  )

  p <- saguiSED::plot_sed_fit_mosaic(
    fit,
    regions = regions[seq_len(min(length(regions), env_int("SAGUI_BAGPIPES_PREVIEW_N", 12L)))],
    ncol = env_int("SAGUI_BAGPIPES_NCOL", 3L),
    normalize = "none",
    transmission_height = env_num("SAGUI_BAGPIPES_TRANS_HEIGHT", 0.20),
    point_size = env_num("SAGUI_BAGPIPES_POINT_SIZE", 2.9),
    base_size = env_num("SAGUI_BAGPIPES_BASE_SIZE", 12)
  )

  saveRDS(fit, file.path(out_dir, "bagpipes_clump_profiles_fit.rds"))
  ggsave(
    file.path(out_dir, "bagpipes_clump_profiles_preview.png"),
    p,
    width = 11.0,
    height = 12.0,
    dpi = 220,
    bg = "white"
  )
  invisible(fit)
}

psf_cube_path <- env_chr("SAGUI_PSF_CUBE")
if (!nzchar(psf_cube_path) || !file.exists(psf_cube_path)) {
  stop("Set SAGUI_PSF_CUBE to a readable PSF-matched FITS cube.")
}

raw_cube_path <- env_chr("SAGUI_RAW_CUBE", psf_cube_path)
if (!file.exists(raw_cube_path)) {
  stop("SAGUI_RAW_CUBE does not exist: ", raw_cube_path)
}

stem <- env_chr("SAGUI_RUN_NAME", safe_stem(psf_cube_path))
out_dir <- env_chr("SAGUI_OUT_DIR", file.path(getwd(), paste0(stem, "_clump_candidates")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cube_psf <- FITSio::readFITS(psf_cube_path)
cube_raw <- FITSio::readFITS(raw_cube_path)
if (!identical(dim(cube_psf$imDat)[1:2], dim(cube_raw$imDat)[1:2])) {
  stop("PSF and raw cubes must have the same spatial dimensions for photometry.")
}

nb <- dim(cube_psf$imDat)[3]
filters <- filters_from_header(cube_raw$header, dim(cube_raw$imDat)[3])
bands <- parse_ints(env_chr("SAGUI_BANDS"), seq_len(nb))
support_mode <- match.arg(env_chr("SAGUI_SUPPORT_MODE", "adaptive"), c("adaptive", "finite", "positive"))
support <- support_from_cube(cube_psf$imDat, bands, mode = support_mode)
exclude_central <- env_logical("SAGUI_EXCLUDE_CENTRAL_BODY", TRUE)
central_reject <- if (exclude_central) {
  central_body_mask(
    cube_psf$imDat,
    bands = bands,
    support = support,
    k = env_num("SAGUI_CENTRAL_BODY_K", 1.0),
    q = env_num("SAGUI_CENTRAL_BODY_Q", 0.90)
  )
} else {
  NULL
}
clump_search <- if (!is.null(central_reject)) !central_reject else NULL
rgb_defaults <- default_rgb_bands(nb)
blue_bands <- parse_ints(env_chr("SAGUI_RGB_BLUE"), rgb_defaults$blue)
green_bands <- parse_ints(env_chr("SAGUI_RGB_GREEN"), rgb_defaults$green)
red_bands <- parse_ints(env_chr("SAGUI_RGB_RED"), rgb_defaults$red)

clump_relax <- env_num("SAGUI_CLUMP_RELAX", 0.5)
seed_table <- read_seed_table(env_chr("SAGUI_CLUMP_SEEDS"), dim(cube_psf$imDat)[1:2])
n_clumps <- env_int_null("SAGUI_N_CLUMPS")
max_k <- env_int_null("SAGUI_MAX_K")

message("Running Sagui clump candidates.")
message("  PSF cube: ", psf_cube_path)
message("  Raw cube: ", raw_cube_path)
message("  Output: ", out_dir)
message("  clump_relax = ", clump_relax)
message("  support_mode = ", support_mode)
message("  exclude_central_body = ", exclude_central)

seg <- sagui::segment_clumps(
  input = cube_psf,
  bands = bands,
  support = support,
  clump_search_mask = clump_search,
  clump_reject_mask = central_reject,
  clump_reject_fraction = env_num("SAGUI_CENTRAL_REJECT_FRACTION", 0.30),
  clump_seeds = seed_table,
  n_clumps = if (is.null(seed_table)) n_clumps else NULL,
  max_clumps = env_int("SAGUI_MAX_CLUMPS", 50L),
  score_quantile = env_num("SAGUI_SCORE_QUANTILE", 0.995),
  min_seed_distance = env_num("SAGUI_MIN_SEED_DISTANCE", 3),
  Ncomp_body = env_int("SAGUI_NCOMP_BODY", 80L),
  clump_levels = env_int("SAGUI_CLUMP_LEVELS", 3L),
  grow_radius = env_int("SAGUI_GROW_RADIUS", 3L),
  max_radius = env_int("SAGUI_MAX_RADIUS", 5L),
  min_pixels = env_int("SAGUI_MIN_PIXELS", 9L),
  footprint_quantile = env_num("SAGUI_FOOTPRINT_QUANTILE", 0.50),
  core_drop_frac = env_num("SAGUI_CORE_DROP_FRAC", 0.35),
  sigma_threshold = env_num("SAGUI_SIGMA_THRESHOLD", 3),
  footprint_mode = env_chr("SAGUI_FOOTPRINT_MODE", "connected"),
  contrast_weight = env_num("SAGUI_CONTRAST_WEIGHT", 0.5),
  sed_weight = env_num("SAGUI_SED_WEIGHT", 0.5),
  small_sigma = env_num("SAGUI_SMALL_SIGMA", 0.7),
  large_sigma = env_num("SAGUI_LARGE_SIGMA", 4),
  apply_quality_filter = FALSE,
  clump_relax = clump_relax,
  body_segment = env_logical("SAGUI_BODY_SEGMENT", TRUE),
  knn_k = env_int("SAGUI_KNN_K", 50L),
  auto_k = env_logical("SAGUI_AUTO_K", FALSE),
  max_k = max_k,
  feature_scale = env_chr("SAGUI_FEATURE_SCALE", "robust_col"),
  spatial_weight = env_num("SAGUI_SPATIAL_WEIGHT", 0.15),
  cluster_pretransform = env_chr("SAGUI_CLUSTER_PRETRANSFORM", "none"),
  verbose = env_logical("SAGUI_VERBOSE", TRUE)
)

quality <- accepted_clump_quality(seg, clump_relax)
accepted_ids <- quality$clump_id[quality$accepted %in% TRUE]
profile_labels <- profile_labels_for_clumps(seg, accepted_ids)
final_labels <- seg$cluster_map

photometry <- sagui::RegionPhotometry(
  cube = cube_raw,
  labels = final_labels,
  band_values = filters,
  error_fallback = env_chr("SAGUI_ERROR_FALLBACK", "mad_sky")
)
profile_photometry <- sagui::RegionPhotometry(
  cube = cube_raw,
  labels = profile_labels,
  band_values = filters,
  error_fallback = env_chr("SAGUI_ERROR_FALLBACK", "mad_sky")
)

p_rgb <- plot_rgb_clump_candidates(
  seg = seg,
  cube = cube_psf$imDat,
  quality = quality,
  clump_relax = clump_relax,
  red_bands = red_bands,
  green_bands = green_bands,
  blue_bands = blue_bands,
  show_labels = env_logical("SAGUI_SHOW_LABELS", FALSE),
  label_size = env_num("SAGUI_LABEL_SIZE", 3.4)
)
p_final <- sagui::plot_clump_diagnostics(seg, mode = "merged")

saveRDS(seg, file.path(out_dir, "clump_segmentation_result.rds"))
utils::write.csv(quality, file.path(out_dir, "clump_quality.csv"), row.names = FALSE)
utils::write.csv(quality[quality$accepted %in% TRUE, ], file.path(out_dir, "clump_candidates_accepted.csv"), row.names = FALSE)
utils::write.csv(as.data.frame(photometry$flux_long), file.path(out_dir, "final_region_fluxes_long.csv"), row.names = FALSE)
utils::write.csv(as.data.frame(profile_photometry$flux_long), file.path(out_dir, "accepted_clump_profile_fluxes_long.csv"), row.names = FALSE)
write_fits_label(seg$clump_footprints$clump_labels, file.path(out_dir, "clump_labels.fits"))
write_fits_label(profile_labels, file.path(out_dir, "accepted_clump_profile_labels.fits"))
write_fits_label(final_labels, file.path(out_dir, "final_labels.fits"))
ggsave(
  file.path(out_dir, "clump_candidates_rgb.png"),
  p_rgb,
  width = env_num("SAGUI_RGB_WIDTH", 10.5),
  height = env_num("SAGUI_RGB_HEIGHT", 5.2),
  dpi = env_int("SAGUI_RGB_DPI", 240L),
  bg = "black"
)
ggsave(
  file.path(out_dir, "final_segmentation.png"),
  p_final,
  width = env_num("SAGUI_SEG_WIDTH", 7.2),
  height = env_num("SAGUI_SEG_HEIGHT", 5.2),
  dpi = env_int("SAGUI_SEG_DPI", 240L),
  bg = "white"
)

run_bagpipes <- env_logical("SAGUI_RUN_BAGPIPES", FALSE)
if (run_bagpipes) {
  run_optional_bagpipes(
    photometry = profile_photometry,
    filters = filters,
    profile_summary = seg$clump_profiles$summary,
    accepted_ids = accepted_ids,
    out_dir = out_dir,
    redshift = env_num("SAGUI_REDSHIFT", NA_real_),
    unit = env_chr("SAGUI_FLUX_UNIT", "Jy")
  )
}

message("Accepted clumps: ", if (length(accepted_ids)) paste(accepted_ids, collapse = ", ") else "none")
message("Wrote RGB candidate map: ", file.path(out_dir, "clump_candidates_rgb.png"))
message("Wrote flux table: ", file.path(out_dir, "accepted_clump_profile_fluxes_long.csv"))
