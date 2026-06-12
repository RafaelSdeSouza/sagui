#' Compute a compact SED-clump score
#'
#' @description
#' Builds a score map for compact structures whose photometric SEDs differ from
#' the local surrounding galaxy light. The score combines a local-contrast term
#' from a collapsed image with a local SED-anomaly term measured relative to a
#' smoothed version of the cube.
#'
#' @param input 3-D array or FITS-like list with `imDat`.
#' @param support Optional logical foreground support mask. If `NULL`, all
#'   finite pixels in the collapsed image are used.
#' @param bands Bands used to build the score. `NULL` uses all bands.
#' @param contrast_weight Weight assigned to the local-contrast term.
#' @param sed_weight Weight assigned to the local SED-anomaly term.
#' @param small_sigma Gaussian smoothing scale for the compact image term.
#' @param large_sigma Gaussian smoothing scale for the local-background term.
#'
#' @return A list with `score`, `contrast`, `sed_anomaly`, and `collapsed`.
#' @export
compute_clump_score <- function(input,
                                support = NULL,
                                bands = NULL,
                                contrast_weight = 0.5,
                                sed_weight = 0.5,
                                small_sigma = 0.7,
                                large_sigma = 4) {
  cube <- .clump_cube(input)
  bands <- .clump_resolve_bands(cube, bands)
  collapsed <- .clump_collapse_bands(cube, bands)

  if (is.null(support)) {
    support <- is.finite(collapsed)
  }
  support <- is.finite(support) & support

  image <- sign(collapsed) * log1p(abs(collapsed))
  small <- .clump_smooth_image(image, sigma = small_sigma)
  large <- .clump_smooth_image(image, sigma = large_sigma)
  contrast <- pmax(small - large, 0)
  contrast <- .clump_normalize01(contrast, mask = support)

  scaled <- .clump_robust_scale_cube(cube[, , bands, drop = FALSE], support = support)
  sed_anomaly <- matrix(0, nrow = dim(cube)[1], ncol = dim(cube)[2])
  n_used <- 0L
  for (j in seq_len(dim(scaled)[3])) {
    band <- scaled[, , j]
    local <- .clump_smooth_image(band, sigma = large_sigma)
    delta <- band - local
    delta[!is.finite(delta)] <- 0
    sed_anomaly <- sed_anomaly + delta^2
    n_used <- n_used + 1L
  }
  sed_anomaly <- sqrt(sed_anomaly / max(1L, n_used))
  sed_anomaly <- .clump_normalize01(sed_anomaly, mask = support)

  contrast_weight <- max(0, as.numeric(contrast_weight))
  sed_weight <- max(0, as.numeric(sed_weight))
  if (!is.finite(contrast_weight)) contrast_weight <- 0.5
  if (!is.finite(sed_weight)) sed_weight <- 0.5
  total_weight <- contrast_weight + sed_weight
  if (!is.finite(total_weight) || total_weight <= 0) {
    contrast_weight <- 0.5
    sed_weight <- 0.5
    total_weight <- 1
  }

  score <- (contrast_weight * contrast + sed_weight * sed_anomaly) / total_weight
  score[!support] <- NA_real_
  score[!is.finite(score)] <- NA_real_

  list(
    score = score,
    contrast = contrast,
    sed_anomaly = sed_anomaly,
    collapsed = collapsed,
    bands = bands,
    small_sigma = small_sigma,
    large_sigma = large_sigma,
    contrast_weight = contrast_weight,
    sed_weight = sed_weight
  )
}

#' Detect compact SED-clump seeds
#'
#' @param score 2-D clump-score map.
#' @param support Optional logical support mask.
#' @param n_clumps Optional exact number of seeds to keep. If `NULL`, seeds are
#'   selected above `score_quantile` up to `max_clumps`.
#' @param max_clumps Maximum number of seeds when `n_clumps = NULL`.
#' @param score_quantile Quantile threshold applied to finite support pixels.
#' @param min_distance Minimum Euclidean distance, in pixels, between seeds.
#'
#' @return A data frame with `clump_id`, `row`, `col`, and `score`.
#' @export
detect_sed_clumps <- function(score,
                              support = NULL,
                              n_clumps = NULL,
                              max_clumps = 50L,
                              score_quantile = 0.995,
                              min_distance = 3) {
  if (!is.matrix(score)) stop("`score` must be a matrix.")
  if (is.null(support)) support <- is.finite(score)
  support <- is.finite(support) & support & is.finite(score)
  if (!any(support, na.rm = TRUE)) {
    return(data.frame(clump_id = integer(), row = integer(), col = integer(), score = numeric()))
  }

  values <- score[support]
  if (!is.null(n_clumps)) {
    candidates <- which(support, arr.ind = TRUE)
  } else {
    score_quantile <- min(max(as.numeric(score_quantile), 0), 1)
    cut <- stats::quantile(values, probs = score_quantile, na.rm = TRUE, names = FALSE)
    candidates <- which(support & score >= cut, arr.ind = TRUE)
    if (!nrow(candidates)) {
      candidates <- which(support, arr.ind = TRUE)
    }
  }
  candidate_score <- score[candidates]
  ord <- order(candidate_score, decreasing = TRUE)
  candidates <- candidates[ord, , drop = FALSE]
  candidate_score <- candidate_score[ord]

  if (!is.null(n_clumps)) {
    max_clumps <- as.integer(n_clumps)
  }
  max_clumps <- max(1L, as.integer(max_clumps))
  min_distance <- max(0, as.numeric(min_distance))

  keep <- matrix(numeric(), ncol = 2)
  keep_score <- numeric()
  for (i in seq_len(nrow(candidates))) {
    point <- candidates[i, ]
    if (nrow(keep)) {
      d2 <- (keep[, 1] - point[1])^2 + (keep[, 2] - point[2])^2
      if (any(d2 < min_distance^2)) next
    }
    keep <- rbind(keep, point)
    keep_score <- c(keep_score, candidate_score[i])
    if (nrow(keep) >= max_clumps) break
  }

  if (!nrow(keep)) {
    return(data.frame(clump_id = integer(), row = integer(), col = integer(), score = numeric()))
  }

  data.frame(
    clump_id = seq_len(nrow(keep)),
    row = as.integer(keep[, 1]),
    col = as.integer(keep[, 2]),
    score = keep_score,
    stringsAsFactors = FALSE
  )
}

#' Grow clump seeds into compact footprints
#'
#' @param clump_seeds Data frame with `clump_id`, `row`, and `col`.
#' @param support Logical foreground support mask.
#' @param score Optional score map used to trim each local footprint.
#' @param grow_radius Initial dilation radius around each seed.
#' @param max_radius Maximum radius tried when a clump is too small.
#' @param min_pixels Minimum footprint size per clump.
#' @param footprint_quantile Optional local score quantile used to keep only
#'   the stronger part of each dilated footprint. Use `NULL` to keep all pixels.
#'
#' @return A list with `clump_labels` and `footprint_summary`.
#' @export
grow_clump_footprints <- function(clump_seeds,
                                  support,
                                  score = NULL,
                                  grow_radius = 3L,
                                  max_radius = 7L,
                                  min_pixels = 9L,
                                  footprint_quantile = 0.10) {
  if (!is.matrix(support)) stop("`support` must be a matrix.")
  support <- is.finite(support) & support
  labels <- matrix(NA_integer_, nrow = nrow(support), ncol = ncol(support))
  if (is.null(clump_seeds) || !nrow(clump_seeds)) {
    return(list(clump_labels = labels, footprint_summary = data.frame()))
  }
  if (!all(c("clump_id", "row", "col") %in% names(clump_seeds))) {
    stop("`clump_seeds` must contain clump_id, row, and col columns.")
  }

  ids <- sort(unique(clump_seeds$clump_id))
  sizes <- vapply(ids, function(id) sum(clump_seeds$clump_id == id), numeric(1))
  ids <- ids[order(sizes, decreasing = TRUE)]
  rows <- list()
  next_id <- 1L

  for (id in ids) {
    pix <- clump_seeds[clump_seeds$clump_id == id, , drop = FALSE]
    radius <- max(0L, as.integer(grow_radius))
    local <- matrix(FALSE, nrow = nrow(support), ncol = ncol(support))
    repeat {
      local <- .clump_dilate_points(pix, dim(support), radius = radius)
      local <- local & support & !is.finite(labels)
      if (!is.null(score) && !is.null(footprint_quantile) && any(local, na.rm = TRUE)) {
        q <- stats::quantile(score[local], probs = footprint_quantile, na.rm = TRUE, names = FALSE)
        if (is.finite(q)) {
          trimmed <- local & is.finite(score) & score >= q
          if (sum(trimmed, na.rm = TRUE) >= min_pixels || radius >= max_radius) {
            local <- trimmed
          }
        }
      }
      if (sum(local, na.rm = TRUE) >= min_pixels || radius >= max_radius) break
      radius <- radius + 1L
    }
    if (!any(local, na.rm = TRUE)) next

    labels[local] <- next_id
    rows[[length(rows) + 1L]] <- data.frame(
      input_clump_id = id,
      clump_id = next_id,
      n_pix = sum(local, na.rm = TRUE),
      grow_radius = radius,
      stringsAsFactors = FALSE
    )
    next_id <- next_id + 1L
  }

  list(
    clump_labels = labels,
    footprint_summary = if (length(rows)) do.call(rbind, rows) else data.frame()
  )
}

#' Split clump footprints into internal profile regions
#'
#' @param clump_labels Integer clump-footprint matrix.
#' @param score Score map used to rank pixels inside each clump.
#' @param clump_levels Number of internal profile regions per clump.
#'
#' @return A list with `profile_labels`, `clump_labels`, and `summary`.
#' @export
segment_clump_profiles <- function(clump_labels,
                                   score,
                                   clump_levels = 3L) {
  if (!is.matrix(clump_labels) || !is.matrix(score)) {
    stop("`clump_labels` and `score` must be matrices.")
  }
  if (!identical(dim(clump_labels), dim(score))) {
    stop("`clump_labels` and `score` must have the same dimensions.")
  }
  profile_labels <- matrix(NA_integer_, nrow = nrow(clump_labels), ncol = ncol(clump_labels))
  ids <- sort(unique(clump_labels[is.finite(clump_labels) & clump_labels > 0]))
  rows <- list()
  next_region <- 1L
  clump_levels <- max(1L, as.integer(clump_levels))

  for (id in ids) {
    idx <- which(clump_labels == id)
    if (!length(idx)) next
    values <- score[idx]
    n_level <- min(clump_levels, length(idx))
    local_level <- ceiling(rank(values, ties.method = "average") / length(values) * n_level)
    local_level <- pmin(n_level, pmax(1L, local_level))
    for (level in seq_len(n_level)) {
      pick <- idx[local_level == level]
      if (!length(pick)) next
      profile_labels[pick] <- next_region
      rows[[length(rows) + 1L]] <- data.frame(
        clump_id = id,
        profile_level = level,
        profile_region = next_region,
        n_pix = length(pick),
        mean_score = mean(score[pick], na.rm = TRUE),
        median_score = stats::median(score[pick], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      next_region <- next_region + 1L
    }
  }

  list(
    profile_labels = profile_labels,
    clump_labels = clump_labels,
    summary = if (length(rows)) do.call(rbind, rows) else data.frame()
  )
}

#' Segment compact SED clumps and the diffuse galaxy body
#'
#' @description
#' Experimental clump-aware segmentation mode. The method builds a loose
#' foreground support, detects compact SED-anomalous seeds (or accepts supplied
#' clump seeds), grows them into clump footprints, splits each clump into
#' internal profile levels, and runs [segment_regions_large()] on the remaining
#' diffuse body. The final `cluster_map` merges clump-profile regions and the
#' diffuse-body segmentation.
#'
#' @param input 3-D array or FITS-like list with `imDat`. This should usually
#'   be the PSF-matched cube used for segmentation.
#' @param Ncomp_body Number of diffuse-body regions.
#' @param bands Bands used for support, clump scoring, and body segmentation.
#' @param support Optional precomputed logical support mask.
#' @param support_info Optional support object, e.g. from
#'   [build_adaptive_support()]. If `NULL`, adaptive support is built.
#' @param support_args Named list passed to [build_adaptive_support()] when
#'   `support_info = NULL` and `support = NULL`.
#' @param clump_seeds Optional data frame with `clump_id`, `row`, and `col`.
#'   If `NULL`, seeds are detected automatically.
#' @param n_clumps Optional exact number of automatic clump seeds.
#' @param max_clumps Maximum number of automatic clump seeds.
#' @param score_quantile Quantile threshold for automatic clump seeds.
#' @param min_seed_distance Minimum seed separation in pixels.
#' @param clump_levels Number of internal profile regions per clump.
#' @param grow_radius Initial clump-footprint growth radius.
#' @param max_radius Maximum clump-footprint growth radius.
#' @param min_pixels Minimum pixels per grown clump footprint.
#' @param footprint_quantile Local score quantile used to trim grown footprints.
#' @param contrast_weight Weight of compact local contrast in the clump score.
#' @param sed_weight Weight of local SED anomaly in the clump score.
#' @param small_sigma Compact smoothing scale for the clump score.
#' @param large_sigma Background smoothing scale for the clump score.
#' @param body_segment Logical; segment the remaining diffuse body.
#' @param knn_k,auto_k,max_k,feature_scale,spatial_weight Sparse-Ward arguments
#'   passed to [segment_regions_large()] for the diffuse body.
#' @param cluster_pretransform Spectral pretransform for diffuse-body
#'   segmentation.
#' @param verbose Logical; print progress messages.
#'
#' @return A Sagui-like list with merged `cluster_map`, clump diagnostics, and
#'   body segmentation products.
#' @export
segment_clumps <- function(input,
                           Ncomp_body = 80,
                           bands = NULL,
                           support = NULL,
                           support_info = NULL,
                           support_args = list(),
                           clump_seeds = NULL,
                           n_clumps = NULL,
                           max_clumps = 50L,
                           score_quantile = 0.995,
                           min_seed_distance = 3,
                           clump_levels = 3L,
                           grow_radius = 3L,
                           max_radius = 7L,
                           min_pixels = 9L,
                           footprint_quantile = 0.10,
                           contrast_weight = 0.5,
                           sed_weight = 0.5,
                           small_sigma = 0.7,
                           large_sigma = 4,
                           body_segment = TRUE,
                           knn_k = 40,
                           auto_k = FALSE,
                           max_k = NULL,
                           feature_scale = c("none", "robust_col"),
                           spatial_weight = 0.10,
                           cluster_pretransform = "none",
                           verbose = TRUE) {
  feature_scale <- match.arg(feature_scale)
  cubedat <- if (is.list(input) && !is.null(input$imDat)) {
    input
  } else {
    list(imDat = input, hdr = NULL, axDat = NULL)
  }
  cube <- cubedat$imDat
  stopifnot(is.array(cube), length(dim(cube)) == 3L)
  bands <- .clump_resolve_bands(cube, bands)

  if (is.null(support_info) && is.null(support)) {
    defaults <- list(
      input = cubedat,
      bands = bands,
      transform = "asinh",
      z_threshold = 2,
      min_band_persistence = min(2L, length(bands)),
      single_band_z = 5,
      smooth_sigma = 0.5
    )
    support_info <- do.call(build_adaptive_support, utils::modifyList(defaults, support_args))
  }
  if (is.null(support)) {
    support <- support_info$mask
  }
  support <- is.finite(support) & support

  if (verbose) {
    message("Computing clump score...")
  }
  score_info <- compute_clump_score(
    input = cubedat,
    support = support,
    bands = bands,
    contrast_weight = contrast_weight,
    sed_weight = sed_weight,
    small_sigma = small_sigma,
    large_sigma = large_sigma
  )

  if (is.null(clump_seeds)) {
    if (verbose) message("Detecting SED-clump seeds...")
    clump_seeds <- detect_sed_clumps(
      score = score_info$score,
      support = support,
      n_clumps = n_clumps,
      max_clumps = max_clumps,
      score_quantile = score_quantile,
      min_distance = min_seed_distance
    )
  } else {
    clump_seeds <- .clump_normalize_seed_table(clump_seeds, dim(support))
  }
  n_seed_groups <- length(unique(clump_seeds$clump_id))

  if (verbose) {
    message(
      "Growing ", n_seed_groups, " clump seed group(s) from ",
      nrow(clump_seeds), " seed pixel(s) into footprints..."
    )
  }
  footprints <- grow_clump_footprints(
    clump_seeds = clump_seeds,
    support = support,
    score = score_info$score,
    grow_radius = grow_radius,
    max_radius = max_radius,
    min_pixels = min_pixels,
    footprint_quantile = footprint_quantile
  )
  profiles <- segment_clump_profiles(
    clump_labels = footprints$clump_labels,
    score = score_info$score,
    clump_levels = clump_levels
  )

  clump_mask <- is.finite(profiles$profile_labels) & profiles$profile_labels > 0
  body_support <- support & !clump_mask
  body_seg <- NULL
  merged <- profiles$profile_labels
  n_clump_regions <- length(unique(profiles$profile_labels[is.finite(profiles$profile_labels) & profiles$profile_labels > 0]))

  if (isTRUE(body_segment) && any(body_support, na.rm = TRUE) && Ncomp_body > 0) {
    if (verbose) message("Segmenting diffuse body with ", Ncomp_body, " requested region(s)...")
    body_input <- cubedat
    body_input$imDat <- mask_cube(cube, body_support, mode = "na")
    body_seg <- segment_regions_large(
      input = body_input,
      Ncomp = as.integer(Ncomp_body),
      use_starlet_mask = FALSE,
      mask_mode = "na",
      collapse_fn = function(x) .clump_collapse_bands(x, bands),
      scale_fn = identity,
      feature_scale = feature_scale,
      cluster_pretransform = cluster_pretransform,
      knn_k = knn_k,
      auto_k = auto_k,
      max_k = max_k,
      spatial_weight = spatial_weight,
      verbose = verbose
    )
    body_labels <- body_seg$cluster_map
    ok <- is.finite(body_labels) & body_labels > 0
    body_labels[ok] <- body_labels[ok] + n_clump_regions
    merged[ok] <- body_labels[ok]
  }

  out <- list(
    cluster_map = merged,
    mask = support,
    collapsed = score_info$collapsed,
    support = list(method = "clump", details = support_info),
    clump_score = score_info,
    clump_seeds = clump_seeds,
    clump_footprints = footprints,
    clump_profiles = profiles,
    body_segmentation = body_seg,
    original_cube = cubedat,
    backend = "clump_aware_sparse_ward",
    backend_info = list(
      Ncomp_body = Ncomp_body,
      n_clump_seeds = n_seed_groups,
      n_clump_seed_pixels = nrow(clump_seeds),
      n_clump_regions = n_clump_regions,
      clump_levels = clump_levels,
      bands = bands,
      feature_scale = feature_scale,
      spatial_weight = spatial_weight,
      cluster_pretransform = if (is.function(cluster_pretransform)) "custom" else cluster_pretransform
    )
  )
  class(out) <- c("sagui_clumps", "list")
  out
}

#' Plot clump-segmentation diagnostics
#'
#' @param x Result from `segment_clumps()`.
#' @param mode Plot mode: `"score"` shows the clump score, `"overlay"` overlays
#'   clump-profile regions on the score, `"profiles"` shows clump profiles only,
#'   and `"merged"` shows the final merged segmentation.
#' @param alpha Clump-profile opacity in overlay mode.
#'
#' @return A `ggplot2` object.
#' @export
plot_clump_diagnostics <- function(x,
                                   mode = c("overlay", "score", "profiles", "merged"),
                                   alpha = 0.72) {
  mode <- match.arg(mode)
  if (!is.list(x) || is.null(x$clump_score)) {
    stop("`x` must be a result returned by `segment_clumps()`.")
  }

  if (mode == "score") {
    return(.clump_plot_scalar(x$clump_score$score, palette = "Inferno"))
  }
  if (mode == "merged") {
    return(.clump_plot_labels(x$cluster_map, muted_body = FALSE))
  }
  if (mode == "profiles") {
    return(.clump_plot_labels(x$clump_profiles$profile_labels, muted_body = FALSE))
  }

  .clump_plot_overlay(
    scalar = x$clump_score$score,
    labels = x$clump_profiles$profile_labels,
    alpha = alpha
  )
}

.clump_cube <- function(input) {
  cube <- if (is.list(input) && !is.null(input$imDat)) input$imDat else input
  stopifnot(is.array(cube), length(dim(cube)) == 3L)
  cube
}

.clump_resolve_bands <- function(cube, bands) {
  n_wave <- dim(cube)[3]
  if (is.null(bands)) return(seq_len(n_wave))
  if (is.numeric(bands)) {
    bands <- as.integer(bands)
    bands <- bands[is.finite(bands) & bands >= 1L & bands <= n_wave]
    if (!length(bands)) stop("No valid band indices were supplied.")
    return(unique(bands))
  }
  if (is.character(bands)) {
    band_names <- dimnames(cube)[[3]]
    if (is.null(band_names)) stop("Character `bands` requires cube dimnames.")
    idx <- match(bands, band_names)
    if (any(is.na(idx))) stop("Some band names could not be matched.")
    return(unique(idx))
  }
  stop("`bands` must be NULL, numeric, or character.")
}

.clump_collapse_bands <- function(cube, bands) {
  Reduce(`+`, lapply(bands, function(j) cube[, , j])) / length(bands)
}

.clump_smooth_image <- function(image, sigma = 0) {
  image <- as.matrix(image)
  if (!is.finite(sigma) || sigma <= 0) return(image)
  if (requireNamespace("EBImage", quietly = TRUE)) {
    out <- EBImage::gblur(EBImage::Image(image), sigma = sigma)
    return(as.matrix(EBImage::imageData(out)))
  }
  if (requireNamespace("imager", quietly = TRUE)) {
    out <- imager::isoblur(imager::as.cimg(image), sigma = sigma)
    return(as.matrix(out[, , 1, 1]))
  }
  image
}

.clump_normalize01 <- function(x, mask = NULL, probs = c(0.01, 0.995)) {
  dims <- dim(x)
  if (is.null(mask)) mask <- is.finite(x)
  values <- x[mask & is.finite(x)]
  if (!length(values)) return(matrix(0, nrow = dims[1], ncol = dims[2]))
  lim <- stats::quantile(values, probs = probs, na.rm = TRUE, names = FALSE)
  y <- pmin(lim[2], pmax(lim[1], x))
  y <- y - min(y[mask & is.finite(y)], na.rm = TRUE)
  den <- max(y[mask & is.finite(y)], na.rm = TRUE)
  if (!is.finite(den) || den <= 0) den <- 1
  y <- y / den
  dim(y) <- dims
  y[!is.finite(y)] <- 0
  y
}

.clump_robust_scale_cube <- function(cube, support) {
  out <- cube
  for (j in seq_len(dim(cube)[3])) {
    band <- cube[, , j]
    values <- band[support & is.finite(band)]
    center <- stats::median(values, na.rm = TRUE)
    scale <- stats::mad(values, center = center, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) scale <- stats::sd(values, na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) scale <- 1
    out[, , j] <- (band - center) / scale
  }
  out[!is.finite(out)] <- 0
  out
}

.clump_dilate_points <- function(points, dims, radius = 0L) {
  mask <- matrix(FALSE, nrow = dims[1], ncol = dims[2])
  radius <- max(0L, as.integer(radius))
  for (i in seq_len(nrow(points))) {
    rr <- points$row[i]
    cc <- points$col[i]
    for (dr in seq.int(-radius, radius)) {
      for (dc in seq.int(-radius, radius)) {
        if (dr * dr + dc * dc > radius * radius) next
        r <- rr + dr
        c <- cc + dc
        if (r < 1L || r > dims[1] || c < 1L || c > dims[2]) next
        mask[r, c] <- TRUE
      }
    }
  }
  mask
}

.clump_normalize_seed_table <- function(seeds, dims) {
  if (!is.data.frame(seeds)) seeds <- as.data.frame(seeds)
  if (!("row" %in% names(seeds)) && "x" %in% names(seeds)) seeds$row <- seeds$x
  if (!("col" %in% names(seeds)) && "y" %in% names(seeds)) seeds$col <- seeds$y
  if (!("clump_id" %in% names(seeds))) seeds$clump_id <- seq_len(nrow(seeds))
  if (!all(c("clump_id", "row", "col") %in% names(seeds))) {
    stop("`clump_seeds` must contain row/col or x/y coordinates.")
  }
  add_one <- any(seeds$row == 0 | seeds$col == 0, na.rm = TRUE)
  out <- data.frame(
    clump_id = seeds$clump_id,
    row = as.integer(seeds$row) + if (add_one) 1L else 0L,
    col = as.integer(seeds$col) + if (add_one) 1L else 0L,
    stringsAsFactors = FALSE
  )
  ok <- out$row >= 1L & out$row <= dims[1] & out$col >= 1L & out$col <= dims[2]
  out[ok, , drop = FALSE]
}

.clump_melt <- function(mat, value = "value") {
  df <- expand.grid(
    Row = seq_len(nrow(mat)),
    Col = seq_len(ncol(mat))
  )
  df[[value]] <- as.vector(mat)
  df
}

.clump_plot_scalar <- function(mat, palette = "Inferno") {
  df <- .clump_melt(.clump_normalize01(mat), "value")
  ggplot2::ggplot(df, ggplot2::aes(Row, Col, fill = value)) +
    ggplot2::geom_tile(width = 1, height = 1) +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::scale_fill_gradientn(colours = grDevices::hcl.colors(256, palette), guide = "none") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "black", colour = NA),
      panel.background = ggplot2::element_rect(fill = "black", colour = NA)
    )
}

.clump_plot_labels <- function(labels, muted_body = FALSE) {
  ids <- sort(unique(labels[is.finite(labels) & labels > 0]))
  cols <- grDevices::colorRampPalette(c("#213E60", "#E68C3A", "#66A7B0", "#D65F59", "#E7C55D", "#315C91"))(max(2L, length(ids)))
  names(cols) <- as.character(ids)
  df <- .clump_melt(labels, "region")
  df$region <- factor(df$region, levels = ids)
  ggplot2::ggplot(df, ggplot2::aes(Row, Col, fill = region)) +
    ggplot2::geom_tile(width = 1, height = 1) +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::scale_fill_manual(values = cols, na.value = "white", guide = "none", drop = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

.clump_plot_overlay <- function(scalar, labels, alpha = 0.72) {
  bg_norm <- .clump_normalize01(scalar)
  bg_idx <- matrix(
    as.integer(cut(bg_norm, breaks = seq(0, 1, length.out = 257), include.lowest = TRUE)),
    nrow = nrow(bg_norm),
    ncol = ncol(bg_norm)
  )
  bg_col <- matrix(
    grDevices::hcl.colors(256, "Inferno")[bg_idx],
    nrow = nrow(bg_norm),
    ncol = ncol(bg_norm)
  )
  bg <- .clump_melt(bg_col, "fill_col")
  lab <- .clump_melt(labels, "region")
  ids <- sort(unique(labels[is.finite(labels) & labels > 0]))
  cols <- grDevices::colorRampPalette(c("#213E60", "#E68C3A", "#66A7B0", "#D65F59", "#E7C55D", "#315C91"))(max(2L, length(ids)))
  names(cols) <- as.character(ids)
  lab$fill_col <- cols[as.character(lab$region)]
  lab <- lab[is.finite(lab$region) & lab$region > 0, , drop = FALSE]

  ggplot2::ggplot() +
    ggplot2::geom_tile(data = bg, ggplot2::aes(Row, Col, fill = fill_col), width = 1, height = 1) +
    ggplot2::geom_tile(
      data = lab,
      ggplot2::aes(Row, Col, fill = fill_col),
      width = 1,
      height = 1,
      alpha = alpha
    ) +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::scale_fill_identity(guide = "none") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "black", colour = NA),
      panel.background = ggplot2::element_rect(fill = "black", colour = NA)
    )
}
