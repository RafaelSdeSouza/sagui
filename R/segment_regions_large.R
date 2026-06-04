.sparse_ward_scale_features <- function(X, scale_fn = NULL) {
  X <- as.matrix(X)

  scaler <- if (is.null(scale_fn)) .safe_scale_spectrum else scale_fn

  scaled <- t(vapply(
    seq_len(nrow(X)),
    function(i) as.numeric(scaler(X[i, ])),
    numeric(ncol(X))
  ))
  scaled[!is.finite(scaled)] <- 0
  scaled
}

.sparse_ward_robust_col_scale <- function(X) {
  X <- as.matrix(X)

  center <- apply(X, 2, stats::median, na.rm = TRUE)
  scale <- apply(X, 2, stats::mad, center = center, constant = 1.4826, na.rm = TRUE)
  fallback <- stats::median(scale[is.finite(scale) & scale > 0], na.rm = TRUE)
  if (!is.finite(fallback) || fallback <= 0) {
    fallback <- 1
  }

  scale[!is.finite(scale) | scale <= 0] <- fallback
  X <- sweep(X, 2, center, "-")
  X <- sweep(X, 2, scale, "/")
  X[!is.finite(X)] <- 0
  X
}

.sparse_ward_cluster_matrix <- function(features,
                                        Ncomp,
                                        knn_k = 40,
                                        auto_k = FALSE,
                                        max_k = NULL,
                                        verbose = FALSE) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required for the sparse-Ward backend.")
  }

  features <- as.matrix(features)
  features[!is.finite(features)] <- 0
  n <- nrow(features)

  if (n < 1L) {
    stop("`features` has no rows.")
  }
  if (!is.finite(Ncomp) || Ncomp < 1L || Ncomp > n) {
    stop("`Ncomp` must be between 1 and the number of feature rows.")
  }
  if (Ncomp == n) {
    return(list(
      labels = seq_len(n),
      knn_k = 0L,
      requested_Ncomp = Ncomp,
      actual_Ncomp = n,
      disconnected = FALSE
    ))
  }

  knn_k <- as.integer(knn_k)
  if (!is.finite(knn_k) || knn_k < 2L) {
    stop("`knn_k` must be an integer >= 2.")
  }
  knn_k <- min(knn_k, n - 1L)

  if (is.null(max_k)) {
    max_k <- min(n - 1L, max(knn_k, 320L))
  } else {
    max_k <- min(as.integer(max_k), n - 1L)
  }
  if (max_k < knn_k) {
    max_k <- knn_k
  }

  labels <- NULL
  actual_k <- NA_integer_
  used_k <- knn_k

  repeat {
    if (verbose) {
      message(sprintf("Building sparse-Ward kNN graph with k = %d...", used_k))
    }

    nn <- RANN::nn2(features, k = used_k + 1L)$nn.idx
    nn <- nn[, -1L, drop = FALSE]

    labels <- sagui_sparse_ward_cut_cpp(features, nn, as.integer(Ncomp))
    actual_k <- length(unique(labels))

    if (verbose) {
      message(sprintf(
        "Sparse Ward returned %d clusters for requested Ncomp = %d.",
        actual_k,
        Ncomp
      ))
    }

    if (actual_k == Ncomp || !isTRUE(auto_k) || used_k >= max_k) {
      break
    }
    used_k <- min(max_k, max(used_k + 1L, ceiling(2 * used_k)))
  }

  disconnected <- actual_k > Ncomp
  if (disconnected) {
    warning(
      "Sparse-Ward graph remained disconnected: requested ",
      Ncomp,
      " clusters but obtained ",
      actual_k,
      ". Try increasing `knn_k`, increasing `max_k`, or adding a small `spatial_weight`."
    )
  }

  list(
    labels = as.integer(labels),
    knn_k = used_k,
    requested_Ncomp = as.integer(Ncomp),
    actual_Ncomp = as.integer(actual_k),
    disconnected = disconnected
  )
}

#' Segment photometric regions from large spectral cubes
#'
#' Large-cube version of [segment_regions()]. It keeps the same photometric
#' mask, pretransform, and output structure, but replaces exact all-pairs Ward
#' clustering with a kNN-restricted sparse-Ward approximation.
#'
#' @inheritParams segment_regions
#' @param knn_k Number of nearest neighbours for the sparse-Ward graph.
#' @param auto_k Logical; if `TRUE`, increase `knn_k` when the graph remains
#'   disconnected before the requested number of regions is reached.
#' @param max_k Maximum `k` allowed when `auto_k = TRUE`.
#' @param feature_scale Optional column-wise feature scaling after row scaling.
#' @param spatial_weight Optional weight for appending normalized x/y pixel
#'   coordinates to the clustering features.
#' @param return_details Logical; include valid indices, features, and labels.
#' @param verbose Logical; print sparse-Ward progress messages.
#'
#' @return A segmentation result list compatible with [segment_regions()] and
#'   downstream photometry helpers such as [RegionPhotometry()].
#' @export
segment_regions_large <- function(input,
                                  Ncomp = 5,
                                  redshift = 0,
                                  pretransform = NULL,
                                  mask_pretransform = "none",
                                  cluster_pretransform = "none",
                                  scale_fn = median_scale,
                                  n_regions = NULL,
                                  use_starlet_mask = TRUE,
                                  support_method = c("starlet", "adaptive", "starlet_contourlet", "contourlet"),
                                  support_args = list(),
                                  collapse_fn = collapse_white_light,
                                  starlet_J = 5,
                                  starlet_scales = 2:5,
                                  include_coarse = FALSE,
                                  denoise_k = 2.5,
                                  mode = c("soft", "hard"),
                                  positive_only = TRUE,
                                  clean_mask = FALSE,
                                  min_mask_area = 1L,
                                  close_size = 1L,
                                  open_size = 1L,
                                  keep_largest = FALSE,
                                  mask_mode = c("na", "zero"),
                                  hclust_method = "ward.D2",
                                  knn_k = 40,
                                  auto_k = FALSE,
                                  max_k = NULL,
                                  feature_scale = c("none", "robust_col"),
                                  spatial_weight = 0,
                                  return_details = FALSE,
                                  verbose = TRUE) {
  mode <- match.arg(mode)
  mask_mode <- match.arg(mask_mode)
  support_method <- match.arg(support_method)
  feature_scale <- match.arg(feature_scale)

  if (!is.null(n_regions)) {
    Ncomp <- n_regions
  }
  if (!is.null(pretransform)) {
    cluster_pretransform <- pretransform
  }
  if (!identical(hclust_method, "ward.D2") && verbose) {
    message("`hclust_method` is ignored by `segment_regions_large()`; sparse Ward is always used.")
  }

  cubedat <- if (is.list(input) && !is.null(input$imDat)) {
    input
  } else {
    list(imDat = input, hdr = NULL, axDat = NULL)
  }
  cube <- cubedat$imDat
  stopifnot(is.array(cube), length(dim(cube)) == 3L)

  n_row <- dim(cube)[1]
  n_col <- dim(cube)[2]
  n_wave <- dim(cube)[3]

  if (isTRUE(use_starlet_mask)) {
    if (support_method == "starlet") {
      mask_info <- build_starlet_mask(
        input = cubedat,
        collapse_fn = collapse_fn,
        pretransform = mask_pretransform,
        starlet_J = starlet_J,
        starlet_scales = starlet_scales,
        include_coarse = include_coarse,
        denoise_k = denoise_k,
        mode = mode,
        positive_only = positive_only,
        clean_mask = clean_mask,
        min_mask_area = min_mask_area,
        close_size = close_size,
        open_size = open_size,
        keep_largest = keep_largest
      )
    } else if (support_method == "adaptive") {
      adaptive_defaults <- list(
        input = cubedat,
        pretransform = mask_pretransform
      )
      if (!is.list(support_args)) {
        stop("`support_args` must be a named list.")
      }
      mask_info <- do.call(
        build_adaptive_support,
        utils::modifyList(adaptive_defaults, support_args)
      )
      mask_info <- .apply_mask_cleanup(
        mask_info,
        clean_mask = clean_mask,
        min_mask_area = min_mask_area,
        close_size = close_size,
        open_size = open_size,
        keep_largest = keep_largest
      )
    } else {
      contourlet_defaults <- list(
        input = cubedat,
        collapse_fn = collapse_fn,
        pretransform = mask_pretransform,
        combine = if (support_method == "starlet_contourlet") "hybrid" else "standalone",
        starlet_J = starlet_J,
        starlet_scales = starlet_scales,
        include_coarse = include_coarse,
        denoise_k = denoise_k,
        mode = mode,
        positive_only = positive_only
      )
      if (!is.list(support_args)) {
        stop("`support_args` must be a named list.")
      }
      mask_info <- do.call(
        build_contourlet_mask,
        utils::modifyList(contourlet_defaults, support_args)
      )
      mask_info <- .apply_mask_cleanup(
        mask_info,
        clean_mask = clean_mask,
        min_mask_area = min_mask_area,
        close_size = close_size,
        open_size = open_size,
        keep_largest = keep_largest
      )
    }
    spatial_mask <- mask_info$mask
  } else {
    collapsed <- collapse_fn(cube)
    spatial_mask <- is.finite(collapsed)
    mask_info <- list(
      collapsed = collapsed,
      decomposition = NULL,
      reconstruction = collapsed,
      mask = spatial_mask
    )
  }

  masked_cube <- mask_cube(cube, spatial_mask, mode = mask_mode)
  masked_input <- cubedat
  masked_input$imDat <- masked_cube
  IFU2D <- cube_to_matrix(masked_input)

  signal <- rowSums(IFU2D, na.rm = TRUE)
  noise <- sqrt(signal)
  signal[is.na(signal) | signal <= 0] <- 0
  noise[is.na(noise) | noise == 0] <- Inf

  finite_counts <- rowSums(is.finite(IFU2D))
  finite_frac <- finite_counts / n_wave
  row_energy <- rowSums(IFU2D^2, na.rm = TRUE)
  row_mad <- apply(IFU2D, 1, function(v) stats::mad(v[is.finite(v)], na.rm = TRUE))

  valid <- signal > 0 &
    finite_counts >= pmin(10L, n_wave) &
    finite_frac >= 0.8 &
    is.finite(row_energy) & row_energy > 0 &
    is.finite(row_mad) & row_mad > 0

  valid_indices <- which(valid)
  if (!length(valid_indices)) {
    stop("No valid spectra remain after masking.")
  }
  if (length(valid_indices) < Ncomp) {
    stop("`Ncomp` is larger than the number of valid spectra.")
  }

  IFU2D_valid <- IFU2D[valid_indices, , drop = FALSE]
  signal_valid <- signal[valid_indices]
  noise_valid <- noise[valid_indices]

  transformed_data <- pretransform_spectra(IFU2D_valid, method = cluster_pretransform)
  features <- .sparse_ward_scale_features(transformed_data, scale_fn = scale_fn)

  if (feature_scale == "robust_col") {
    features <- .sparse_ward_robust_col_scale(features)
  }

  if (spatial_weight > 0) {
    ij <- arrayInd(valid_indices, .dim = c(n_row, n_col))
    x_sd <- stats::sd(ij[, 2])
    y_sd <- stats::sd(ij[, 1])
    if (!is.finite(x_sd) || x_sd <= 0) x_sd <- 1
    if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 1
    xy <- cbind(
      x = (ij[, 2] - mean(ij[, 2])) / x_sd,
      y = (ij[, 1] - mean(ij[, 1])) / y_sd
    )
    features <- cbind(features, spatial_weight * xy)
  }

  if (verbose) {
    message(sprintf("Valid pixels for sparse Ward: %d", length(valid_indices)))
  }
  fit <- .sparse_ward_cluster_matrix(
    features = features,
    Ncomp = Ncomp,
    knn_k = knn_k,
    auto_k = auto_k,
    max_k = max_k,
    verbose = verbose
  )
  clusters <- fit$labels

  cluster_map <- matrix(NA_integer_, nrow = n_row, ncol = n_col)
  cluster_map[valid_indices] <- clusters

  cluster_snr <- vapply(sort(unique(clusters)), function(cluster_id) {
    idx <- which(clusters == cluster_id)
    sum(signal_valid[idx]) / sqrt(sum(noise_valid[idx]^2))
  }, numeric(1))

  out <- list(
    cluster_map = cluster_map,
    mask = spatial_mask,
    collapsed = mask_info$collapsed,
    starlet = list(
      decomposition = if (!is.null(mask_info$decomposition)) {
        mask_info$decomposition
      } else if (!is.null(mask_info$base_support)) {
        mask_info$base_support$decomposition
      } else {
        NULL
      },
      reconstruction = if (!is.null(mask_info$reconstruction)) {
        mask_info$reconstruction
      } else if (!is.null(mask_info$base_support)) {
        mask_info$base_support$reconstruction
      } else {
        NULL
      }
    ),
    support = list(
      method = if (isTRUE(use_starlet_mask)) support_method else "none",
      details = mask_info
    ),
    masked_cube = masked_cube,
    cluster_snr = cluster_snr,
    header = cubedat$hdr,
    axDat = cubedat$axDat,
    original_cube = cubedat,
    redshift = redshift,
    mask_pretransform = if (is.function(mask_pretransform)) "custom" else mask_pretransform,
    support_method = if (isTRUE(use_starlet_mask)) support_method else "none",
    cluster_pretransform = if (is.function(cluster_pretransform)) "custom" else cluster_pretransform,
    backend = "sparse_ward",
    backend_info = list(
      algorithm = "sparse_ward",
      knn_k = fit$knn_k,
      requested_Ncomp = fit$requested_Ncomp,
      actual_Ncomp = fit$actual_Ncomp,
      disconnected = fit$disconnected,
      scale_fn = if (is.null(scale_fn)) "none" else "custom",
      feature_scale = feature_scale,
      spatial_weight = spatial_weight,
      valid_pixels = length(valid_indices)
    )
  )

  if (return_details) {
    out <- c(out, list(
      valid_indices = valid_indices,
      features = features,
      labels = clusters
    ))
  }

  out
}
