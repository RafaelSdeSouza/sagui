#' Build a starlet-based segmentation mask from a spectral cube
#'
#' @param input 3-D array or FITS-like list with `imDat`.
#' @param collapse_fn Function used to collapse the cube to a 2-D image.
#' @param pretransform Optional spectral pretransform applied to the cube before
#'   the white-light collapse and starlet decomposition. This is useful when the
#'   photometric mask should be derived from a transformed representation rather
#'   than the original flux cube.
#' @param starlet_J Number of starlet scales.
#' @param starlet_scales Scales to reconstruct.
#' @param include_coarse Logical; include the coarse plane in the reconstruction.
#' @param denoise_k Optional starlet denoising threshold.
#' @param mode Thresholding mode for starlet reconstruction.
#' @param positive_only Logical; keep only positive reconstruction values in the mask.
#' @return A list with the collapsed image, decomposition, reconstruction, and logical mask.
#' @export
build_starlet_mask <- function(input,
                               collapse_fn = collapse_white_light,
                               pretransform = "none",
                               starlet_J = 5,
                               starlet_scales = 2:5,
                               include_coarse = FALSE,
                               denoise_k = 0,
                               mode = c("soft", "hard"),
                               positive_only = TRUE) {
  mode <- match.arg(mode)
  cube <- if (is.list(input) && !is.null(input$imDat)) input$imDat else input
  stopifnot(is.array(cube), length(dim(cube)) == 3L)

  transformed_cube <- pretransform_cube(cube, method = pretransform)
  collapsed <- collapse_fn(transformed_cube)
  decomposition <- starlet_mask(collapsed, J = starlet_J)
  reconstruction <- starlet_reconstruct(
    decomposition,
    keep_scales = starlet_scales,
    include_coarse = include_coarse,
    denoise_k = denoise_k,
    mode = mode
  )

  mask <- is.finite(reconstruction)
  if (isTRUE(positive_only)) {
    mask <- mask & reconstruction > 0
  }

  list(
    collapsed = collapsed,
    decomposition = decomposition,
    reconstruction = reconstruction,
    mask = mask
  )
}

# Robust row scaling for spectral vectors before clustering.
.safe_scale_spectrum <- function(x) {
  xf <- x[is.finite(x)]
  if (!length(xf)) {
    return(rep(0, length(x)))
  }

  center <- stats::median(xf)
  scale <- stats::mad(xf, center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    scale <- stats::sd(xf, na.rm = TRUE)
  }
  if (!is.finite(scale) || scale <= 0) {
    scale <- 1
  }

  out <- (x - center) / scale
  out[!is.finite(out)] <- 0
  out
}

.compute_distance_matrix <- function(x) {
  if (requireNamespace("torch", quietly = TRUE)) {
    xt <- torch::torch_tensor(x, dtype = torch::torch_float())
    dmat <- as.matrix(torch::as_array(torch::torch_cdist(xt, xt, p = 2)$cpu()))
    return(stats::as.dist(dmat))
  }

  stats::dist(x)
}

.compute_hclust <- function(d, method) {
  if (requireNamespace("fastcluster", quietly = TRUE)) {
    return(fastcluster::hclust(d, method = method))
  }

  stats::hclust(d, method = method)
}

#' Segment photometric regions from a spectral cube
#'
#' @param input 3-D array or FITS-like list with `imDat`.
#' @param Ncomp Number of output regions.
#' @param redshift Reserved compatibility argument carried over from `capivara`.
#'   Currently unused.
#' @param pretransform Deprecated alias for `cluster_pretransform`.
#' @param mask_pretransform Optional spectral pretransform applied before the
#'   white-light collapse and starlet mask. This lets the segmentation mask be
#'   built on a transformed photometric representation while leaving the
#'   clustering stage untouched.
#' @param cluster_pretransform Optional spectral pretransform applied to the
#'   valid spectra matrix before row-wise scaling and clustering. May be one of
#'   `"none"`, `"asinh"`, `"log1p"`, `"signed_log1p"`, `"copula_uniform"`,
#'   `"copula_gaussian"`, or a custom function returning a matrix with the same
#'   dimensions as the input.
#' @param scale_fn Per-spectrum scaling function applied row-wise before
#'   clustering.
#' @param n_regions Deprecated alias for `Ncomp`.
#' @param use_starlet_mask Logical; if `TRUE`, derive a photometric mask before clustering.
#' @param collapse_fn Function used to collapse the cube to a 2-D image.
#' @param starlet_J Number of starlet scales.
#' @param starlet_scales Scales to keep when reconstructing the starlet image.
#' @param include_coarse Logical; include the coarse plane in the starlet reconstruction.
#' @param denoise_k Optional starlet denoising threshold.
#' @param mode Starlet thresholding mode.
#' @param positive_only Logical; keep only positive reconstruction values in the mask.
#' @param mask_mode Mask fill mode passed to `mask_cube()`.
#' @param hclust_method Linkage method passed to `hclust()`.
#' @return A segmentation result list containing the cluster map, mask products, and metadata.
#' @export
segment_regions <- function(input,
                            Ncomp = 5,
                            redshift = 0,
                            pretransform = NULL,
                            mask_pretransform = "none",
                            cluster_pretransform = "none",
                            scale_fn = median_scale,
                            n_regions = NULL,
                            use_starlet_mask = TRUE,
                            collapse_fn = collapse_white_light,
                            starlet_J = 5,
                            starlet_scales = 2:5,
                            include_coarse = FALSE,
                            denoise_k = 0,
                            mode = c("soft", "hard"),
                            positive_only = TRUE,
                            mask_mode = c("na", "zero"),
                            hclust_method = "ward.D2") {
  mode <- match.arg(mode)
  mask_mode <- match.arg(mask_mode)
  if (!is.null(n_regions)) {
    Ncomp <- n_regions
  }
  if (!is.null(pretransform)) {
    cluster_pretransform <- pretransform
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
  spectra <- cube_to_matrix(cubedat)

  if (isTRUE(use_starlet_mask)) {
    mask_info <- build_starlet_mask(
      input = cubedat,
      collapse_fn = collapse_fn,
      pretransform = mask_pretransform,
      starlet_J = starlet_J,
      starlet_scales = starlet_scales,
      include_coarse = include_coarse,
      denoise_k = denoise_k,
      mode = mode,
      positive_only = positive_only
    )
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
  scaler <- if (is.null(scale_fn)) .safe_scale_spectrum else scale_fn
  scaled_data <- t(apply(transformed_data, 1, scaler))
  scaled_data[!is.finite(scaled_data)] <- 0

  distance_matrix <- .compute_distance_matrix(scaled_data)
  hc <- .compute_hclust(distance_matrix, method = hclust_method)
  clusters <- stats::cutree(hc, k = Ncomp)

  cluster_map <- matrix(NA_integer_, nrow = n_row, ncol = n_col)
  cluster_map[valid_indices] <- clusters

  cluster_snr <- vapply(sort(unique(clusters)), function(cluster_id) {
    idx <- which(clusters == cluster_id)
    sum(signal_valid[idx]) / sqrt(sum(noise_valid[idx]^2))
  }, numeric(1))

  list(
    cluster_map = cluster_map,
    mask = spatial_mask,
    collapsed = mask_info$collapsed,
    starlet = list(
      decomposition = mask_info$decomposition,
      reconstruction = mask_info$reconstruction
    ),
    masked_cube = masked_cube,
    cluster_snr = cluster_snr,
    header = cubedat$hdr,
    axDat = cubedat$axDat,
    original_cube = cubedat,
    redshift = redshift,
    mask_pretransform = if (is.function(mask_pretransform)) "custom" else mask_pretransform,
    cluster_pretransform = if (is.function(cluster_pretransform)) "custom" else cluster_pretransform,
    hclust = hc
  )
}
