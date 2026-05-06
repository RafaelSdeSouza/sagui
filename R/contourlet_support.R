#' Build a contourlet-like directional support mask from a spectral cube
#'
#' This function extends the default starlet support with a directional,
#' multiscale arm-sensitive support map. The implementation is contourlet-like
#' rather than a full contourlet transform: it uses a bank of oriented,
#' elongated Gaussian filters across multiple scales, combines them into a
#' directional response map, and optionally grows the starlet support into
#' nearby arm-like pixels.
#'
#' @param input 3-D array or FITS-like list with `imDat`.
#' @param collapse_fn Function used to collapse the cube to a 2-D image.
#' @param pretransform Optional spectral pretransform applied to the cube before
#'   the white-light collapse.
#' @param combine Either `"hybrid"` (recommended; start from the starlet support
#'   and extend it directionally) or `"standalone"` (use only the directional
#'   support).
#' @param base_mask Optional precomputed base mask. When `combine = "hybrid"`
#'   and `base_mask` is `NULL`, the base support is built from
#'   [build_starlet_mask()].
#' @param starlet_J Number of starlet scales used when building the hybrid base
#'   support.
#' @param starlet_scales Scales to reconstruct for the hybrid base support.
#' @param include_coarse Logical; include the coarse plane in the hybrid base
#'   support reconstruction.
#' @param denoise_k Optional starlet denoising threshold used in the hybrid base
#'   support.
#' @param mode Thresholding mode for the hybrid base support reconstruction.
#' @param positive_only Logical; keep only positive reconstruction values in the
#'   hybrid base support.
#' @param angles_deg Vector of filter orientations in degrees.
#' @param scales List of `(sigma_major, sigma_minor)` pairs defining the
#'   directional filter bank.
#' @param background_sigma Smoothing scale used to estimate the local background
#'   before the directional response is measured.
#' @param score_low_q Lower quantile used to define the allowed directional
#'   growth region.
#' @param score_high_q Upper quantile used to define strong directional seeds.
#' @param contrast_low_q Lower quantile used to define the allowed contrast
#'   growth region.
#' @param contrast_high_q Upper quantile used to define strong contrast seeds.
#' @param max_growth_iter Maximum number of constrained growth iterations.
#' @param close_size Size of the binary closing brush.
#' @param open_size Size of the binary opening brush.
#' @return A list containing the collapsed image, the optional base support, the
#'   final logical mask, the directional score map, the local contrast map, and
#'   threshold metadata.
#' @export
build_contourlet_mask <- function(input,
                                  collapse_fn = collapse_white_light,
                                  pretransform = "none",
                                  combine = c("hybrid", "standalone"),
                                  base_mask = NULL,
                                  starlet_J = 5,
                                  starlet_scales = 2:5,
                                  include_coarse = FALSE,
                                  denoise_k = 2.5,
                                  mode = c("soft", "hard"),
                                  positive_only = TRUE,
                                  angles_deg = seq(0, 165, by = 15),
                                  scales = list(c(4.5, 1.2), c(7.5, 1.8), c(11.0, 2.6)),
                                  background_sigma = 8,
                                  score_low_q = 0.92,
                                  score_high_q = 0.98,
                                  contrast_low_q = 0.80,
                                  contrast_high_q = 0.93,
                                  max_growth_iter = 8L,
                                  close_size = 3L,
                                  open_size = 3L) {
  combine <- match.arg(combine)
  mode <- match.arg(mode)

  cube <- if (is.list(input) && !is.null(input$imDat)) input$imDat else input
  stopifnot(is.array(cube), length(dim(cube)) == 3L)

  transformed_cube <- pretransform_cube(cube, method = pretransform)
  collapsed <- collapse_fn(transformed_cube)

  nr <- nrow(collapsed)
  nc <- ncol(collapsed)

  if (!is.null(base_mask) && !identical(dim(base_mask), c(nr, nc))) {
    stop("`base_mask` must have the same spatial dimensions as the collapsed image.")
  }

  base_support <- NULL
  if (combine == "hybrid") {
    if (is.null(base_mask)) {
      base_support <- build_starlet_mask(
        input = input,
        collapse_fn = collapse_fn,
        pretransform = pretransform,
        starlet_J = starlet_J,
        starlet_scales = starlet_scales,
        include_coarse = include_coarse,
        denoise_k = denoise_k,
        mode = mode,
        positive_only = positive_only
      )
      seed_mask <- base_support$mask
    } else {
      seed_mask <- as.matrix(base_mask)
    }
  } else {
    seed_mask <- matrix(FALSE, nrow = nr, ncol = nc)
  }

  white_light_prepped <- asinh(pmax(collapsed, 0))
  background_map <- .directional_gaussian_blur(white_light_prepped, sigma = background_sigma)
  contrast_map <- pmax(white_light_prepped - background_map, 0)
  arm_score <- .contourlet_like_response(contrast_map, angles_deg = angles_deg, scales = scales)

  seed_dilated <- .binary_dilate(seed_mask, size = 5L)
  outside_seed <- !seed_dilated & is.finite(arm_score)
  if (!any(outside_seed, na.rm = TRUE)) {
    outside_seed <- is.finite(arm_score)
  }

  score_low <- .positive_quantile(arm_score[outside_seed], score_low_q, fallback_k = 1.6)
  score_high <- .positive_quantile(arm_score[outside_seed], score_high_q, fallback_k = 2.5)
  contrast_low <- .positive_quantile(contrast_map[outside_seed], contrast_low_q, fallback_k = 1.2)
  contrast_high <- .positive_quantile(contrast_map[outside_seed], contrast_high_q, fallback_k = 1.8)

  strong_arm <- (arm_score > score_high) & (contrast_map > contrast_high)
  directional_candidate <- (arm_score > score_low) & (contrast_map > contrast_low)

  if (combine == "hybrid") {
    allowed_growth <- seed_mask | directional_candidate
    seed_for_growth <- seed_mask | strong_arm
  } else {
    allowed_growth <- directional_candidate
    seed_for_growth <- strong_arm
    if (!any(seed_for_growth, na.rm = TRUE)) {
      seed_for_growth <- directional_candidate
    }
  }

  allowed_growth <- .binary_close_open(
    allowed_growth,
    close_size = close_size,
    open_size = open_size
  )

  grown <- .grow_mask_constrained(
    seed = seed_for_growth,
    allowed = allowed_growth,
    max_iter = max_growth_iter
  )
  grown <- .binary_close_open(
    grown,
    close_size = close_size,
    open_size = open_size
  )

  list(
    collapsed = collapsed,
    base_support = base_support,
    mask = grown,
    arm_score = arm_score,
    contrast_map = contrast_map,
    added_pixels = grown & !seed_mask,
    score_low = score_low,
    score_high = score_high,
    contrast_low = contrast_low,
    contrast_high = contrast_high,
    combine = combine
  )
}

.directional_support_backend <- function() {
  if (requireNamespace("EBImage", quietly = TRUE)) {
    return("EBImage")
  }
  if (requireNamespace("imager", quietly = TRUE)) {
    return("imager")
  }
  stop(
    "Directional contourlet-like support requires the optional package ",
    "`EBImage` or `imager`."
  )
}

.image_from_matrix <- function(x, backend) {
  if (backend == "EBImage") {
    return(EBImage::Image(x))
  }

  imager::as.cimg(t(x), x = ncol(x), y = nrow(x), cc = 1, z = 1)
}

.matrix_from_image <- function(x, backend) {
  arr <- as.array(x)
  if (backend == "EBImage") {
    if (length(dim(arr)) <= 2L) {
      return(arr)
    }
    return(arr[, , 1])
  }

  t(as.matrix(x))
}

.disc_brush <- function(size = 3L) {
  size <- as.integer(size)
  if (!is.finite(size) || size < 1L) {
    size <- 1L
  }
  if (size %% 2L == 0L) {
    size <- size + 1L
  }

  radius <- (size - 1L) / 2
  grid <- seq(-radius, radius)
  xy <- expand.grid(x = grid, y = grid)
  mat <- matrix((xy$x^2 + xy$y^2) <= (radius^2 + 1e-8), nrow = size, ncol = size)
  mat * 1
}

.binary_dilate <- function(mask, size = 3L) {
  backend <- .directional_support_backend()
  brush <- .disc_brush(size)

  if (backend == "EBImage") {
    out <- EBImage::dilate(EBImage::Image(mask * 1), brush)
    return(.matrix_from_image(out, backend) > 0)
  }

  out <- imager::dilate(
    .image_from_matrix(mask * 1, backend),
    .image_from_matrix(brush, backend)
  )
  .matrix_from_image(out, backend) > 0
}

.binary_close_open <- function(mask, close_size = 3L, open_size = 3L) {
  backend <- .directional_support_backend()
  out <- mask * 1

  if (close_size > 1L) {
    brush <- .disc_brush(close_size)
    if (backend == "EBImage") {
      out <- .matrix_from_image(
        EBImage::closing(EBImage::Image(out), brush),
        backend
      )
    } else {
      out <- .matrix_from_image(
        imager::erode(
          imager::dilate(.image_from_matrix(out, backend), .image_from_matrix(brush, backend)),
          .image_from_matrix(brush, backend)
        ),
        backend
      )
    }
  }

  if (open_size > 1L) {
    brush <- .disc_brush(open_size)
    if (backend == "EBImage") {
      out <- .matrix_from_image(
        EBImage::opening(EBImage::Image(out), brush),
        backend
      )
    } else {
      out <- .matrix_from_image(
        imager::dilate(
          imager::erode(.image_from_matrix(out, backend), .image_from_matrix(brush, backend)),
          .image_from_matrix(brush, backend)
        ),
        backend
      )
    }
  }

  out > 0
}

.directional_gaussian_blur <- function(image, sigma) {
  backend <- .directional_support_backend()
  if (backend == "EBImage") {
    return(.matrix_from_image(EBImage::gblur(EBImage::Image(image), sigma = sigma), backend))
  }

  .matrix_from_image(
    imager::isoblur(.image_from_matrix(image, backend), sigma = sigma, neumann = TRUE),
    backend
  )
}

.oriented_gaussian_kernel <- function(theta_deg, sigma_major, sigma_minor) {
  radius <- ceiling(3 * max(sigma_major, sigma_minor))
  grid <- seq(-radius, radius)
  xy <- expand.grid(x = grid, y = grid)
  theta <- theta_deg * pi / 180
  xr <- xy$x * cos(theta) + xy$y * sin(theta)
  yr <- -xy$x * sin(theta) + xy$y * cos(theta)
  ker <- exp(-0.5 * ((xr / sigma_major)^2 + (yr / sigma_minor)^2))
  ker <- matrix(ker, nrow = length(grid), ncol = length(grid), byrow = FALSE)
  ker <- ker / sum(ker)
  ker - mean(ker)
}

.directional_filter_response <- function(image, kernel) {
  backend <- .directional_support_backend()
  if (backend == "EBImage") {
    return(EBImage::filter2(image, kernel, boundary = "replicate"))
  }

  .matrix_from_image(
    imager::correlate(
      .image_from_matrix(image, backend),
      .image_from_matrix(kernel, backend),
      dirichlet = FALSE,
      normalise = FALSE
    ),
    backend
  )
}

.contourlet_like_response <- function(image, angles_deg, scales) {
  nr <- nrow(image)
  nc <- ncol(image)
  best <- matrix(0, nrow = nr, ncol = nc)

  for (scale in scales) {
    sigma_major <- scale[1]
    sigma_minor <- scale[2]
    for (theta in angles_deg) {
      ker <- .oriented_gaussian_kernel(theta, sigma_major, sigma_minor)
      resp <- .directional_filter_response(image, ker)
      resp[!is.finite(resp)] <- 0
      best <- pmax(best, resp, 0, na.rm = TRUE)
    }
  }

  best
}

.robust_threshold <- function(x, k) {
  vals <- x[is.finite(x)]
  if (!length(vals)) {
    return(0)
  }
  med <- stats::median(vals, na.rm = TRUE)
  madv <- stats::mad(vals, center = med, constant = 1, na.rm = TRUE)
  med + k * madv
}

.positive_quantile <- function(x, prob, fallback_k = 2) {
  vals <- x[is.finite(x) & x > 0]
  if (length(vals) >= 10) {
    return(as.numeric(stats::quantile(vals, probs = prob, na.rm = TRUE, names = FALSE)))
  }
  .robust_threshold(x, fallback_k)
}

.grow_mask_constrained <- function(seed, allowed, max_iter = 8L) {
  grown <- seed
  max_iter <- as.integer(max_iter)
  if (!is.finite(max_iter) || max_iter < 1L) {
    return(grown)
  }

  for (iter in seq_len(max_iter)) {
    next_mask <- .binary_dilate(grown, size = 3L) & allowed
    next_mask <- next_mask | seed
    if (identical(next_mask, grown)) {
      break
    }
    grown <- next_mask
  }

  grown
}
