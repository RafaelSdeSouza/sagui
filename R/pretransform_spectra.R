#' Apply a spectral pretransform before segmentation
#'
#' @description
#' Applies a column-wise pretransform to a spectra matrix before row-wise
#' scaling and clustering. This is useful for benchmarking clustering behavior
#' under simple nonlinear transforms such as `asinh()`, `log1p()`, or
#' rank-based copula mappings.
#'
#' @param x Numeric matrix with one spectrum per row and one band per column.
#' @param method Either a string naming a built-in transform or a function that
#'   takes a numeric matrix and returns a numeric matrix with the same
#'   dimensions. Built-in options are `"none"`, `"asinh"`, `"log1p"`,
#'   `"signed_log1p"`, `"copula_uniform"`, and `"copula_gaussian"`.
#'
#' @return A numeric matrix with the same dimensions as `x`.
#' @export
pretransform_spectra <- function(
    x,
    method = c(
      "none",
      "asinh",
      "log1p",
      "signed_log1p",
      "copula_uniform",
      "copula_gaussian"
    )
) {
  if (!is.matrix(x)) {
    stop("`x` must be a numeric matrix with spectra in rows.")
  }
  if (!is.numeric(x)) {
    stop("`x` must be numeric.")
  }

  if (is.function(method)) {
    out <- method(x)
    if (!is.matrix(out) || !is.numeric(out) || !identical(dim(out), dim(x))) {
      stop("A custom pretransform must return a numeric matrix with the same dimensions as `x`.")
    }
    return(out)
  }

  method <- match.arg(method)

  out <- switch(
    method,
    none = x,
    asinh = asinh(x),
    log1p = {
      xf <- x[is.finite(x)]
      if (any(xf < 0)) {
        stop("`log1p` pretransform requires all finite values to be non-negative.")
      }
      log1p(x)
    },
    signed_log1p = sign(x) * log1p(abs(x)),
    copula_uniform = .empirical_copula_transform(x, gaussian = FALSE),
    copula_gaussian = .empirical_copula_transform(x, gaussian = TRUE)
  )

  storage.mode(out) <- "double"
  out
}

#' Apply a spectral pretransform to a full cube
#'
#' @param x 3-D array or FITS-like list with `imDat`.
#' @param method Either a built-in transform name or a custom function accepted
#'   by `pretransform_spectra()`.
#'
#' @return A 3-D numeric array with the same dimensions as the input cube.
#' @export
# Apply the same pretransform to a full cube by flattening to spectra,
# transforming band-wise, then reshaping back to the original cube geometry.
pretransform_cube <- function(x, method = "none") {
  cube <- if (is.list(x) && !is.null(x$imDat)) x$imDat else x
  if (!is.array(cube) || length(dim(cube)) != 3L) {
    stop("`x` must be a 3-D array or a FITS-like list with `imDat`.")
  }

  transformed <- pretransform_spectra(cube_to_matrix(cube), method = method)
  array(transformed, dim = dim(cube), dimnames = dimnames(cube))
}

.empirical_copula_transform <- function(x, gaussian = FALSE) {
  out <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))

  for (j in seq_len(ncol(x))) {
    column <- x[, j]
    ok <- is.finite(column)
    n_ok <- sum(ok)
    if (!n_ok) {
      next
    }

    ranks <- rank(column[ok], ties.method = "average")
    u <- (ranks - 0.5) / n_ok
    if (isTRUE(gaussian)) {
      u <- stats::qnorm(u)
    }

    out[ok, j] <- u
  }

  out
}
