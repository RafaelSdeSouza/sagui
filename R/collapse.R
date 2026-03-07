#' Robust per-wavelength background and scatter
#'
#' Computes the spatial median (background) and MAD-based scatter per wavelength
#' slice of a 3-D cube.
#'
#' @param cube 3-D numeric array with dimensions `[nx, ny, nlambda]`.
#' @return List with numeric vectors `bkg` and `mad`, both of length `nlambda`.
#' @export
estimate_bkg_mad_per_lambda <- function(cube) {
  stopifnot(length(dim(cube)) == 3L)
  nx <- dim(cube)[1]; ny <- dim(cube)[2]; nlam <- dim(cube)[3]
  M  <- matrix(cube, nrow = nx * ny, ncol = nlam)
  
  bkg <- matrixStats::colMedians(M, na.rm = TRUE)
  
  # MAD-based scatter per λ (scaled once to be sigma-consistent)
  mad <- matrixStats::colMads(M, na.rm = TRUE)   # <- remove 1.4826 *
  
  # guard against zero/NA MAD
  repl <- stats::median(mad[is.finite(mad) & mad > 0], na.rm = TRUE)
  if (!is.finite(repl) || repl <= 0) repl <- 1
  mad[!is.finite(mad) | mad <= 0] <- repl
  
  list(bkg = bkg, mad = mad)
}

#' Collapse a spectral cube into a white-light image (robust, weighted)
#'
#' Steps:
#' 1) subtract spatial median per-wavelength,
#' 2) clip very negative outliers at -k * MAD,
#' 3) (optional) inverse-variance weights 1/MAD^2,
#' 4) weighted sum and reshape to a 2-D image,
#' 5) force non-negative for visualization safety.
#'
#' @param cube 3-D numeric array with dimensions `[nx, ny, nlambda]`.
#' @param kclip numeric; negative clip in units of MAD (default 2)
#' @param use_weights logical; if TRUE, use 1/MAD^2 weights (default TRUE)
#' @return 2-D matrix with dimensions `[nx, ny]`.
#' @export
collapse_white_light <- function(cube, kclip = 2, use_weights = TRUE) {
  stopifnot(length(dim(cube)) == 3L)
  est <- estimate_bkg_mad_per_lambda(cube)
  bkg <- est$bkg; mad <- est$mad
  
  nx <- dim(cube)[1]; ny <- dim(cube)[2]; nlam <- dim(cube)[3]
  
  # subtract background per λ
  for (k in seq_len(nlam)) cube[,,k] <- cube[,,k] - bkg[k]
  
  # clip strong negatives; non-finite -> 0
  for (k in seq_len(nlam)) {
    thr <- -kclip * mad[k]
    sli <- cube[,,k]
    sli[!is.finite(sli)] <- 0
    sli[sli < thr] <- thr            # <- was 0 (minimal but important change)
    cube[,,k] <- sli
  }
  
  # weights
  if (isTRUE(use_weights)) {
    w <- 1 / (mad^2)
    w <- w / stats::median(w[is.finite(w) & w > 0], na.rm = TRUE)
  } else {
    w <- rep(1, length(mad))
  }
  
  M <- matrix(cube, nrow = nx * ny, ncol = nlam)
  wl <- as.numeric(M %*% w)
  wl <- matrix(wl, nrow = nx, ncol = ny)
  
  # optional; you can keep it for “visualization safety”
  wl[wl < 0] <- 0
  wl
}
