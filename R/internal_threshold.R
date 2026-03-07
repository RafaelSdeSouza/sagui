#' @keywords ring
#' @noRd
.ring_mad <- function(m, cx=NULL, cy=NULL, nbins=8, r_in_frac=0.15, r_out_frac=0.95){
  stopifnot(is.matrix(m))
  nY <- nrow(m); nX <- ncol(m)
  if (is.null(cx)) cx <- (nX+1)/2; if (is.null(cy)) cy <- (nY+1)/2
  xx <- matrix(rep(1:nX, each=nY), nY, nX)
  yy <- matrix(rep(1:nY, nX), nY, nX)
  r  <- sqrt((xx-cx)^2 + (yy-cy)^2)
  r  <- r / max(r, na.rm=TRUE)

  brks <- seq(r_in_frac, r_out_frac, length.out=nbins+1)
  madv <- numeric(nbins); wgt <- numeric(nbins)
  for(i in seq_len(nbins)){
    sel <- is.finite(m) & r >= brks[i] & r < brks[i+1]
    if (any(sel)) {
      madv[i] <- stats::mad(m[sel], center=stats::median(m[sel]), constant=1.4826, na.rm=TRUE)
      wgt[i]  <- sum(sel)
    } else { madv[i] <- NA_real_; wgt[i] <- 0 }
  }
  # média ponderada dos MADs de anel
  mad <- sum(madv*wgt, na.rm=TRUE) / max(1,sum(wgt[wgt>0]))
  if(!is.finite(mad) || mad<=0) mad <- stats::mad(as.numeric(m[is.finite(m)]),
                                                  center=stats::median(m[is.finite(m)]),
                                                  constant=1.4826, na.rm=TRUE)
  mad
}
#' @keywords internal
#' @noRd
.soft_thresh <- function(x, t) sign(x) * pmax(abs(x) - t, 0)

#' @keywords internal
#' @noRd
.hard_thresh <- function(x, t) ifelse(abs(x) >= t, x, 0)

#' Robust per-plane noise estimate (MAD / 0.6745)
#' @keywords internal
#' @noRd
.estimate_sigma <- function(w) stats::mad(w, center = 0, constant = 1, na.rm = TRUE)


#' Robust sigma via MAD
#' @param x numeric vector/matrix
#' @return scalar robust sigma
#' @keywords internal
#' @noRd
.sigma_MAD <- function(x) 1.4826 * stats::median(abs(x - stats::median(x, na.rm=TRUE)), na.rm=TRUE)

#' Binary mask by k·sigma threshold on a matrix
#' @param img numeric matrix (band-pass reconstruction)
#' @param k numeric threshold multiplier (default 3.5)
#' @param positive_only logical; keep only positive structures
#' @return logical matrix mask
#' @keywords internal
#' @noRd
.make_sigma_mask <- function(img, k = 3.5, positive_only = TRUE) {
  sig <- .sigma_MAD(img)
  if (!is.finite(sig) || sig <= 0) return(matrix(FALSE, nrow(img), ncol(img)))
  if (positive_only) {
    (img >  k * sig) & is.finite(img)
  } else {
    (abs(img) > k * sig) & is.finite(img)
  }
}
