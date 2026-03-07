#' Reconstruct and threshold starlet scales to produce a region mask
#'
#' @param dec `sagui_starlet` object (or a matrix; if matrix, J inferred from `keep_scales`)
#' @param keep_scales integer vector (e.g., 3:6)
#' @param include_coarse logical
#' @param denoise_k NULL or numeric (pre-recon per-scale)
#' @param mode "soft" or "hard" (per-scale)
#' @param threshold "mad" for k·MAD, or "abs" for absolute tau
#' @param k numeric; multiplier for MAD (unused; kept for compat)
#' @param tau numeric; absolute threshold when threshold == "abs"
#' @param positive_only logical; TRUE uses rec>thr; FALSE uses |rec|>thr
#' @return list(soft_rec, hard_rec, mask, sigma, threshold, seeds, candidate)
#' @export
threshold_starlet_regions <- function(
    dec,
    keep_scales    = 2:5,
    include_coarse = FALSE,
    denoise_k      = 3,
    mode           = c("soft","hard"),
    threshold      = c("mad","abs"),
    k              = 3.5,
    tau            = NULL,
    positive_only  = TRUE,
    per_scale_positive = TRUE,
    k_hi = 5.0, k_lo = 3.0,
    area_min = 12,
    keep_negatives = FALSE,
    ring_nbins = 8
){
  sigma_MAD <- function(x){
    xf <- x[is.finite(x)]
    if (!length(xf)) return(NA_real_)
    stats::mad(xf, center = stats::median(xf), constant = 1.4826, na.rm = TRUE)
  }

  mode <- match.arg(mode); threshold <- match.arg(threshold)
  if (is.matrix(dec)) dec <- starlet_mask(dec, J = max(keep_scales))

  dec_local <- dec
  if (per_scale_positive && !is.null(dec_local$w)) {
    dec_local$w <- lapply(dec_local$w, function(w){ w[!is.finite(w)] <- NA_real_; pmax(w,0) })
  }

  rec <- starlet_reconstruct(dec_local, keep_scales, include_coarse,
                             denoise_k, mode, na_policy="preserve")
  rec[!is.finite(rec)] <- NA_real_

  rec_pos <- if (keep_negatives) rec else pmax(rec, 0)
  mag     <- if (positive_only) rec_pos else abs(rec)

  # --- thresholds
  if (threshold == "abs") {
    thr_lo <- if (is.null(tau)) 0 else tau
    thr_hi <- max(thr_lo*1.5, thr_lo + 1e-12)
    sig    <- NA_real_
  } else {
    sig_bg <- .ring_mad(rec_pos, nbins = ring_nbins)
    if (!is.finite(sig_bg) || sig_bg <= 0) sig_bg <- sigma_MAD(rec_pos)
    thr_lo <- k_lo * sig_bg
    thr_hi <- k_hi * sig_bg
    sig    <- sig_bg
  }

  candidate <- (mag > thr_lo) & is.finite(mag)
  seeds     <- (mag > thr_hi) & is.finite(mag)

  keep_mask <- candidate
  if (requireNamespace("imager", quietly = TRUE)) {
    im <- asNamespace("imager")
    to_cimg <- function(m) im$as.cimg(t(m), x=ncol(m), y=nrow(m), cc=1, z=1)
    to_mat  <- function(ci) t(as.matrix(ci))

    # crescer apenas componentes com sementes acima de thr_hi
    lab <- to_mat(im$label(to_cimg(candidate)))
    seed_labels <- unique(lab[seeds & candidate])
    seed_labels <- seed_labels[is.finite(seed_labels) & seed_labels > 0]
    keep_mask <- candidate & (lab %in% seed_labels)

    # remove manchas pequenas
    lab2 <- to_mat(im$label(to_cimg(keep_mask)))
    if (length(unique(lab2)) > 1L) {
      tab <- table(lab2[lab2 > 0])
      small <- as.integer(names(tab[tab < area_min]))
      if (length(small)) keep_mask[ lab2 %in% small ] <- FALSE
    }

    # --- safe closing: prefer imager::closing; else dilate→erode; else no-op
    ci  <- to_cimg(keep_mask)
    nxk <- 3; nyk <- 3
    closing_ok <- "closing" %in% getNamespaceExports("imager")
    if (closing_ok) {
      ci2 <- imager::closing(ci, nx = nxk, ny = nyk)
    } else if (all(c("dilate_square","erode_square") %in% getNamespaceExports("imager"))) {
      s   <- max(nxk, nyk)
      ci2 <- im$erode_square(im$dilate_square(ci, size = s), size = s)
    } else if (requireNamespace("imagerExtra", quietly = TRUE)) {
      ci2 <- imagerExtra::mclosing(ci, size = max(nxk, nyk))
    } else {
      ci2 <- ci
    }
    keep_mask <- to_mat(im$threshold(ci2, thr = 0.5)) > 0.5
  } else {
    keep_mask[is.na(keep_mask)] <- FALSE
  }

  hard_rec <- rec_pos
  hard_rec[!keep_mask] <- 0

  list(
    soft_rec   = rec,
    hard_rec   = hard_rec,
    mask       = keep_mask,
    sigma      = sig,
    threshold  = c(low=thr_lo, high=thr_hi),
    seeds      = seeds,
    candidate  = candidate
  )
}
