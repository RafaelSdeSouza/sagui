#' Reconstruct from a subset of starlet scales
#' @param dec \code{sagui_starlet} object
#' @param keep_scales integer vector of scales to include (e.g., 3:6)
#' @param include_coarse logical; add coarse plane \code{cJ}
#' @param denoise_k NULL or numeric (k·sigma MAD per-plane)
#' @param mode "soft" or "hard" thresholding
#' @param na_policy "preserve" or "zero" for NA handling during sum
#' @export
starlet_reconstruct <- function(dec,
                                keep_scales = c(2,3),
                                include_coarse = FALSE,
                                denoise_k = NULL,          # scalar or vector length J
                                mode = c("soft","hard"),
                                na_policy = c("preserve","zero")) {
  mode <- match.arg(mode)
  na_policy <- match.arg(na_policy)

  wlist <- dec$w
  J <- length(wlist)

  stopifnot(is.numeric(keep_scales), all(keep_scales >= 1), all(keep_scales <= J))

  # 1) Optional denoising per scale (do NOT move NAs)
  if (!is.null(denoise_k)) {
    kvec <- if (length(denoise_k) == 1) rep(denoise_k, J) else denoise_k
    stopifnot(length(kvec) == J)

    for (j in seq_len(J)) {
      w <- wlist[[j]]
      idx <- is.finite(w)
      if (any(idx)) {
        sig <- .estimate_sigma(w[idx])
        thr <- kvec[j] * sig
        if (is.finite(thr) && thr > 0) {
          if (mode == "soft") {
            w[idx] <- .soft_thresh(w[idx], thr)
          } else {
            w[idx] <- w[idx] * (abs(w[idx]) >= thr)
          }
        }
      }
      wlist[[j]] <- w
    }
  }

  # 2) Select planes
  selected <- wlist[keep_scales]

  # 3) Sum with NA handling
  if (na_policy == "zero") {
    # present counts how many components contribute finite values per pixel
    present <- Reduce(`+`, lapply(selected, function(w) as.numeric(is.finite(w))))
    selected0 <- lapply(selected, function(w) { w2 <- w; w2[!is.finite(w2)] <- 0; w2 })
    rec <- Reduce(`+`, selected0)

    if (isTRUE(include_coarse)) {
      cJ <- dec$cJ
      present <- present + as.numeric(is.finite(cJ))
      cJ0 <- cJ; cJ0[!is.finite(cJ0)] <- 0
      rec <- rec + cJ0
    }

    rec[present == 0] <- NA_real_
  } else {
    rec <- Reduce(`+`, selected)
    if (isTRUE(include_coarse)) rec <- rec + dec$cJ
  }

  # 4) Enforce original invalid pixels (assume logical validity mask if present)
  if (!is.null(dec$mask)) {
    if (is.logical(dec$mask)) rec[!dec$mask] <- NA_real_ else rec[dec$mask == 0] <- NA_real_
  }

  rec
}
