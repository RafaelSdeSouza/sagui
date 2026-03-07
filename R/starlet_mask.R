#' Starlet (à trous) decomposition
#' @param img 2D numeric matrix
#' @param J number of scales (>=1)
#' @return object of class \code{sagui_starlet}: list(w = list w1..wJ, cJ = coarse)
#' @export
starlet_mask  <- function(img, J = 5) {
  stopifnot(is.matrix(img))
#  mask <- is.finite(img) * 1.0
  mask_valid <- is.finite(img)

  img0 <- img
  img0[!mask_valid] <- 0
#  img0[!is.finite(img0)] <- 0

  k0 <- b3_kernel()
  c_j <- img0
  wlist <- vector("list", J)

  for (j in seq_len(J)) {
    step <- 2^(j - 1)
    k_j  <- upsample_kernel(k0, step)
    smooth <- .smooth_sep_masked(c_j, k_j, mask_valid)

     w     <- c_j - smooth
    wlist[[j]] <- w
    c_j <- smooth
  }

  # keep original invalid pixels invalid in outputs (optional)
  for (j in seq_len(J)) wlist[[j]][!mask_valid] <- NA_real_
  c_j[!mask_valid] <- NA_real_
  structure(
    list(w = wlist, cJ = c_j, mask = mask_valid),
    class = c("sagui_starlet", "list")
  )

}
