#' Mask a 3D cube by a spatial label map
#'
#' Sets pixels outside selected regions to zero or NA across all bands.
#'
#' @param cube 3-D numeric array with dimensions `[H, W, B]`.
#' @param labels 2-D matrix `[H, W]`; values greater than zero are kept.
#' @param mode "zero" or "na" — how to fill masked pixels.
#' @return Masked cube with same dimensions as input.
#' @export
mask_cube <- function(cube, labels, mode = c("zero","na")) {
  mode <- match.arg(mode)
  stopifnot(is.array(cube), length(dim(cube)) == 3L, is.matrix(labels))
  H <- dim(cube)[1]; W <- dim(cube)[2]; B <- dim(cube)[3]
  stopifnot(identical(dim(labels), c(H, W)))
  m <- (labels > 0)
  m3 <- array(m, dim = c(H, W, B))
  out <- cube
  if (mode == "zero") {
    out[!m3] <- 0
  } else {
    storage.mode(out) <- "double"
    out[!m3] <- NA_real_
  }
  out
}
