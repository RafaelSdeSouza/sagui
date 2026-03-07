#' Flatten a spectral cube into a pixel-by-band matrix
#'
#' @param x 3-D array or FITS-like list with `imDat`.
#' @return Numeric matrix with one row per spatial pixel and one column per band.
#' @export
cube_to_matrix <- function(x) {
  cube <- if (is.list(x) && !is.null(x$imDat)) x$imDat else x
  stopifnot(is.array(cube), length(dim(cube)) == 3L)

  nx <- dim(cube)[1]
  ny <- dim(cube)[2]
  nb <- dim(cube)[3]

  array(cube, dim = c(nx * ny, nb))
}
