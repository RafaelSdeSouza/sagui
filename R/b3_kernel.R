#' B3-spline low-pass kernel (1D)
#' @return numeric vector of length 5
#' @export
b3_kernel <- function() c(1, 4, 6, 4, 1) / 16
