#' Center a spectrum by its median
#'
#' @param x Numeric vector.
#' @return Numeric vector centered on the median of `x`.
#' @export
median_scale <- function(x) {
  if (!is.numeric(x)) {
    stop("`x` must be numeric.")
  }

  x - stats::median(x, na.rm = TRUE)
}
