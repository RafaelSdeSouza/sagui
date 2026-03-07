#' @export
asinh_stretch <- function(x, qlo = 0.001, qhi = 0.999, scale = NULL,
                          nonneg = TRUE, na.rm = TRUE) {
  z <- x

  if (na.rm)
    z[!is.finite(z)] <- NA_real_

  if (nonneg)
    z[z < 0] <- 0

  lo <- as.numeric(stats::quantile(z, probs = qlo, na.rm = TRUE, names = FALSE))
  hi <- as.numeric(stats::quantile(z, probs = qhi, na.rm = TRUE, names = FALSE))
  if (!is.finite(lo)) lo <- 0
  if (!is.finite(hi) || hi <= lo) hi <- max(z, na.rm = TRUE)

  z <- pmax(z - lo, 0)

  if (is.null(scale)) {
    scale <- (hi - lo) / 20
    if (!is.finite(scale) || scale == 0) scale <- 1
  }

  z2 <- asinh(z / scale)
  z2 <- z2 / max(z2, na.rm = TRUE)

  z2[!is.finite(z2)] <- NA_real_
  z2
}
