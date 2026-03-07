#' Insert zeros between kernel taps
#' @param k base kernel (numeric)
#' @param step integer upsampling factor (>=1)
#' @return upsampled kernel
#' @export
upsample_kernel <- function(k, step) {
  if (step <= 1) return(k)
  out <- rep(0, length(k) + (length(k) - 1) * (step - 1))
  out[seq(1, length(out), by = step)] <- k
  out
}
