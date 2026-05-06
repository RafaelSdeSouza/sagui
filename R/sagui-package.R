#' sagui
#'
#' Photometric analysis tools for spectral cubes.
#'
#' @keywords internal
#' @useDynLib sagui, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats setNames var
"_PACKAGE"

utils::globalVariables(c(
  "band", "Cluster", "Col", "flux", "flux_err", "geometry", "n_eff",
  "n_pix", "ok", "region", "Row", "Value"
))
