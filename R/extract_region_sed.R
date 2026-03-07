#' Extract integrated SEDs for segmented regions
#'
#' @param cube A 3-D array with dimensions `[nx, ny, nb]` or a FITS-like list
#'   with `imDat`.
#' @param labels A matrix `[nx, ny]` with positive region identifiers.
#' @param bkg Deprecated compatibility argument. Ignored.
#' @param var_cube Variance cube with the same shape as `cube`, or a FITS-like
#'   list.
#' @param sigma_band Per-band standard deviation vector used when `var_cube` is
#'   `NULL`.
#' @param band_values Optional band labels or wavelengths. If omitted and `cube`
#'   is FITS-like, `FITSio::axVec()` is used.
#' @param digits_lambda_colnames Significant digits used when numeric band values
#'   are converted to column names.
#' @param return_painted_cube Logical; if `TRUE`, also return a cube painted
#'   with region fluxes.
#' @param error_fallback Fallback error model used when no variance cube is
#'   provided.
#' @return A list with `flux_long`, `flux_wide`, `painted_cube`, `bands`, and
#'   `band_names`.
#' @export
extract_region_sed <- function(cube,
                               labels,
                               bkg = NULL,
                               var_cube = NULL,
                               sigma_band = NULL,
                               band_values = NULL,
                               digits_lambda_colnames = 6,
                               return_painted_cube = FALSE,
                               error_fallback = c("none","flux_over_sqrt_n","poisson","mad_sky")) {
  RegionPhotometry(
    cube = cube,
    labels = labels,
    bkg = bkg,
    var_cube = var_cube,
    sigma_band = sigma_band,
    band_values = band_values,
    digits_lambda_colnames = digits_lambda_colnames,
    return_painted_cube = return_painted_cube,
    error_fallback = error_fallback
  )
}

#' Summarize integrated SEDs for segmented regions
#'
#' @inheritParams extract_region_sed
#' @return A list with `flux_long`, `flux_wide`, `painted_cube`, `bands`, and
#'   `band_names`.
#' @export
summarize_region_sed <- function(cube,
                                 labels,
                                 bkg = NULL,
                                 var_cube = NULL,
                                 sigma_band = NULL,
                                 band_values = NULL,
                                 digits_lambda_colnames = 6,
                                 return_painted_cube = FALSE,
                                 error_fallback = c("none","flux_over_sqrt_n","poisson","mad_sky")) {
  extract_region_sed(
    cube = cube,
    labels = labels,
    bkg = bkg,
    var_cube = var_cube,
    sigma_band = sigma_band,
    band_values = band_values,
    digits_lambda_colnames = digits_lambda_colnames,
    return_painted_cube = return_painted_cube,
    error_fallback = error_fallback
  )
}
