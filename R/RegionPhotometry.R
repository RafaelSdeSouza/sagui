#' Region-level photometry on hyperspectral cubes (no background subtraction)
#'
#' @description
#' Aggregates flux by region and band for a hyperspectral or IFU cube.
#' Variance-aware uncertainties are used when available, with simple fallback
#' error models otherwise. A numeric `lambda` column is also added when band
#' values can be interpreted as wavelengths.
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
#'
#' @return A list with `flux_long`, `flux_wide`, `painted_cube`, `bands`, and
#'   `band_names`.
#'
#' @export
#' @importFrom dplyr group_by summarise mutate select arrange left_join rename_with first count
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble as_tibble tibble
#' @importFrom stats mad median
#' @importFrom FITSio axVec
RegionPhotometry <- function(
    cube, labels,
    bkg = NULL,                 # mantido, porém ignorado
    var_cube = NULL,
    sigma_band = NULL,
    band_values = NULL,
    digits_lambda_colnames = 6,
    return_painted_cube = FALSE,
    error_fallback = c("none","flux_over_sqrt_n","poisson","mad_sky")
) {
  error_fallback <- match.arg(error_fallback)

  # --- helpers ---
  get_imdat <- function(x) if (is.list(x) && !is.null(x$imDat)) x$imDat else x

  # [nx,ny,nb] -> [n_pix, nb], colunas = bandas
  to_mat <- function(A3) {
    nb <- dim(A3)[3]
    t(matrix(aperm(A3, c(3, 1, 2)), nrow = nb))
  }

  # --- inputs & bands ---
  M <- get_imdat(cube)
  stopifnot(is.array(M), length(dim(M)) == 3L)
  nx <- dim(M)[1]; ny <- dim(M)[2]; nb <- dim(M)[3]

  # labels -> [nx,ny]
  if (is.array(labels) && length(dim(labels)) == 3L && dim(labels)[3] == 1L) {
    labels <- labels[,,1]
  }
  stopifnot(is.matrix(labels), all(dim(labels)[1:2] == c(nx, ny)))

  # bands (valores + rótulos)
  if (!is.null(band_values)) {
    bands <- band_values
  } else if (is.list(cube) && !is.null(cube$axDat)) {
    bands <- FITSio::axVec(3, cube$axDat)
  } else {
    bands <- seq_len(nb)
  }
  stopifnot(length(bands) == nb)

  # rótulos seguros (sem padding) + mapeamento numérico para lambda
  if (is.numeric(bands)) {
    band_names <- sprintf("%.*g", digits_lambda_colnames, bands)
  } else {
    band_names <- as.character(bands)
  }
  band_names <- trimws(band_names)
  band_names <- make.unique(band_names, sep = "_")

  band_map <- tibble::tibble(
    band   = band_names,
    lambda = suppressWarnings(as.numeric(bands))  # numérico se possível
  )

  # --- background: ignorado (compatibilidade)
  if (!is.null(bkg)) {
    warning("`bkg` is deprecated and ignored in `RegionPhotometry()`.")
  }

  # --- variância ---
  use_var <- FALSE
  if (!is.null(var_cube)) {
    Vc <- get_imdat(var_cube); stopifnot(is.array(Vc), all(dim(Vc) == dim(M)))
    use_var <- TRUE
  } else if (!is.null(sigma_band)) {
    stopifnot(length(sigma_band) == nb)
    Vc <- array(0, dim = dim(M))
    for (j in seq_len(nb)) Vc[,,j] <- sigma_band[j]^2
    use_var <- TRUE
  } else {
    Vc <- NULL
  }

  # --- flatten ---
  X   <- to_mat(M)                         # [n_pix_total, nb]
  cls <- as.vector(labels)                 # [n_pix_total]
  valid <- is.finite(cls) & (cls > 0)

  if (!any(valid)) {
    warning("No valid pixels found (`labels > 0` and finite). Returning empty tables.")
    empty <- tibble::tibble(region = integer(), band = character(), flux = numeric(),
                            flux_err = numeric(), n_eff = integer(), n_pix = integer(),
                            lambda = numeric())
    return(list(
      flux_long    = empty,
      flux_wide    = tibble::tibble(),
      painted_cube = if (isTRUE(return_painted_cube)) array(NA_real_, dim = dim(M)) else NULL,
      bands        = bands,
      band_names   = band_names
    ))
  }

  Xv   <- X[valid, , drop = FALSE]
  clsv <- cls[valid]
  pix  <- seq_len(nrow(Xv))

  flux_df <- tibble::as_tibble(Xv, .name_repair = "minimal"); names(flux_df) <- band_names
  flux_long0 <- flux_df |>
    dplyr::mutate(pix = pix) |>
    tidyr::pivot_longer(-pix, names_to = "band", values_to = "flux")

  if (use_var) {
    V   <- to_mat(Vc)[valid, , drop = FALSE]
    var_df <- tibble::as_tibble(V, .name_repair = "minimal"); names(var_df) <- band_names
    var_long <- var_df |>
      dplyr::mutate(pix = pix) |>
      tidyr::pivot_longer(-pix, names_to = "band", values_to = "var")
    df <- dplyr::left_join(flux_long0, var_long, by = c("pix","band"))
  } else {
    df <- flux_long0 |>
      dplyr::mutate(var = NA_real_)
  }

  # região + flag ok
  df2 <- df |>
    dplyr::mutate(
      region = clsv[pix],
      ok = if (use_var) (is.finite(flux) & is.finite(var)) else is.finite(flux)
    )

  # --- agregação por região × banda ---
  df_sum <- df2 |>
    dplyr::group_by(region, band) |>
    dplyr::summarise(
      n_eff   = sum(ok),
      flux    = sum(ifelse(ok, flux, 0), na.rm = TRUE),
      var_sum = if (use_var) sum(ifelse(ok, var, 0), na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      flux_err = dplyr::case_when(
        use_var ~ sqrt(pmax(var_sum, 0)),
        !use_var & identical(error_fallback, "flux_over_sqrt_n") ~ abs(flux) / sqrt(pmax(n_eff, 1)),
        !use_var & identical(error_fallback, "poisson")          ~ sqrt(pmax(flux, 0)),
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::select(region, band, flux, flux_err, n_eff)

  # fallback "mad_sky" (estima sigma por banda em pixels de céu)
  if (!use_var && identical(error_fallback, "mad_sky")) {
    sky <- !is.finite(cls) | (cls <= 0)
    if (!any(sky)) {
      warning("Sem pixels de céu (labels <= 0). Caindo para 'flux_over_sqrt_n'.")
      df_sum <- df_sum |>
        dplyr::mutate(flux_err = ifelse(is.finite(flux_err), flux_err, abs(flux)/sqrt(pmax(n_eff,1))))
    } else {
      Xsky <- X[sky, , drop = FALSE]
      sky_sigma <- apply(
        Xsky, 2L,
        function(v) 1.4826 * stats::mad(v, center = stats::median(v, na.rm = TRUE), na.rm = TRUE)
      )
      names(sky_sigma) <- band_names
      df_sum <- df_sum |>
        dplyr::mutate(flux_err = sky_sigma[band] * sqrt(pmax(n_eff, 1)))
    }
  }

  # anexar lambda (numérica) para plot
  df_sum <- dplyr::left_join(df_sum, band_map, by = "band")

  # contagem de pixels por região
  n_pix_region <- tibble::tibble(region = clsv) |>
    dplyr::count(region, name = "n_pix")

  # long final
  flux_long <- df_sum |>
    dplyr::left_join(n_pix_region, by = "region") |>
    dplyr::arrange(region, band)

  # --- WIDE tables ---
  flux_wide <- flux_long |>
    dplyr::select(region, n_pix, band, flux) |>
    tidyr::pivot_wider(names_from = band, values_from = flux, values_fn = dplyr::first) |>
    dplyr::arrange(region)

  err_wide <- flux_long |>
    dplyr::select(region, band, flux_err) |>
    tidyr::pivot_wider(names_from = band, values_from = flux_err, values_fn = dplyr::first) |>
    dplyr::rename_with(~ paste0(.x, "_err"), -region)

  neff_wide <- flux_long |>
    dplyr::select(region, band, n_eff) |>
    tidyr::pivot_wider(names_from = band, values_from = n_eff, values_fn = dplyr::first) |>
    dplyr::rename_with(~ paste0(.x, "_n_eff"), -region)

  flux_wide <- flux_wide |>
    dplyr::left_join(neff_wide, by = "region") |>
    dplyr::left_join(err_wide, by = "region")

  # --- painted cube (opcional)
  painted_cube <- NULL
  if (isTRUE(return_painted_cube)) {
    painted_cube <- array(NA_real_, dim = dim(M))
    reg_vec <- as.vector(labels)
    look <- flux_long |>
      dplyr::select(region, band, flux)

    for (j in seq_len(nb)) {
      bname <- band_names[j]
      sel <- look$band == bname
      regs_b <- look$region[sel]
      flux_b <- look$flux[sel]
      idx <- match(reg_vec, regs_b)
      layer <- flux_b[idx]
      layer[!(is.finite(reg_vec) & reg_vec > 0)] <- NA_real_
      painted_cube[,,j] <- matrix(layer, nrow = nx, ncol = ny)
    }
    dimnames(painted_cube) <- list(NULL, NULL, band_names)
  }

  list(
    flux_long    = flux_long,   # tem 'lambda' numérica
    flux_wide    = flux_wide,
    painted_cube = painted_cube,
    bands        = bands,
    band_names   = band_names
  )
}
