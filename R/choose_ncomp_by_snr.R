.finite_summary <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(list(min = NA_real_, median = NA_real_, q16 = NA_real_, q84 = NA_real_))
  }

  list(
    min = min(x),
    median = stats::median(x),
    q16 = unname(stats::quantile(x, 0.16, na.rm = TRUE)),
    q84 = unname(stats::quantile(x, 0.84, na.rm = TRUE))
  )
}

.compute_region_snr_from_sed <- function(flux_long,
                                         wavelength_range = NULL,
                                         snr_stat = c("integrated", "median_per_wavelength"),
                                         variance_inflation = 1) {
  snr_stat <- match.arg(snr_stat)

  required <- c("region", "flux", "flux_err")
  missing <- setdiff(required, names(flux_long))
  if (length(missing)) {
    stop("`flux_long` is missing required columns: ", paste(missing, collapse = ", "))
  }

  df <- flux_long
  if (!is.null(wavelength_range)) {
    if (!("lambda" %in% names(df))) {
      stop("`wavelength_range` requires a `lambda` column in `flux_long`.")
    }
    if (length(wavelength_range) != 2L) {
      stop("`wavelength_range` must have length 2.")
    }
    lo <- min(wavelength_range)
    hi <- max(wavelength_range)
    df <- df[is.finite(df$lambda) & df$lambda >= lo & df$lambda <= hi, , drop = FALSE]
    if (!nrow(df)) {
      stop("No wavelengths fall inside `wavelength_range`.")
    }
  }

  split_df <- split(df, df$region)
  vapply(split_df, function(one) {
    good <- is.finite(one$flux) & is.finite(one$flux_err) & one$flux_err > 0
    if (!any(good)) {
      return(NA_real_)
    }
    one <- one[good, , drop = FALSE]

    if (snr_stat == "integrated") {
      total_flux <- sum(one$flux, na.rm = TRUE)
      total_var <- variance_inflation * sum(one$flux_err^2, na.rm = TRUE)
      if (!is.finite(total_var) || total_var <= 0) {
        return(NA_real_)
      }
      abs(total_flux) / sqrt(total_var)
    } else {
      stats::median(abs(one$flux) / sqrt(variance_inflation * one$flux_err^2), na.rm = TRUE)
    }
  }, numeric(1))
}

.resolve_sagui_snr_target <- function(min_snr,
                                      min_SNR,
                                      median_snr,
                                      target_snr,
                                      screen_stat) {
  supplied <- c(
    min_snr = !is.null(min_snr),
    min_SNR = !is.null(min_SNR),
    median_snr = !is.null(median_snr),
    target_snr = !is.null(target_snr)
  )

  if (!any(supplied)) {
    stop("Supply `min_snr` for a minimum regional S/N target.")
  }
  if (sum(supplied) > 1L) {
    stop("Use only one of `min_snr`, `min_SNR`, `median_snr`, or `target_snr`.")
  }

  if (!is.null(min_snr)) {
    threshold <- min_snr
    criterion <- "min"
    label <- "min_snr"
  } else if (!is.null(min_SNR)) {
    threshold <- min_SNR
    criterion <- "min"
    label <- "min_SNR"
  } else if (!is.null(median_snr)) {
    threshold <- median_snr
    criterion <- "median"
    label <- "median_snr"
  } else {
    threshold <- target_snr
    criterion <- match.arg(screen_stat)
    label <- "target_snr"
  }

  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) || threshold <= 0) {
    stop("The requested S/N threshold must be a positive finite number.")
  }

  list(
    threshold = threshold,
    criterion = criterion,
    label = label
  )
}

#' Choose the Number of SAGUI Regions from a Target S/N
#'
#' This helper runs SAGUI over a grid of candidate region counts and keeps the
#' finest segmentation satisfying a regional S/N cut. Use `min_snr` for the
#' common case where every region must exceed a minimum S/N. Unlike the quick
#' `cluster_snr` diagnostic returned by [segment_regions()], this function can
#' evaluate the S/N from the flux-conserving regional SED table produced by
#' [extract_region_sed()].
#'
#' @param input A 3-D array or FITS-like list with `imDat`.
#' @param min_snr Minimum regional S/N. This is the recommended user-facing
#'   argument.
#' @param min_SNR Alias for `min_snr`, accepted for readability.
#' @param median_snr Optional median regional S/N target. This is less strict
#'   than `min_snr`.
#' @param target_snr Deprecated compatibility alias. Use `min_snr` or
#'   `median_snr`.
#' @param var_cube Optional variance cube with the same dimensions as `input`.
#' @param k_values Candidate numbers of regions to test.
#' @param wavelength_range Optional numeric interval used to select wavelengths
#'   for the S/N calculation.
#' @param band_values Optional band labels or wavelengths passed to
#'   [extract_region_sed()].
#' @param backend Clustering backend: `"exact"` calls [segment_regions()] and
#'   `"sparse"` calls [segment_regions_large()].
#' @param snr_stat Either integrated S/N across the selected wavelengths or the
#'   median per-wavelength S/N.
#' @param screen_stat Deprecated compatibility option used only with
#'   `target_snr`. Prefer `min_snr` or `median_snr`.
#' @param variance_inflation Multiplicative factor applied to propagated
#'   variances.
#' @param error_fallback Fallback error model used by [extract_region_sed()]
#'   when no `var_cube` is provided.
#' @param verbose Logical; print progress messages.
#' @param ... Additional arguments passed to [segment_regions()] or
#'   [segment_regions_large()].
#'
#' @return A SAGUI segmentation result for the selected number of regions. The
#'   result also contains `region_snr`, `snr_grid`, `snr_target`, and
#'   `snr_selection` entries.
#' @export
choose_ncomp_by_snr <- function(input,
                                min_snr = NULL,
                                min_SNR = NULL,
                                median_snr = NULL,
                                target_snr = NULL,
                                var_cube = NULL,
                                k_values = c(5, 10, 15, 20, 30, 40, 60),
                                wavelength_range = NULL,
                                band_values = NULL,
                                backend = c("exact", "sparse"),
                                snr_stat = c("integrated", "median_per_wavelength"),
                                screen_stat = c("min", "median"),
                                variance_inflation = 1,
                                error_fallback = c("mad_sky", "none", "flux_over_sqrt_n", "poisson"),
                                verbose = TRUE,
                                ...) {
  backend <- match.arg(backend)
  snr_stat <- match.arg(snr_stat)
  error_fallback <- match.arg(error_fallback)
  snr_target <- .resolve_sagui_snr_target(
    min_snr = min_snr,
    min_SNR = min_SNR,
    median_snr = median_snr,
    target_snr = target_snr,
    screen_stat = screen_stat
  )

  k_values <- sort(unique(as.integer(k_values)), decreasing = TRUE)
  k_values <- k_values[is.finite(k_values) & k_values >= 1L]
  if (!length(k_values)) {
    stop("`k_values` must contain at least one positive integer.")
  }

  segment_fun <- switch(
    backend,
    exact = segment_regions,
    sparse = segment_regions_large
  )

  fits <- vector("list", length(k_values))
  names(fits) <- as.character(k_values)
  region_snr_list <- vector("list", length(k_values))
  names(region_snr_list) <- as.character(k_values)

  rows <- lapply(seq_along(k_values), function(i) {
    k <- k_values[i]
    if (isTRUE(verbose)) {
      message("Testing Ncomp = ", k)
    }

    fit <- tryCatch(
      segment_fun(input = input, Ncomp = k, ...),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(data.frame(
        Ncomp = k,
        actual_Ncomp = NA_integer_,
        min_region_snr = NA_real_,
        median_region_snr = NA_real_,
        q16_region_snr = NA_real_,
        q84_region_snr = NA_real_,
        screen_snr = NA_real_,
        all_regions_finite = FALSE,
        pass = FALSE,
        error = conditionMessage(fit)
      ))
    }

    sed <- extract_region_sed(
      cube = input,
      labels = fit$cluster_map,
      var_cube = var_cube,
      band_values = band_values,
      error_fallback = if (is.null(var_cube)) error_fallback else "none"
    )

    region_snr <- .compute_region_snr_from_sed(
      flux_long = sed$flux_long,
      wavelength_range = wavelength_range,
      snr_stat = snr_stat,
      variance_inflation = variance_inflation
    )
    region_snr_list[[i]] <<- region_snr
    fits[[i]] <<- fit

    summary <- .finite_summary(region_snr)
    finite <- is.finite(region_snr)
    screen_snr <- if (snr_target$criterion == "min") summary$min else summary$median

    data.frame(
      Ncomp = k,
      actual_Ncomp = length(unique(stats::na.omit(as.vector(fit$cluster_map)))),
      min_region_snr = summary$min,
      median_region_snr = summary$median,
      q16_region_snr = summary$q16,
      q84_region_snr = summary$q84,
      screen_snr = screen_snr,
      all_regions_finite = all(finite),
      pass = is.finite(screen_snr) && screen_snr >= snr_target$threshold,
      error = NA_character_
    )
  })

  snr_grid <- do.call(rbind, rows)
  ok <- which(snr_grid$pass)
  if (!length(ok)) {
    stop("No SAGUI configuration satisfies the requested S/N threshold.")
  }

  chosen_i <- ok[1]
  chosen_k <- snr_grid$Ncomp[chosen_i]
  out <- fits[[as.character(chosen_k)]]
  out$Ncomp <- chosen_k
  out$target_snr <- snr_target$threshold
  out$snr_target <- snr_target$threshold
  out$min_snr <- if (snr_target$criterion == "min") snr_target$threshold else NULL
  out$median_snr <- if (snr_target$criterion == "median") snr_target$threshold else NULL
  out$region_snr <- region_snr_list[[as.character(chosen_k)]]
  out$snr_grid <- snr_grid
  out$snr_selection <- list(
    backend = backend,
    snr_stat = snr_stat,
    criterion = snr_target$criterion,
    requested_with = snr_target$label,
    wavelength_range = wavelength_range,
    variance_inflation = variance_inflation,
    error_fallback = if (is.null(var_cube)) error_fallback else "variance_cube"
  )
  out
}
