# Choose the Number of SAGUI Regions from a Target S/N

This helper runs SAGUI over a grid of candidate region counts and keeps
the finest segmentation satisfying a regional S/N cut. Use `min_snr` for
the common case where every region must exceed a minimum S/N. Unlike the
quick `cluster_snr` diagnostic returned by
[`segment_regions()`](https://rafaelsdesouza.com.br/sagui/reference/segment_regions.md),
this function can evaluate the S/N from the flux-conserving regional SED
table produced by
[`extract_region_sed()`](https://rafaelsdesouza.com.br/sagui/reference/extract_region_sed.md).

## Usage

``` r
choose_ncomp_by_snr(
  input,
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
  ...
)
```

## Arguments

- input:

  A 3-D array or FITS-like list with `imDat`.

- min_snr:

  Minimum regional S/N. This is the recommended user-facing argument.

- min_SNR:

  Alias for `min_snr`, accepted for readability.

- median_snr:

  Optional median regional S/N target. This is less strict than
  `min_snr`.

- target_snr:

  Deprecated compatibility alias. Use `min_snr` or `median_snr`.

- var_cube:

  Optional variance cube with the same dimensions as `input`.

- k_values:

  Candidate numbers of regions to test.

- wavelength_range:

  Optional numeric interval used to select wavelengths for the S/N
  calculation.

- band_values:

  Optional band labels or wavelengths passed to
  [`extract_region_sed()`](https://rafaelsdesouza.com.br/sagui/reference/extract_region_sed.md).

- backend:

  Clustering backend: `"exact"` calls
  [`segment_regions()`](https://rafaelsdesouza.com.br/sagui/reference/segment_regions.md)
  and `"sparse"` calls
  [`segment_regions_large()`](https://rafaelsdesouza.com.br/sagui/reference/segment_regions_large.md).

- snr_stat:

  Either integrated S/N across the selected wavelengths or the median
  per-wavelength S/N.

- screen_stat:

  Deprecated compatibility option used only with `target_snr`. Prefer
  `min_snr` or `median_snr`.

- variance_inflation:

  Multiplicative factor applied to propagated variances.

- error_fallback:

  Fallback error model used by
  [`extract_region_sed()`](https://rafaelsdesouza.com.br/sagui/reference/extract_region_sed.md)
  when no `var_cube` is provided.

- verbose:

  Logical; print progress messages.

- ...:

  Additional arguments passed to
  [`segment_regions()`](https://rafaelsdesouza.com.br/sagui/reference/segment_regions.md)
  or
  [`segment_regions_large()`](https://rafaelsdesouza.com.br/sagui/reference/segment_regions_large.md).

## Value

A SAGUI segmentation result for the selected number of regions. The
result also contains `region_snr`, `snr_grid`, `snr_target`, and
`snr_selection` entries.
