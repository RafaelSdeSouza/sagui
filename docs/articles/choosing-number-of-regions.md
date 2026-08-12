# Choose Region Count

The number of regions is an analysis parameter. A finer map is not
automatically more informative when regional photometric uncertainty is
large.

## S/N screen

[`choose_ncomp_by_snr()`](https://rafaelsdesouza.com.br/sagui/reference/choose_ncomp_by_snr.md)
evaluates candidate counts and keeps the finest segmentation satisfying
the requested regional signal-to-noise rule.

``` r
selected <- choose_ncomp_by_snr(
  input = cube,
  var_cube = var_cube,
  band_values = bands,
  k_values = c(8, 12, 16, 20, 24),
  min_snr = 10,
  backend = "exact",
  cluster_pretransform = "none"
)

selected$snr_grid
```

The S/N depends on the variance cube and wavelength range. A fallback
error estimate is useful for inspection, but it does not replace
calibrated uncertainties.

## Stability

Compare neighbouring candidate counts and check whether the structures
of interest persist.
