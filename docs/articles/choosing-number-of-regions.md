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

The variance cube and wavelength range define the meaning of the S/N
threshold; preserve both with the result. A fallback error estimate is a
quick diagnostic, not a replacement for calibrated uncertainty
propagation.

## Stability

Inspect neighbouring candidate counts and check whether the structures
used in the scientific interpretation persist. Record the accepted rule
before examining downstream model fits where practical.

The real Sagui-10 target-S/N products are listed in [Paper
Examples](https://rafaelsdesouza.com.br/sagui/articles/paper-examples-reproduction.html#target-sn-selection).
