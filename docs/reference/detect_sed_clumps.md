# Detect compact SED-clump seeds

Detect compact SED-clump seeds

## Usage

``` r
detect_sed_clumps(
  score,
  support = NULL,
  n_clumps = NULL,
  max_clumps = 50L,
  score_quantile = 0.995,
  min_distance = 3
)
```

## Arguments

- score:

  2-D clump-score map.

- support:

  Optional logical support mask.

- n_clumps:

  Optional exact number of seeds to keep. If `NULL`, seeds are selected
  above `score_quantile` up to `max_clumps`.

- max_clumps:

  Maximum number of seeds when `n_clumps = NULL`.

- score_quantile:

  Quantile threshold applied to finite support pixels.

- min_distance:

  Minimum Euclidean distance, in pixels, between seeds.

## Value

A data frame with `clump_id`, `row`, `col`, and `score`.
