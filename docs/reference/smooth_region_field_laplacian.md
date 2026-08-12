# Smooth region values over an adjacency graph

Smooth region values over an adjacency graph

## Usage

``` r
smooth_region_field_laplacian(
  region_id_mat,
  region_values,
  outside_ids = c(NA_integer_, 0L),
  adjacency = c("rook", "queen"),
  lambda = 5,
  keep_na_outside = TRUE
)
```

## Arguments

- region_id_mat:

  matrix of region identifiers.

- region_values:

  named numeric vector or data frame with `region_id` and `value`.

- outside_ids:

  unused compatibility placeholder for non-region values.

- adjacency:

  neighborhood definition: `"rook"` or `"queen"`.

- lambda:

  smoothing strength.

- keep_na_outside:

  logical; keep non-region pixels as `NA`.

## Value

A list with the interpolated pixel matrix and region-level predictions.
