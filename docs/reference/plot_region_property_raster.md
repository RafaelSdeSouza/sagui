# Plot a region property map on a pixel grid

Plot a region property map on a pixel grid

## Usage

``` r
plot_region_property_raster(
  cluster_map,
  values,
  palette = "magma",
  value_label = "Value",
  na_color = "white",
  show_legend = TRUE
)
```

## Arguments

- cluster_map:

  matrix of region identifiers.

- values:

  named numeric vector of region values.

- palette:

  viridis option name or palette vector.

- value_label:

  legend title.

- na_color:

  fill color for missing values.

- show_legend:

  logical; whether to display the legend.

## Value

A `ggplot2` object.
