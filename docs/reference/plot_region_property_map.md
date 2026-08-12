# Plot a region property map in Voronoi-style polygons

Plot a region property map in Voronoi-style polygons

## Usage

``` r
plot_region_property_map(
  cluster_data,
  values,
  palette = "magma",
  value_label = "Value",
  border_color = "black",
  border_linewidth = 0.8,
  background_color = "black",
  na_color = "transparent",
  show_legend = TRUE,
  legend_position = "right"
)
```

## Arguments

- cluster_data:

  list containing `cluster_map`.

- values:

  named numeric vector of region values.

- palette:

  viridis option name or palette vector.

- value_label:

  legend title.

- border_color:

  polygon border color.

- border_linewidth:

  polygon border width.

- background_color:

  plot background color.

- na_color:

  fill color for missing values.

- show_legend:

  logical; whether to display the legend.

- legend_position:

  legend placement used when `show_legend` is `TRUE`.

## Value

A `ggplot2` object.
