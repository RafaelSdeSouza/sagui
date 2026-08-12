# Plot a segmented region map with dissolved Voronoi-style polygons

Plot a segmented region map with dissolved Voronoi-style polygons

## Usage

``` r
plot_region_map(
  cluster_data,
  palette = "magma",
  border_color = "black",
  border_linewidth = 0.8,
  background_color = "black"
)
```

## Arguments

- cluster_data:

  List containing a `cluster_map` matrix.

- palette:

  Viridis option name or a vector of colors.

- border_color:

  Boundary line color.

- border_linewidth:

  Boundary line width.

- background_color:

  Background fill color.

## Value

A `ggplot2` object.
