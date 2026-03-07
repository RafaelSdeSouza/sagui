#' Plot a region property map in Voronoi-style polygons
#'
#' @param cluster_data list containing `cluster_map`.
#' @param values named numeric vector of region values.
#' @param palette viridis option name or palette vector.
#' @param value_label legend title.
#' @param border_color polygon border color.
#' @param border_linewidth polygon border width.
#' @param background_color plot background color.
#' @param na_color fill color for missing values.
#' @param show_legend logical; whether to display the legend.
#' @param legend_position legend placement used when `show_legend` is `TRUE`.
#' @return A `ggplot2` object.
#' @export
plot_region_property_map <- function(cluster_data,
                                     values,
                                     palette = "magma",
                                     value_label = "Value",
                                     border_color = "black",
                                     border_linewidth = 0.8,
                                     background_color = "black",
                                     na_color = "transparent",
                                     show_legend = TRUE,
                                     legend_position = "right") {

  m <- cluster_data$cluster_map
  stopifnot(is.matrix(m))
  nr <- nrow(m); nc <- ncol(m)

  # ----------------------------
  # Values handling
  # ----------------------------
  if (is.null(names(values))) {
    names(values) <- as.character(seq_along(values))
  }

  # ----------------------------
  # Build pixel grid (x = Row, y = Col)
  # ----------------------------
  bbox <- sf::st_bbox(
    c(xmin = 0.5, ymin = 0.5, xmax = nr + 0.5, ymax = nc + 0.5),
    crs = NA
  )

  grid <- sf::st_make_grid(
    sf::st_as_sfc(bbox),
    n = c(nr, nc),
    what = "polygons"
  )

  df <- data.frame(
    Row     = rep(seq_len(nr), times = nc),
    Col     = rep(seq_len(nc), each  = nr),
    Cluster = as.vector(m)
  )

  sf_cells <- sf::st_sf(df, geometry = grid) |>
    dplyr::filter(!is.na(Cluster)) |>
    dplyr::mutate(
      Cluster = as.integer(Cluster),
      Value   = unname(values[as.character(Cluster)])
    )

  # ----------------------------
  # Dissolve pixels → regions (same as your style)
  # ----------------------------
  sf_regions <- sf_cells |>
    dplyr::group_by(Cluster) |>
    dplyr::summarise(
      Value = dplyr::first(Value),
      geometry = sf::st_union(geometry),
      .groups = "drop"
    )

  # ----------------------------
  # Single boundary network
  # ----------------------------
  boundary_sf <- sf::st_sf(
    geometry = sf::st_boundary(sf::st_union(sf_regions))
  )

  # ----------------------------
  # Palette handling for continuous values
  # ----------------------------
  if (is.character(palette) && length(palette) == 1) {
    fill_scale <- ggplot2::scale_fill_viridis_c(
      option = palette,
      name = value_label,
      na.value = na_color
    )
  } else if (is.character(palette) && length(palette) > 1) {
    fill_scale <- ggplot2::scale_fill_gradientn(
      colours = palette,
      name = value_label,
      na.value = na_color
    )
  } else {
    stop("palette must be a viridis option string or a vector of colors")
  }

  # ----------------------------
  # Plot
  # ----------------------------
  ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = sf_regions,
      ggplot2::aes(fill = Value),
      color = NA
    ) +
    ggplot2::geom_sf(
      data = boundary_sf,
      color = border_color,
      linewidth = border_linewidth,
      linejoin = "round",
      lineend  = "round"
    ) +
    fill_scale +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position  = if (show_legend) legend_position else "none",
      plot.background  = ggplot2::element_rect(fill = background_color, color = NA),
      panel.background = ggplot2::element_rect(fill = background_color, color = NA),
      legend.background = ggplot2::element_rect(fill = background_color, color = NA),
      legend.key        = ggplot2::element_rect(fill = background_color, color = NA),
      legend.text  = ggplot2::element_text(color = "white"),
      legend.title = ggplot2::element_text(color = "white")
    )
}
