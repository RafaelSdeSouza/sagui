#' Plot a segmented region map with dissolved Voronoi-style polygons
#'
#' @param cluster_data List containing a `cluster_map` matrix.
#' @param palette Viridis option name or a vector of colors.
#' @param border_color Boundary line color.
#' @param border_linewidth Boundary line width.
#' @param background_color Background fill color.
#' @return A `ggplot2` object.
#' @export
plot_region_map <- function(cluster_data,
                            palette = "magma",
                            border_color = "black",
                            border_linewidth = 0.8,
                            background_color = "black") {

  m <- cluster_data$cluster_map
  stopifnot(is.matrix(m))
  nr <- nrow(m); nc <- ncol(m)

  # ----------------------------
  # Palette handling
  # ----------------------------
  cl_levels <- sort(unique(as.integer(m[!is.na(m)])))
  n_clusters <- length(cl_levels)

  if (is.character(palette) && length(palette) == 1) {
    colors <- viridis::viridis(n_clusters, option = palette)
  } else if (is.character(palette) && length(palette) > 1) {
    colors <- if (length(palette) < n_clusters)
      grDevices::colorRampPalette(palette)(n_clusters)
    else
      palette[1:n_clusters]
  } else {
    stop("palette must be a viridis option string or a vector of colors")
  }
  names(colors) <- as.character(cl_levels)

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
    dplyr::mutate(Cluster = factor(Cluster))

  # ----------------------------
  # Dissolve pixels → regions
  # ----------------------------
  sf_regions <- sf_cells |>
    dplyr::group_by(Cluster) |>
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")


#   sf_regions <- sf_regions |>
#    dplyr::mutate(geometry = sf::st_simplify(geometry, dTolerance = 0.5))

  # ----------------------------
  # Single boundary network
  # ----------------------------
  boundary_sf <- sf::st_sf(
    geometry = sf::st_boundary(sf::st_union(sf_regions))
  )

  # ----------------------------
  # Plot
  # ----------------------------
  ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = sf_regions,
      ggplot2::aes(fill = Cluster),
      color = NA
    ) +
    ggplot2::geom_sf(
      data = boundary_sf,
      color = border_color,
      linewidth = border_linewidth,
      linejoin = "round",
      lineend  = "round"
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.background  = ggplot2::element_rect(fill = background_color, color = NA),
      panel.background = ggplot2::element_rect(fill = background_color, color = NA)
    )
}
