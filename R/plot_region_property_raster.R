#' Plot a region property map on a pixel grid
#'
#' @param cluster_map matrix of region identifiers.
#' @param values named numeric vector of region values.
#' @param palette viridis option name or palette vector.
#' @param value_label legend title.
#' @param na_color fill color for missing values.
#' @param show_legend logical; whether to display the legend.
#' @return A `ggplot2` object.
#' @export
plot_region_property_raster <- function(cluster_map,
                                        values,
                                        palette = "magma",
                                        value_label = "Value",
                                        na_color = "white",
                                        show_legend = TRUE) {
  stopifnot(is.matrix(cluster_map), is.numeric(values))

  # If values is not named, assume regions are 1:length(values)
  if (is.null(names(values))) {
    names(values) <- as.character(seq_along(values))
  }
  
  # Flatten cluster_map to a vector of region IDs
  region_vec <- as.vector(cluster_map)
  
  # Match region IDs to values (handle NA)
  value_vec <- rep(NA_real_, length(region_vec))
  ok <- !is.na(region_vec)
  value_vec[ok] <- values[as.character(region_vec[ok])]
  
  # Build long data frame
  n_row <- nrow(cluster_map)
  n_col <- ncol(cluster_map)
  
  prop_df <- tibble::tibble(
    Row   = rep(seq_len(n_row), times = n_col),
    Col   = rep(seq_len(n_col), each  = n_row),
    Value = value_vec
  )

  fill_scale <- if (is.character(palette) && length(palette) == 1) {
    viridis::scale_fill_viridis(option = palette, na.value = na_color, name = value_label)
  } else if (is.character(palette) && length(palette) > 1) {
    ggplot2::scale_fill_gradientn(colours = palette, na.value = na_color, name = value_label)
  } else {
    stop("palette must be a viridis option string or a vector of colors")
  }

  ggplot2::ggplot(prop_df, ggplot2::aes(x = Row, y = Col, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    fill_scale +
    ggplot2::theme_void() +
    ggplot2::theme(
      axis.text        = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_blank(),
      legend.position  = if (show_legend) "right" else "none"
    )
}
