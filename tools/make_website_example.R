#!/usr/bin/env Rscript

if (!requireNamespace("sagui", quietly = TRUE)) {
  stop("Install sagui before running this script.")
}
if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("cowplot", quietly = TRUE)) {
  stop("This website figure requires ggplot2 and cowplot.")
}

set.seed(42)

bands <- c(0.90, 1.15, 1.50, 1.82, 2.00, 2.77, 3.56, 4.10, 4.44)
nx <- ny <- 24L
xy <- expand.grid(x = seq_len(nx), y = seq_len(ny))
dx <- xy$x - 12.5
dy <- xy$y - 12.5
radius <- sqrt((dx / 1.10)^2 + (dy / 0.78)^2)
support <- matrix((dx / 11)^2 + (dy / 8.3)^2 <= 1, nx, ny)

profiles <- rbind(
  bulge = c(0.58, 0.78, 1.02, 1.16, 1.14, 0.98, 0.82, 0.72, 0.66),
  disc  = c(0.82, 0.94, 1.04, 1.08, 1.05, 0.92, 0.84, 0.78, 0.74),
  knot1 = c(1.18, 1.11, 1.02, 0.96, 0.93, 0.88, 0.98, 0.91, 0.86),
  knot2 = c(1.08, 1.04, 0.99, 0.95, 0.94, 0.91, 1.05, 0.99, 0.93)
)

weights <- cbind(
  exp(-(radius / 3.0)^2),
  exp(-radius / 8.5),
  exp(-((dx - 4.5)^2 + (dy + 1.8)^2) / 5),
  exp(-((dx + 4.0)^2 + (dy - 2.6)^2) / 6)
)
weights <- weights / rowSums(weights)
surface_brightness <- 0.12 + 1.35 * exp(-radius / 8)
flux <- surface_brightness * (weights %*% profiles)
flux <- flux + matrix(stats::rnorm(length(flux), sd = 0.012), nrow(flux))
cube <- array(flux, dim = c(nx, ny, length(bands)))
for (b in seq_along(bands)) cube[, , b][!support] <- NA_real_

seg <- sagui::segment_regions(
  input = cube,
  Ncomp = 6,
  use_starlet_mask = FALSE,
  cluster_pretransform = "none"
)
regional <- sagui::extract_region_sed(
  cube = cube,
  labels = seg$cluster_map,
  band_values = bands,
  sigma_band = rep(0.012, length(bands))
)

palette <- c("#173f66", "#277da1", "#43aa8b", "#f2bd43", "#ec6a5c", "#b83b4a")
map_df <- transform(xy, region = factor(as.vector(seg$cluster_map)))
support_df <- transform(xy, support = as.numeric(as.vector(support)))

p_map <- ggplot2::ggplot(map_df, ggplot2::aes(x, y, fill = region)) +
  ggplot2::geom_tile(width = 1, height = 1) +
  ggplot2::geom_contour(
    data = support_df,
    ggplot2::aes(x = x, y = y, z = support),
    breaks = 0.5,
    colour = "#17233b",
    linewidth = 0.75,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_fill_manual(values = palette, na.value = "white", drop = FALSE) +
  ggplot2::coord_equal(expand = FALSE) +
  ggplot2::labs(title = "Photometric regions", subtitle = "Categorical labels", x = NULL, y = NULL) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(colour = "#17233b", face = "bold", size = 15),
    plot.subtitle = ggplot2::element_text(colour = "#5c6878", size = 10),
    legend.position = "none",
    plot.margin = ggplot2::margin(8, 10, 8, 8)
  )

sed <- regional$flux_long
anchor <- sed$flux[sed$lambda == 2.00]
names(anchor) <- sed$region[sed$lambda == 2.00]
normaliser <- unname(anchor[as.character(sed$region)])
sed$normalised_flux <- sed$flux / normaliser
sed$normalised_error <- sed$flux_err / normaliser
sed$region <- factor(sed$region)

p_sed <- ggplot2::ggplot(
  sed,
  ggplot2::aes(lambda, normalised_flux, colour = region, group = region)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = normalised_flux - normalised_error,
      ymax = normalised_flux + normalised_error,
      fill = region
    ),
    alpha = 0.10,
    colour = NA,
    show.legend = FALSE
  ) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_point(size = 2.5, colour = "white") +
  ggplot2::geom_point(size = 1.6) +
  ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
  ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
  ggplot2::scale_x_continuous(breaks = c(0.9, 1.5, 2, 2.8, 3.6, 4.4)) +
  ggplot2::labs(
    title = "Regional SEDs",
    subtitle = "Integrated flux, normalised at 2 µm",
    x = "Observed wavelength (µm)",
    y = expression(F[nu] / F[nu] * " at 2 µm"),
    colour = "Region"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.line = ggplot2::element_line(colour = "#17233b", linewidth = 0.5),
    axis.text = ggplot2::element_text(colour = "#33445a"),
    axis.title = ggplot2::element_text(colour = "#17233b"),
    legend.position = "top",
    legend.justification = "left",
    legend.title = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(colour = "#17233b", face = "bold", size = 15),
    plot.subtitle = ggplot2::element_text(colour = "#5c6878", size = 10),
    plot.margin = ggplot2::margin(8, 8, 8, 10)
  )

combined <- cowplot::plot_grid(p_map, p_sed, nrow = 1, rel_widths = c(0.82, 1.42))
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
  "man/figures/sagui-first-run-synthetic-sed.png",
  combined,
  width = 10,
  height = 4.35,
  dpi = 220,
  bg = "white"
)

dir.create("docs/website-audit", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  regional$flux_long,
  "docs/website-audit/sagui-first-run-regional-seds.csv",
  row.names = FALSE
)

message("Wrote the synthetic website figure and regional SED table.")
