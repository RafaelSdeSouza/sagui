#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(grid)
  library(png)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1] %||% file.path(
  "/Users/rd23aag/Documents/GitHub/crp8_segmentation",
  "results/figures/paper_repro/target_snr_selection/sagui10"
)
output_pdf <- args[2] %||% file.path(input_dir, "sagui10_target_snr_static.pdf")

targets <- c(5, 10, 15, 20)
files <- file.path(input_dir, sprintf("sagui10_target_snr_%02d.png", targets))
missing <- files[!file.exists(files)]
if (length(missing) > 0) {
  stop("Missing target-S/N PNG(s):\n", paste(missing, collapse = "\n"), call. = FALSE)
}

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

panel_w <- 1 / length(files)
title_h <- 0.24
image_h <- 0.76

pdf(output_pdf, width = 8.2, height = 2.65, useDingbats = FALSE, family = "Times")
grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.965, height = 0.92))

for (i in seq_along(files)) {
  x0 <- (i - 1) * panel_w
  xc <- x0 + panel_w / 2

  grid.text(
    label = sprintf("S/N = %d", targets[i]),
    x = unit(xc, "npc"),
    y = unit(image_h + title_h * 0.62, "npc"),
    gp = gpar(fontsize = 24, col = "black", fontfamily = "Times", fontface = "italic")
  )

  img <- png::readPNG(files[i])
  pushViewport(viewport(
    x = unit(xc, "npc"),
    y = unit(image_h / 2, "npc"),
    width = unit(panel_w, "npc"),
    height = unit(image_h, "npc")
  ))
  grid.raster(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
  grid.rect(gp = gpar(fill = NA, col = "grey35", lwd = 0.55))
  popViewport()
}

popViewport()
dev.off()

message("Wrote ", output_pdf)
