# sagui <img align="right" src="images/sagui_upscaled.png" width="190" alt="sagui logo">

Photometric SED-based segmentation for IFU data cubes.

## Website

The package website is the main documentation entry point:

- [sagui website](https://rafaelsdesouza.github.io/sagui/)

## Installation

```r
install.packages("remotes")
remotes::install_github("RafaelSdeSouza/sagui")
library(sagui)
```

## Minimal example

```r
library(sagui)
library(FITSio)

x <- FITSio::readFITS("manga-8135-12701-LOGCUBE.fits")

seg <- segment_regions(
  input = x,
  Ncomp = 30,
  use_starlet_mask = TRUE,
  starlet_J = 5,
  starlet_scales = 2:5,
  pretransform = "asinh"
)

sed <- extract_region_sed(
  cube = x,
  labels = seg$cluster_map,
  band_values = FITSio::axVec(3, x$axDat)
)
```

For larger cubes where exact Ward clustering becomes memory-limited, use the
sparse-Ward approximation with the same photometric mask and output structure:

```r
seg <- segment_regions_large(
  input = x,
  Ncomp = 30,
  use_starlet_mask = TRUE,
  pretransform = "asinh",
  knn_k = 40,
  scale_fn = identity
)
```

`segment_regions_large()` mirrors `segment_regions()` but avoids constructing
the full all-pairs distance matrix. Use `scale_fn = identity` for raw
photometry, or `scale_fn = median_scale` to match the default exact workflow.

## Scope

`sagui` focuses on:

- photometric masking with starlet reconstruction
- region segmentation of IFU cubes
- sparse-Ward segmentation for large photometric cubes
- optional spectral pretransforms for clustering benchmarks
- integrated region SEDs with uncertainties
- region-level visualization for downstream fitting
