# Segment photometric regions from large spectral cubes

Large-cube version of
[`segment_regions()`](https://rafaelsdesouza.com.br/sagui/reference/segment_regions.md).
It keeps the same photometric mask, pretransform, and output structure,
but replaces exact all-pairs Ward clustering with a kNN-restricted
sparse-Ward approximation.

## Usage

``` r
segment_regions_large(
  input,
  Ncomp = 5,
  redshift = 0,
  pretransform = NULL,
  mask_pretransform = "none",
  cluster_pretransform = "none",
  scale_fn = median_scale,
  n_regions = NULL,
  use_starlet_mask = TRUE,
  support_method = c("starlet", "adaptive", "starlet_contourlet", "contourlet"),
  support_args = list(),
  collapse_fn = collapse_white_light,
  starlet_J = 5,
  starlet_scales = 2:5,
  include_coarse = FALSE,
  denoise_k = 2.5,
  mode = c("soft", "hard"),
  positive_only = TRUE,
  clean_mask = FALSE,
  min_mask_area = 1L,
  close_size = 1L,
  open_size = 1L,
  keep_largest = FALSE,
  mask_mode = c("na", "zero"),
  hclust_method = "ward.D2",
  knn_k = 40,
  auto_k = FALSE,
  max_k = NULL,
  feature_scale = c("none", "robust_col"),
  spatial_weight = 0,
  return_details = FALSE,
  verbose = TRUE
)
```

## Arguments

- input:

  3-D array or FITS-like list with `imDat`.

- Ncomp:

  Number of output regions.

- redshift:

  Reserved compatibility argument carried over from `capivara`.
  Currently unused.

- pretransform:

  Deprecated alias for `cluster_pretransform`.

- mask_pretransform:

  Optional spectral pretransform applied before the white-light collapse
  and starlet mask. This lets the segmentation mask be built on a
  transformed photometric representation while leaving the clustering
  stage untouched.

- cluster_pretransform:

  Optional spectral pretransform applied to the valid spectra matrix
  before row-wise scaling and clustering. May be one of `"none"`,
  `"asinh"`, `"log1p"`, `"signed_log1p"`, `"copula_uniform"`,
  `"copula_gaussian"`, or a custom function returning a matrix with the
  same dimensions as the input.

- scale_fn:

  Per-spectrum scaling function applied row-wise before clustering. Use
  `identity` for no row scaling.

- n_regions:

  Deprecated alias for `Ncomp`.

- use_starlet_mask:

  Logical; if `TRUE`, derive a photometric mask before clustering.

- support_method:

  Foreground support builder used when `use_starlet_mask = TRUE`.
  Options are `"starlet"` (default), `"adaptive"` (multi-band
  statistical support), `"starlet_contourlet"` (directional refinement),
  and `"contourlet"` (standalone directional support; more
  experimental).

- support_args:

  Optional named list passed to the selected support builder. For
  `"adaptive"`, arguments are passed to
  [`build_adaptive_support()`](https://rafaelsdesouza.com.br/sagui/reference/build_adaptive_support.md).
  For contourlet methods, arguments are passed to
  [`build_contourlet_mask()`](https://rafaelsdesouza.com.br/sagui/reference/build_contourlet_mask.md).

- collapse_fn:

  Function used to collapse the cube to a 2-D image.

- starlet_J:

  Number of starlet scales.

- starlet_scales:

  Scales to keep when reconstructing the starlet image.

- include_coarse:

  Logical; include the coarse plane in the starlet reconstruction.

- denoise_k:

  Optional starlet denoising threshold.

- mode:

  Starlet thresholding mode.

- positive_only:

  Logical; keep only positive reconstruction values in the mask.

- clean_mask:

  Logical; apply optional support-mask cleanup after the foreground
  support is built.

- min_mask_area:

  Minimum connected-component area retained when `clean_mask = TRUE`.

- close_size:

  Binary closing brush size used when `clean_mask = TRUE`. Use `1` to
  skip closing.

- open_size:

  Binary opening brush size used when `clean_mask = TRUE`. Use `1` to
  skip opening.

- keep_largest:

  Logical; keep only the largest connected component after area
  filtering when `clean_mask = TRUE`.

- mask_mode:

  Mask fill mode passed to
  [`mask_cube()`](https://rafaelsdesouza.com.br/sagui/reference/mask_cube.md).

- hclust_method:

  Linkage method passed to
  [`hclust()`](https://rdrr.io/r/stats/hclust.html).

- knn_k:

  Number of nearest neighbours for the sparse-Ward graph.

- auto_k:

  Logical; if `TRUE`, increase `knn_k` when the graph remains
  disconnected before the requested number of regions is reached.

- max_k:

  Maximum `k` allowed when `auto_k = TRUE`.

- feature_scale:

  Optional column-wise feature scaling after row scaling.

- spatial_weight:

  Optional weight for appending normalized x/y pixel coordinates to the
  clustering features.

- return_details:

  Logical; include valid indices, features, and labels.

- verbose:

  Logical; print sparse-Ward progress messages.

## Value

A segmentation result list compatible with
[`segment_regions()`](https://rafaelsdesouza.com.br/sagui/reference/segment_regions.md)
and downstream photometry helpers such as
[`RegionPhotometry()`](https://rafaelsdesouza.com.br/sagui/reference/RegionPhotometry.md).
