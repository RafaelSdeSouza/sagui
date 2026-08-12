# Build a contourlet-like directional support mask from a spectral cube

This function extends the default starlet support with a directional,
multiscale arm-sensitive support map. The implementation is
contourlet-like rather than a full contourlet transform: it uses a bank
of oriented, elongated Gaussian filters across multiple scales, combines
them into a directional response map, and optionally grows the starlet
support into nearby arm-like pixels.

## Usage

``` r
build_contourlet_mask(
  input,
  collapse_fn = collapse_white_light,
  pretransform = "none",
  combine = c("hybrid", "standalone"),
  base_mask = NULL,
  starlet_J = 5,
  starlet_scales = 2:5,
  include_coarse = FALSE,
  denoise_k = 2.5,
  mode = c("soft", "hard"),
  positive_only = TRUE,
  angles_deg = seq(0, 165, by = 15),
  scales = list(c(4.5, 1.2), c(7.5, 1.8), c(11, 2.6)),
  background_sigma = 8,
  score_low_q = 0.92,
  score_high_q = 0.98,
  contrast_low_q = 0.8,
  contrast_high_q = 0.93,
  max_growth_iter = 8L,
  close_size = 3L,
  open_size = 3L
)
```

## Arguments

- input:

  3-D array or FITS-like list with `imDat`.

- collapse_fn:

  Function used to collapse the cube to a 2-D image.

- pretransform:

  Optional spectral pretransform applied to the cube before the
  white-light collapse.

- combine:

  Either `"hybrid"` (recommended; start from the starlet support and
  extend it directionally) or `"standalone"` (use only the directional
  support).

- base_mask:

  Optional precomputed base mask. When `combine = "hybrid"` and
  `base_mask` is `NULL`, the base support is built from
  [`build_starlet_mask()`](https://rafaelsdesouza.com.br/sagui/reference/build_starlet_mask.md).

- starlet_J:

  Number of starlet scales used when building the hybrid base support.

- starlet_scales:

  Scales to reconstruct for the hybrid base support.

- include_coarse:

  Logical; include the coarse plane in the hybrid base support
  reconstruction.

- denoise_k:

  Optional starlet denoising threshold used in the hybrid base support.

- mode:

  Thresholding mode for the hybrid base support reconstruction.

- positive_only:

  Logical; keep only positive reconstruction values in the hybrid base
  support.

- angles_deg:

  Vector of filter orientations in degrees.

- scales:

  List of `(sigma_major, sigma_minor)` pairs defining the directional
  filter bank.

- background_sigma:

  Smoothing scale used to estimate the local background before the
  directional response is measured.

- score_low_q:

  Lower quantile used to define the allowed directional growth region.

- score_high_q:

  Upper quantile used to define strong directional seeds.

- contrast_low_q:

  Lower quantile used to define the allowed contrast growth region.

- contrast_high_q:

  Upper quantile used to define strong contrast seeds.

- max_growth_iter:

  Maximum number of constrained growth iterations.

- close_size:

  Size of the binary closing brush.

- open_size:

  Size of the binary opening brush.

## Value

A list containing the collapsed image, the optional base support, the
final logical mask, the directional score map, the local contrast map,
and threshold metadata.
