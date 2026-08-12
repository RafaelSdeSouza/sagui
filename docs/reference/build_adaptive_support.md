# Build an adaptive multi-band support mask

Builds a foreground support from band-wise detection evidence rather
than from a single collapsed white image. Each selected band is
optionally transformed, converted to a robust sky-normalized
significance map, and then thresholded. The final support keeps pixels
detected in at least `min_band_persistence` bands, with an optional
high-significance single-band rescue channel.

This is intended for difficult broadband/medium-band cubes where a
single white-light starlet mask can either cut diffuse emission or keep
too much structured background.

## Usage

``` r
build_adaptive_support(
  input,
  bands = NULL,
  pretransform = "none",
  transform = c("none", "asinh", "signed_log1p", "log1p", "copula_uniform",
    "copula_gaussian"),
  sky_method = c("border", "all"),
  border_fraction = 0.1,
  z_threshold = 3,
  min_band_persistence = NULL,
  single_band_z = Inf,
  smooth_sigma = 0
)
```

## Arguments

- input:

  3-D array or FITS-like list with `imDat`.

- bands:

  Bands used to build the support. May be `NULL` for all bands, a
  numeric vector of band indices, or a character vector matching cube
  dimnames along the third dimension.

- pretransform:

  Optional cube pretransform applied before support construction. Uses
  [`pretransform_cube()`](https://rafaelsdesouza.com.br/sagui/reference/pretransform_cube.md).

- transform:

  Per-band image transform before robust sky normalization. Options are
  `"none"`, `"asinh"`, `"signed_log1p"`, `"log1p"`, `"copula_uniform"`,
  and `"copula_gaussian"`.

- sky_method:

  Pixels used for robust sky/noise estimation. `"border"` uses image
  borders; `"all"` uses all finite pixels.

- border_fraction:

  Fraction of rows/columns used as image border when
  `sky_method = "border"`.

- z_threshold:

  Band-wise detection threshold in robust sigma units.

- min_band_persistence:

  Minimum number of bands above `z_threshold`. If `NULL`, defaults to
  two bands when possible.

- single_band_z:

  Optional single-band rescue threshold. Pixels above this robust
  significance in any selected band are kept even if they fail
  `min_band_persistence`. Use `Inf` to disable.

- smooth_sigma:

  Optional Gaussian smoothing applied to each significance map before
  thresholding. Requires `EBImage` or `imager`.

## Value

A list with `collapsed`, `reconstruction`, `mask`, `evidence`,
`evidence_count`, `z_maps`, `band_stats`, and metadata.
