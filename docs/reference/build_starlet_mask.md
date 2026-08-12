# Build a starlet-based segmentation mask from a spectral cube

Build a starlet-based segmentation mask from a spectral cube

## Usage

``` r
build_starlet_mask(
  input,
  collapse_fn = collapse_white_light,
  pretransform = "none",
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
  keep_largest = FALSE
)
```

## Arguments

- input:

  3-D array or FITS-like list with `imDat`.

- collapse_fn:

  Function used to collapse the cube to a 2-D image.

- pretransform:

  Optional spectral pretransform applied to the cube before the
  white-light collapse and starlet decomposition. This is useful when
  the photometric mask should be derived from a transformed
  representation rather than the original flux cube.

- starlet_J:

  Number of starlet scales.

- starlet_scales:

  Scales to reconstruct.

- include_coarse:

  Logical; include the coarse plane in the reconstruction.

- denoise_k:

  Optional starlet denoising threshold.

- mode:

  Thresholding mode for starlet reconstruction.

- positive_only:

  Logical; keep only positive reconstruction values in the mask.

- clean_mask:

  Logical; apply optional support-mask cleanup.

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

## Value

A list with the collapsed image, decomposition, reconstruction, and
logical mask.
