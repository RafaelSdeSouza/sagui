# Reconstruct and threshold starlet scales to produce a region mask

Reconstruct and threshold starlet scales to produce a region mask

## Usage

``` r
threshold_starlet_regions(
  dec,
  keep_scales = 2:5,
  include_coarse = FALSE,
  denoise_k = 3,
  mode = c("soft", "hard"),
  threshold = c("mad", "abs"),
  k = 3.5,
  tau = NULL,
  positive_only = TRUE,
  per_scale_positive = TRUE,
  k_hi = 5,
  k_lo = 3,
  area_min = 12,
  keep_negatives = FALSE,
  ring_nbins = 8
)
```

## Arguments

- dec:

  `sagui_starlet` object (or a matrix; if matrix, J inferred from
  `keep_scales`)

- keep_scales:

  integer vector (e.g., 3:6)

- include_coarse:

  logical

- denoise_k:

  NULL or numeric (pre-recon per-scale)

- mode:

  "soft" or "hard" (per-scale)

- threshold:

  "mad" for k·MAD, or "abs" for absolute tau

- k:

  numeric; multiplier for MAD (unused; kept for compat)

- tau:

  numeric; absolute threshold when threshold == "abs"

- positive_only:

  logical; TRUE uses rec\>thr; FALSE uses \|rec\|\>thr

- per_scale_positive:

  logical; clamp each reconstructed starlet scale to positive values
  before thresholding.

- k_hi:

  numeric; high MAD threshold used to define seed pixels.

- k_lo:

  numeric; low MAD threshold used to define candidate pixels.

- area_min:

  integer; minimum connected-component area to keep.

- keep_negatives:

  logical; retain negative reconstruction values.

- ring_nbins:

  integer; number of radial bins for background MAD.

## Value

list(soft_rec, hard_rec, mask, sigma, threshold, seeds, candidate)
