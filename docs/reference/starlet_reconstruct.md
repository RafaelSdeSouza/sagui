# Reconstruct from a subset of starlet scales

Reconstruct from a subset of starlet scales

## Usage

``` r
starlet_reconstruct(
  dec,
  keep_scales = c(2, 3),
  include_coarse = FALSE,
  denoise_k = NULL,
  mode = c("soft", "hard"),
  na_policy = c("preserve", "zero")
)
```

## Arguments

- dec:

  `sagui_starlet` object

- keep_scales:

  integer vector of scales to include (e.g., 3:6)

- include_coarse:

  logical; add coarse plane `cJ`

- denoise_k:

  NULL or numeric (k·sigma MAD per-plane)

- mode:

  "soft" or "hard" thresholding

- na_policy:

  "preserve" or "zero" for NA handling during sum
