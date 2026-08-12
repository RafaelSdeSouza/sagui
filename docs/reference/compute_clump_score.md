# Compute a compact SED-clump score

Builds a score map for compact structures whose photometric SEDs differ
from the local surrounding galaxy light. The score combines a
local-contrast term from a collapsed image with a local SED-anomaly term
measured relative to a smoothed version of the cube.

## Usage

``` r
compute_clump_score(
  input,
  support = NULL,
  bands = NULL,
  contrast_weight = 0.5,
  sed_weight = 0.5,
  small_sigma = 0.7,
  large_sigma = 4
)
```

## Arguments

- input:

  3-D array or FITS-like list with `imDat`.

- support:

  Optional logical foreground support mask. If `NULL`, all finite pixels
  in the collapsed image are used.

- bands:

  Bands used to build the score. `NULL` uses all bands.

- contrast_weight:

  Weight assigned to the local-contrast term.

- sed_weight:

  Weight assigned to the local SED-anomaly term.

- small_sigma:

  Gaussian smoothing scale for the compact image term.

- large_sigma:

  Gaussian smoothing scale for the local-background term.

## Value

A list with `score`, `contrast`, `sed_anomaly`, and `collapsed`.
