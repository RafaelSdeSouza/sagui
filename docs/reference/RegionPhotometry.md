# Region-level photometry on hyperspectral cubes (no background subtraction)

Aggregates flux by region and band for a hyperspectral or IFU cube.
Variance-aware uncertainties are used when available, with simple
fallback error models otherwise. A numeric `lambda` column is also added
when band values can be interpreted as wavelengths.

## Usage

``` r
RegionPhotometry(
  cube,
  labels,
  bkg = NULL,
  var_cube = NULL,
  sigma_band = NULL,
  band_values = NULL,
  digits_lambda_colnames = 6,
  return_painted_cube = FALSE,
  error_fallback = c("none", "flux_over_sqrt_n", "poisson", "mad_sky")
)
```

## Arguments

- cube:

  A 3-D array with dimensions `[nx, ny, nb]` or a FITS-like list with
  `imDat`.

- labels:

  A matrix `[nx, ny]` with positive region identifiers.

- bkg:

  Deprecated compatibility argument. Ignored.

- var_cube:

  Variance cube with the same shape as `cube`, or a FITS-like list.

- sigma_band:

  Per-band standard deviation vector used when `var_cube` is `NULL`.

- band_values:

  Optional band labels or wavelengths. If omitted and `cube` is
  FITS-like,
  [`FITSio::axVec()`](https://rdrr.io/pkg/FITSio/man/axVec.html) is
  used.

- digits_lambda_colnames:

  Significant digits used when numeric band values are converted to
  column names.

- return_painted_cube:

  Logical; if `TRUE`, also return a cube painted with region fluxes.

- error_fallback:

  Fallback error model used when no variance cube is provided.

## Value

A list with `flux_long`, `flux_wide`, `painted_cube`, `bands`, and
`band_names`.
