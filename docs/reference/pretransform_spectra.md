# Apply a spectral pretransform before segmentation

Applies a column-wise pretransform to a spectra matrix before row-wise
scaling and clustering. This is useful for benchmarking clustering
behavior under simple nonlinear transforms such as
[`asinh()`](https://rdrr.io/r/base/Hyperbolic.html),
[`log1p()`](https://rdrr.io/r/base/Log.html), or rank-based copula
mappings.

## Usage

``` r
pretransform_spectra(
  x,
  method = c("none", "asinh", "log1p", "signed_log1p", "copula_uniform",
    "copula_gaussian")
)
```

## Arguments

- x:

  Numeric matrix with one spectrum per row and one band per column.

- method:

  Either a string naming a built-in transform or a function that takes a
  numeric matrix and returns a numeric matrix with the same dimensions.
  Built-in options are `"none"`, `"asinh"`, `"log1p"`, `"signed_log1p"`,
  `"copula_uniform"`, and `"copula_gaussian"`.

## Value

A numeric matrix with the same dimensions as `x`.
