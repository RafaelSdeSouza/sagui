# Asinh-stretch an image or numeric array

Asinh-stretch an image or numeric array

## Usage

``` r
asinh_stretch(
  x,
  qlo = 0.001,
  qhi = 0.999,
  scale = NULL,
  nonneg = TRUE,
  na.rm = TRUE
)
```

## Arguments

- x:

  Numeric vector, matrix, or array.

- qlo:

  Lower quantile used as the display floor.

- qhi:

  Upper quantile used as the display ceiling.

- scale:

  Optional asinh scale. If `NULL`, it is estimated from the quantile
  range.

- nonneg:

  Logical; clamp negative values to zero before stretching.

- na.rm:

  Logical; ignore non-finite values when estimating quantiles.

## Value

`x` transformed to an asinh-stretched numeric object with the same
dimensions.
