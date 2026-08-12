# Flatten a spectral cube into a pixel-by-band matrix

Flatten a spectral cube into a pixel-by-band matrix

## Usage

``` r
cube_to_matrix(x)
```

## Arguments

- x:

  3-D array or FITS-like list with `imDat`.

## Value

Numeric matrix with one row per spatial pixel and one column per band.
