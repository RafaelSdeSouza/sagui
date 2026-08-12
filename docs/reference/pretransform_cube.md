# Apply a spectral pretransform to a full cube

Apply a spectral pretransform to a full cube

## Usage

``` r
pretransform_cube(x, method = "none")
```

## Arguments

- x:

  3-D array or FITS-like list with `imDat`.

- method:

  Either a built-in transform name or a custom function accepted by
  [`pretransform_spectra()`](https://rafaelsdesouza.com.br/sagui/reference/pretransform_spectra.md).

## Value

A 3-D numeric array with the same dimensions as the input cube.
