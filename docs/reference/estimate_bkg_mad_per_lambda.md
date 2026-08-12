# Robust per-wavelength background and scatter

Computes the spatial median (background) and MAD-based scatter per
wavelength slice of a 3-D cube.

## Usage

``` r
estimate_bkg_mad_per_lambda(cube)
```

## Arguments

- cube:

  3-D numeric array with dimensions `[nx, ny, nlambda]`.

## Value

List with numeric vectors `bkg` and `mad`, both of length `nlambda`.
