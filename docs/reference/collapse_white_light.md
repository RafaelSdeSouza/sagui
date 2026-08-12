# Collapse a spectral cube into a white-light image (robust, weighted)

Steps:

1.  subtract spatial median per-wavelength,

2.  clip very negative outliers at -k \* MAD,

3.  (optional) inverse-variance weights 1/MAD^2,

4.  weighted sum and reshape to a 2-D image,

5.  force non-negative for visualization safety.

## Usage

``` r
collapse_white_light(cube, kclip = 2, use_weights = TRUE)
```

## Arguments

- cube:

  3-D numeric array with dimensions `[nx, ny, nlambda]`.

- kclip:

  numeric; negative clip in units of MAD (default 2)

- use_weights:

  logical; if TRUE, use 1/MAD^2 weights (default TRUE)

## Value

2-D matrix with dimensions `[nx, ny]`.
