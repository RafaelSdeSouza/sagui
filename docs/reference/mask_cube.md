# Mask a 3D cube by a spatial label map

Sets pixels outside selected regions to zero or NA across all bands.

## Usage

``` r
mask_cube(cube, labels, mode = c("zero", "na"))
```

## Arguments

- cube:

  3-D numeric array with dimensions `[H, W, B]`.

- labels:

  2-D matrix `[H, W]`; values greater than zero are kept.

- mode:

  "zero" or "na" — how to fill masked pixels.

## Value

Masked cube with same dimensions as input.
