# Support Masks

A support mask decides which pixels are eligible for segmentation. It is
not a region map and it does not encode physical ordering.

## Starlet support

The default path collapses the cube, reconstructs selected starlet
scales, and thresholds the result.

``` r
support <- build_starlet_mask(
  input = cube,
  starlet_J = 5,
  starlet_scales = 2:5,
  denoise_k = 2.5,
  positive_only = TRUE
)
image(support$mask)
```

Record the collapse function, scale range, threshold, and cleanup
settings. Changing them changes the eligible footprint.

## Adaptive support

[`build_adaptive_support()`](https://rafaelsdesouza.com.br/sagui/reference/build_adaptive_support.md)
combines evidence across bands. It is useful when no single collapse
represents the source adequately.

``` r
support <- build_adaptive_support(cube, pretransform = "none")
```

Directional contourlet support is available but remains a more
experimental choice. Inspect its directional response and compare it
with a neutral support before drawing scientific conclusions.

## Background values

Use `FALSE` in a logical support and `NA` or zero only according to the
downstream function’s documented convention. The segmentation functions
return unassigned background as `NA` in `cluster_map`.
