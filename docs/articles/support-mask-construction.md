# Support Masks

A support mask decides which pixels are eligible for segmentation. It is
not a segmentation map.

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

The support depends on the collapse function, scale range, threshold,
and cleanup settings.

## Adaptive support

[`build_adaptive_support()`](https://rafaelsdesouza.com.br/sagui/reference/build_adaptive_support.md)
combines evidence across bands. It is useful when no single collapse
represents the source adequately.

``` r
support <- build_adaptive_support(cube, pretransform = "none")
```

For directional structures,
[`build_contourlet_mask()`](https://rafaelsdesouza.com.br/sagui/reference/build_contourlet_mask.md)
provides an alternative to starlets. Compare both masks before
segmentation.

## Background values

Use `FALSE` in a logical support. The segmentation functions return
unassigned background as `NA` in `cluster_map`.
