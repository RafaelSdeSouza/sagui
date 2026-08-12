# Plot clump-segmentation diagnostics

Plot clump-segmentation diagnostics

## Usage

``` r
plot_clump_diagnostics(
  x,
  mode = c("overlay", "score", "profiles", "merged"),
  alpha = 0.72
)
```

## Arguments

- x:

  Result from
  [`segment_clumps()`](https://rafaelsdesouza.com.br/sagui/reference/segment_clumps.md).

- mode:

  Plot mode: `"score"` shows the clump score, `"overlay"` overlays
  clump-profile regions on the score, `"profiles"` shows clump profiles
  only, and `"merged"` shows the final merged segmentation.

- alpha:

  Clump-profile opacity in overlay mode.

## Value

A `ggplot2` object.
