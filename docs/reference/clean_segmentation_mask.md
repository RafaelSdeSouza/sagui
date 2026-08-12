# Clean a logical segmentation support mask

Optionally applies binary close/open morphology and removes small
connected components from a 2-D support mask.

## Usage

``` r
clean_segmentation_mask(
  mask,
  min_mask_area = 1L,
  close_size = 1L,
  open_size = 1L,
  keep_largest = FALSE,
  connectivity = c(8L, 4L)
)
```

## Arguments

- mask:

  Logical or numeric matrix. Positive/`TRUE` finite pixels are kept.

- min_mask_area:

  Minimum connected-component area to retain.

- close_size:

  Size of the optional binary closing brush. Use `1` to skip.

- open_size:

  Size of the optional binary opening brush. Use `1` to skip.

- keep_largest:

  Logical; keep only the largest connected component after area
  filtering.

- connectivity:

  Pixel connectivity for component labeling. Use `8` for diagonal
  connectivity or `4` for edge-only connectivity.

## Value

A logical matrix with the same dimensions as `mask`.
