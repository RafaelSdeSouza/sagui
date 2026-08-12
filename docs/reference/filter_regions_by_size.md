# Filter regions by pixel count

Filter regions by pixel count

## Usage

``` r
filter_regions_by_size(labels, min_size = 30L, max_size = Inf)
```

## Arguments

- labels:

  Matrix of integer labels.

- min_size:

  Minimum region size to keep.

- max_size:

  Maximum region size to keep.

## Value

Relabeled matrix containing only regions within the requested size
range.
