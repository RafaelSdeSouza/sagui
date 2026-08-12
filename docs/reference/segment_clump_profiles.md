# Split clump footprints into internal profile regions

Split clump footprints into internal profile regions

## Usage

``` r
segment_clump_profiles(clump_labels, score, clump_levels = 3L)
```

## Arguments

- clump_labels:

  Integer clump-footprint matrix.

- score:

  Score map used to rank pixels inside each clump.

- clump_levels:

  Number of internal profile regions per clump.

## Value

A list with `profile_labels`, `clump_labels`, and `summary`.
