# Reject background-like regions after segmentation

Applies a simple, explicit region-level background classifier to a SAGUI
segmentation. This is designed to be used after an intentionally loose
support and optional oversegmentation: let the clustering isolate
candidate structures, then remove regions with weak multi-band support.

## Usage

``` r
reject_background_regions(
  segmentation,
  evidence = NULL,
  evidence_count = NULL,
  min_median_evidence = NULL,
  min_median_persistence = 1,
  low_evidence_quantile = 0.2,
  reject_edge_low_evidence = TRUE,
  reject_large_low_evidence = TRUE,
  large_area_quantile = 0.75
)
```

## Arguments

- segmentation:

  A result returned by
  [`segment_regions()`](https://rafaelsdesouza.com.br/sagui/reference/segment_regions.md)
  or
  [`segment_regions_large()`](https://rafaelsdesouza.com.br/sagui/reference/segment_regions_large.md).

- evidence:

  Optional 2-D evidence map. If `NULL`, the function tries to use
  `segmentation$support$details$evidence`.

- evidence_count:

  Optional 2-D band-persistence map. If `NULL`, the function tries to
  use `segmentation$support$details$evidence_count`.

- min_median_evidence:

  Minimum median evidence required for a region. If `NULL`, uses
  `low_evidence_quantile` of region median evidence values.

- min_median_persistence:

  Minimum median band-persistence count required for a region.

- low_evidence_quantile:

  Quantile used to define low-evidence regions when
  `min_median_evidence = NULL`.

- reject_edge_low_evidence:

  Logical; reject edge-connected regions with low median evidence.

- reject_large_low_evidence:

  Logical; reject large regions with low median evidence.

- large_area_quantile:

  Quantile of region areas used to identify large regions when
  `reject_large_low_evidence = TRUE`.

## Value

The input segmentation with `cluster_map` updated and with
`background_diagnostics` and `background_regions` appended.
