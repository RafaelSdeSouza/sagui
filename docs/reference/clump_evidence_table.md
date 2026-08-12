# Summarise evidence for candidate clumps

Builds a compact evidence table for candidate clump footprints. A
candidate is marked as accepted when it is spatially resolved and has
sufficient local compactness and local SED-anomaly evidence.

## Usage

``` r
clump_evidence_table(
  clump_labels,
  score_info,
  min_pixels = 9L,
  min_peak_score = 0.35,
  min_median_score = NULL,
  min_peak_contrast = 0.15,
  min_peak_sed_anomaly = 0.15,
  high_score_threshold = 0.35,
  min_high_score_pixels = 3L,
  probable_peak_score = 0.15,
  probable_peak_contrast = 0.15,
  probable_median_contrast = 0.02,
  probable_peak_sed_anomaly = 0.1,
  clump_relax = 0.5
)
```

## Arguments

- clump_labels:

  Integer clump-footprint matrix.

- score_info:

  Result from
  [`compute_clump_score()`](https://rafaelsdesouza.com.br/sagui/reference/compute_clump_score.md).

- min_pixels:

  Minimum number of pixels in a clump footprint.

- min_peak_score:

  Minimum peak combined clump score. Use `NULL` to skip.

- min_median_score:

  Minimum median combined clump score. Use `NULL` to skip.

- min_peak_contrast:

  Minimum peak local-contrast score. Use `NULL` to skip.

- min_peak_sed_anomaly:

  Minimum peak local SED-anomaly score. Use `NULL` to skip.

- high_score_threshold:

  Threshold used to count high-score pixels.

- min_high_score_pixels:

  Minimum number of pixels above `high_score_threshold`.

- probable_peak_score, probable_peak_contrast, probable_median_contrast,
  probable_peak_sed_anomaly:

  Looser thresholds used to mark candidates as `"probable"` clumps.
  `probable_median_contrast` is a compactness/coherence floor that helps
  reject broad galaxy-profile pieces with coherent SEDs but little local
  brightness contrast.

- clump_relax:

  Single global permissiveness factor for the evidence cuts, bounded to
  `[0, 1]`. `0` is conservative, `0.5` keeps the supplied thresholds
  close to their nominal values, and `1` is permissive. Internally the
  score/contrast thresholds are rescaled smoothly around `0.5`.

## Value

A data frame with clump evidence metrics, `quality`, and `accepted`.
