# Segment compact SED clumps and the diffuse galaxy body

Experimental clump-aware segmentation mode. The method builds a loose
foreground support, detects compact SED-anomalous seeds (or accepts
supplied clump seeds), grows them into clump footprints, splits each
clump into internal profile levels, and runs
[`segment_regions_large()`](https://rafaelsdesouza.com.br/sagui/reference/segment_regions_large.md)
on the remaining diffuse body. The final `cluster_map` merges
clump-profile regions and the diffuse-body segmentation.

## Usage

``` r
segment_clumps(
  input,
  Ncomp_body = 80,
  bands = NULL,
  support = NULL,
  support_info = NULL,
  support_args = list(),
  clump_search_mask = NULL,
  clump_reject_mask = NULL,
  clump_reject_fraction = 0.5,
  clump_seeds = NULL,
  n_clumps = NULL,
  max_clumps = 50L,
  score_quantile = 0.995,
  min_seed_distance = 3,
  clump_levels = 3L,
  grow_radius = 3L,
  max_radius = 7L,
  min_pixels = 9L,
  footprint_quantile = 0.1,
  core_drop_frac = 0.2,
  sigma_threshold = 3,
  footprint_mode = c("connected", "radial"),
  connectivity = 8L,
  contrast_weight = 0.5,
  sed_weight = 0.5,
  small_sigma = 0.7,
  large_sigma = 4,
  apply_quality_filter = FALSE,
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
  clump_relax = 0.5,
  body_segment = TRUE,
  knn_k = 40,
  auto_k = FALSE,
  max_k = NULL,
  feature_scale = c("none", "robust_col"),
  spatial_weight = 0.1,
  cluster_pretransform = "none",
  verbose = TRUE
)
```

## Arguments

- input:

  3-D array or FITS-like list with `imDat`. This should usually be the
  PSF-matched cube used for segmentation.

- Ncomp_body:

  Number of diffuse-body regions.

- bands:

  Bands used for support, clump scoring, and body segmentation.

- support:

  Optional precomputed logical support mask.

- support_info:

  Optional support object, e.g. from
  [`build_adaptive_support()`](https://rafaelsdesouza.com.br/sagui/reference/build_adaptive_support.md).
  If `NULL`, adaptive support is built.

- support_args:

  Named list passed to
  [`build_adaptive_support()`](https://rafaelsdesouza.com.br/sagui/reference/build_adaptive_support.md)
  when `support_info = NULL` and `support = NULL`.

- clump_search_mask:

  Optional logical mask with the same spatial dimensions as the cube.
  When supplied, automatic clump seeds are searched only inside this
  mask, while the diffuse body can still be segmented over the full
  foreground support. Use this to target context-specific compact
  structures, e.g. off-disc jellyfish knots, in-disc HII regions, or
  halo compact-source candidates.

- clump_reject_mask:

  Optional logical mask with the same spatial dimensions as the cube.
  Candidate clump footprints whose overlap with this mask exceeds
  `clump_reject_fraction` are removed before profile splitting. This is
  a generic context filter, not a science label.

- clump_reject_fraction:

  Maximum allowed footprint overlap with `clump_reject_mask`.

- clump_seeds:

  Optional data frame with `clump_id`, `row`, and `col`. If `NULL`,
  seeds are detected automatically.

- n_clumps:

  Optional exact number of automatic clump seeds.

- max_clumps:

  Maximum number of automatic clump seeds.

- score_quantile:

  Quantile threshold for automatic clump seeds.

- min_seed_distance:

  Minimum seed separation in pixels.

- clump_levels:

  Number of internal profile regions per clump.

- grow_radius:

  Initial clump-footprint growth radius.

- max_radius:

  Maximum clump-footprint growth radius.

- min_pixels:

  Minimum pixels per grown clump footprint.

- footprint_quantile:

  Local score quantile used to trim grown footprints.

- core_drop_frac:

  For connected footprints, keep connected pixels above this fraction of
  the seed-core brightness contrast relative to the local background.

- sigma_threshold:

  For connected footprints, keep connected pixels above this robust
  local-background significance threshold. Use `NULL` to disable.

- footprint_mode:

  Footprint-growth mode. `"connected"` uses the seed only as an anchor
  and returns the connected high-score component, allowing elongated or
  irregular clump shapes. `"radial"` keeps the older compact
  radius-based growth.

- connectivity:

  Pixel connectivity used by `footprint_mode = "connected"`.

- contrast_weight:

  Weight of compact local contrast in the clump score.

- sed_weight:

  Weight of local SED anomaly in the clump score.

- small_sigma:

  Compact smoothing scale for the clump score.

- large_sigma:

  Background smoothing scale for the clump score.

- apply_quality_filter:

  Logical; if `TRUE`, discard candidate clumps that fail the evidence
  cuts before profile splitting.

- min_peak_score, min_median_score, min_peak_contrast,
  min_peak_sed_anomaly:

  Evidence cuts passed to
  [`clump_evidence_table()`](https://rafaelsdesouza.com.br/sagui/reference/clump_evidence_table.md).

- high_score_threshold, min_high_score_pixels:

  High-score-pixel criterion passed to
  [`clump_evidence_table()`](https://rafaelsdesouza.com.br/sagui/reference/clump_evidence_table.md).

- probable_peak_score, probable_peak_contrast, probable_median_contrast,
  probable_peak_sed_anomaly:

  Looser evidence thresholds used for `"probable"` clumps.

- clump_relax:

  Single global permissiveness factor for the clump evidence cuts,
  bounded to `[0, 1]`. `0` is conservative, `0.5` keeps the nominal
  thresholds, and `1` accepts fainter/weaker compact candidates. This
  does not use SED-fitting results.

- body_segment:

  Logical; segment the remaining diffuse body.

- knn_k, auto_k, max_k, feature_scale, spatial_weight:

  Sparse-Ward arguments passed to
  [`segment_regions_large()`](https://rafaelsdesouza.com.br/sagui/reference/segment_regions_large.md)
  for the diffuse body.

- cluster_pretransform:

  Spectral pretransform for diffuse-body segmentation.

- verbose:

  Logical; print progress messages.

## Value

A Sagui-like list with merged `cluster_map`, clump diagnostics, and body
segmentation products.
