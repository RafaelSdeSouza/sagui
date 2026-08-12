# Grow clump seeds into compact footprints

Grow clump seeds into compact footprints

## Usage

``` r
grow_clump_footprints(
  clump_seeds,
  support,
  score = NULL,
  brightness = NULL,
  grow_radius = 3L,
  max_radius = 7L,
  min_pixels = 9L,
  footprint_quantile = 0.1,
  core_drop_frac = 0.2,
  sigma_threshold = 3,
  footprint_mode = c("connected", "radial"),
  connectivity = 8L
)
```

## Arguments

- clump_seeds:

  Data frame with `clump_id`, `row`, and `col`.

- support:

  Logical foreground support mask.

- score:

  Optional score map used to trim each local footprint.

- brightness:

  Optional local-brightness/contrast map used by connected footprints to
  stop growth once pixels fall back to the local background.

- grow_radius:

  Initial dilation radius around each seed.

- max_radius:

  Maximum radius tried when a clump is too small.

- min_pixels:

  Minimum footprint size per clump.

- footprint_quantile:

  Optional local score quantile used to keep only the stronger part of
  each dilated footprint. Use `NULL` to keep all pixels.

- core_drop_frac:

  For `footprint_mode = "connected"`, keep connected pixels above this
  fraction of the seed-core brightness contrast relative to the local
  background.

- sigma_threshold:

  For `footprint_mode = "connected"`, keep connected pixels above this
  robust local-background significance threshold. Use `NULL` to disable.

- footprint_mode:

  Footprint-growth mode. `"connected"` uses the seed only as an anchor
  and returns the connected high-score component, allowing elongated or
  irregular clump shapes. `"radial"` keeps compact radius-based growth.

- connectivity:

  Pixel connectivity used by `footprint_mode = "connected"`.

## Value

A list with `clump_labels` and `footprint_summary`.
