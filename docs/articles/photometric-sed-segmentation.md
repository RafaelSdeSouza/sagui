# SED Segmentation

Sagui clusters one photometric SED per supported pixel. The support mask
and clustering transform are chosen before measuring the regional
fluxes.

## Exact Ward

``` r
seg <- segment_regions(
  input = cube,
  Ncomp = 20,
  use_starlet_mask = TRUE,
  cluster_pretransform = "none"
)
```

Exact Ward clustering stores all pairwise distances. Use it only when
the supported pixel count fits available memory.

## Sparse Ward

``` r
seg <- segment_regions_large(
  input = cube,
  Ncomp = 20,
  use_starlet_mask = TRUE,
  cluster_pretransform = "none",
  knn_k = 40,
  auto_k = FALSE
)
```

The sparse backend restricts candidate merges to a nearest-neighbour
graph. Its result depends on `knn_k`, `auto_k`, the support mask, and
the transform.

## Choose transforms

Use `cluster_pretransform`, not the deprecated `pretransform` alias. The
default `"none"` preserves the input representation. Alternatives such
as `"copula_gaussian"` change the clustering geometry; compare them
against `"none"`.

Region labels are categorical. Their integer values do not represent a
physical ordering.
