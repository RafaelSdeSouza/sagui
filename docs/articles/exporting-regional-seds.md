# Export Regional SEDs

Measure regional SEDs on the flux cube intended for science, using the
categorical label map learned from the segmentation representation.

## Extract

``` r
regional <- extract_region_sed(
  cube = science_cube,
  labels = seg$cluster_map,
  var_cube = science_variance,
  band_values = bands
)
```

`flux_long` contains one row per region and band. `flux_wide` provides
one row per region for tools that expect tabular SEDs. `flux_err`, when
a variance cube is supplied, follows the propagated regional variance.

## Write tables

``` r
write.csv(regional$flux_long, "regional_seds_long.csv", row.names = FALSE)
write.csv(regional$flux_wide, "regional_seds_wide.csv", row.names = FALSE)
```

Store a provenance record beside the tables: input cube, variance cube,
band axis, PSF target, support settings, segmentation settings, random
seed where applicable, package version, and output paths.

## Keep labels intact

Do not renumber labels between the map and the table unless the mapping
is also written. Region identifiers are join keys, not measurements.
