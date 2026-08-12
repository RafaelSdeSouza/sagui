# Export Regional SEDs

Measure regional SEDs from the calibrated flux cube, using the
segmentation map to group pixels.

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

Save the band axis, PSF target, support mask, and segmentation
parameters with the tables.

## Keep labels intact

Do not renumber labels between the map and the table unless the mapping
is also saved. The label numbers only identify regions.
