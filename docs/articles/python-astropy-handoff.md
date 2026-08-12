# Python and Astropy

Sagui returns a categorical map and regional photometry tables. CSV is
portable; FITS preserves the spatial WCS with the map.

## Export from R

``` r
write.csv(regional$flux_long, "regional_seds_long.csv", row.names = FALSE)

label_image <- seg$cluster_map
label_image[is.na(label_image)] <- 0L
FITSio::writeFITSim(label_image, file = "region_labels.fits")
```

Here, zero means outside support or unassigned; positive integers are
categorical region identifiers. Copy the original FITS header to
preserve the WCS.

## Read with Python

``` python
import pandas as pd
from astropy.io import fits

sed = pd.read_csv("regional_seds_long.csv")
labels = fits.getdata("region_labels.fits")

assert (sed["region"] > 0).all()
assert labels.ndim == 2
```

Use the `region` column, rather than row order, to join SED-fitting
results back to the map.

## Fit downstream

Sagui does not fit physical SED models. Pass its regional fluxes and
uncertainties, together with the filters and redshift, to the fitting
code.
