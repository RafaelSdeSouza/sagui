# Paper Examples

This page indexes the real examples already present in the repository.
It does not imply that every upstream survey cutout can be redistributed
or reproduced without its original data access and calibration context.

## Sagui-10

![Sagui-10 photometric SED segmentation. Colours are categorical region
labels.](../reference/figures/sagui10_n20.png)

Sagui-10 photometric SED segmentation. Colours are categorical region
labels.

The site-ready image is `images/examples/sagui10_n20.png`. The regional
SED fit shown below is retained as a verified repository asset; consult
the provenance table before treating it as reproducible from public
inputs alone.

![Regional SED fit retained with the Sagui-10
example.](../reference/figures/sagui10_sedfit_region_1.png)

Regional SED fit retained with the Sagui-10 example.

## Faint bridge

![Sagui-11 low-surface-brightness composite and copula-pretransformed
segmentation.](../reference/figures/sagui11_lsb_composite.png)

Sagui-11 low-surface-brightness composite and copula-pretransformed
segmentation.

![Sagui-11 low-surface-brightness composite and copula-pretransformed
segmentation.](../reference/figures/sagui11_lsb_copula.png)

Sagui-11 low-surface-brightness composite and copula-pretransformed
segmentation.

The copula transform is an analysis choice, not a neutral display
transform. Compare it with an untransformed run before interpreting the
bridge.

## Target S/N

![Sagui-10 candidate selected under a minimum regional S/N threshold of
10.](../reference/figures/sagui10_target_snr_10.png)

Sagui-10 candidate selected under a minimum regional S/N threshold of
10.

The producer is `scripts/make_sagui10_target_snr_panel.R`;
machine-readable candidate and selection tables live under
`images/examples/target_snr/sagui10/`. The threshold must be read
together with the variance model and candidate grid.

## Bagpipes hand-off

![Bagpipes model photometry for Sagui-8 region 8, generated from a Sagui
regional SED.](../reference/figures/sagui8_region8_bagpipes_sedfit.png)

Bagpipes model photometry for Sagui-8 region 8, generated from a Sagui
regional SED.

The producer is `scripts/make_sagui8_bagpipes_sedfit_example.py`; its
generated model photometry is in `outputs/website_sedfit/`. Bagpipes
remains an optional downstream interface rather than part of Sagui’s
segmentation algorithm.
