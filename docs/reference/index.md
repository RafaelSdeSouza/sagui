# Package index

## Core segmentation

Exact and sparse photometric SED segmentation, plus region-count
selection.

- [`segment_regions()`](https://rafaelsdesouza.com.br/sagui/reference/segment_regions.md)
  : Segment photometric regions from a spectral cube
- [`segment_regions_large()`](https://rafaelsdesouza.com.br/sagui/reference/segment_regions_large.md)
  : Segment photometric regions from large spectral cubes
- [`choose_ncomp_by_snr()`](https://rafaelsdesouza.com.br/sagui/reference/choose_ncomp_by_snr.md)
  : Choose the Number of SAGUI Regions from a Target S/N

## Photometry preparation

Cube conversion, display transforms, masking, and spectrum transforms.

- [`cube_to_matrix()`](https://rafaelsdesouza.com.br/sagui/reference/cube_to_matrix.md)
  : Flatten a spectral cube into a pixel-by-band matrix
- [`collapse_white_light()`](https://rafaelsdesouza.com.br/sagui/reference/collapse_white_light.md)
  : Collapse a spectral cube into a white-light image (robust, weighted)
- [`mask_cube()`](https://rafaelsdesouza.com.br/sagui/reference/mask_cube.md)
  : Mask a 3D cube by a spatial label map
- [`pretransform_cube()`](https://rafaelsdesouza.com.br/sagui/reference/pretransform_cube.md)
  : Apply a spectral pretransform to a full cube
- [`pretransform_spectra()`](https://rafaelsdesouza.com.br/sagui/reference/pretransform_spectra.md)
  : Apply a spectral pretransform before segmentation
- [`median_scale()`](https://rafaelsdesouza.com.br/sagui/reference/median_scale.md)
  : Center a spectrum by its median
- [`asinh_stretch()`](https://rafaelsdesouza.com.br/sagui/reference/asinh_stretch.md)
  : Asinh-stretch an image or numeric array

## Support construction

Multiscale support builders and mask refinement.

- [`build_starlet_mask()`](https://rafaelsdesouza.com.br/sagui/reference/build_starlet_mask.md)
  : Build a starlet-based segmentation mask from a spectral cube
- [`starlet_mask()`](https://rafaelsdesouza.com.br/sagui/reference/starlet_mask.md)
  : Starlet (à trous) decomposition
- [`starlet_reconstruct()`](https://rafaelsdesouza.com.br/sagui/reference/starlet_reconstruct.md)
  : Reconstruct from a subset of starlet scales
- [`b3_kernel()`](https://rafaelsdesouza.com.br/sagui/reference/b3_kernel.md)
  : B3-spline low-pass kernel (1D)
- [`upsample_kernel()`](https://rafaelsdesouza.com.br/sagui/reference/upsample_kernel.md)
  : Insert zeros between kernel taps
- [`build_adaptive_support()`](https://rafaelsdesouza.com.br/sagui/reference/build_adaptive_support.md)
  : Build an adaptive multi-band support mask
- [`build_contourlet_mask()`](https://rafaelsdesouza.com.br/sagui/reference/build_contourlet_mask.md)
  : Build a contourlet-like directional support mask from a spectral
  cube
- [`clean_segmentation_mask()`](https://rafaelsdesouza.com.br/sagui/reference/clean_segmentation_mask.md)
  : Clean a logical segmentation support mask
- [`threshold_starlet_regions()`](https://rafaelsdesouza.com.br/sagui/reference/threshold_starlet_regions.md)
  : Reconstruct and threshold starlet scales to produce a region mask

## Regional SED products

Flux-preserving regional photometry and background summaries.

- [`RegionPhotometry()`](https://rafaelsdesouza.com.br/sagui/reference/RegionPhotometry.md)
  : Region-level photometry on hyperspectral cubes (no background
  subtraction)
- [`extract_region_sed()`](https://rafaelsdesouza.com.br/sagui/reference/extract_region_sed.md)
  : Extract integrated SEDs for segmented regions
- [`summarize_region_sed()`](https://rafaelsdesouza.com.br/sagui/reference/summarize_region_sed.md)
  : Summarize integrated SEDs for segmented regions
- [`estimate_bkg_mad_per_lambda()`](https://rafaelsdesouza.com.br/sagui/reference/estimate_bkg_mad_per_lambda.md)
  : Robust per-wavelength background and scatter

## Maps and post-processing

Categorical maps, region properties, relabelling, and smoothing.

- [`plot_region_map()`](https://rafaelsdesouza.com.br/sagui/reference/plot_region_map.md)
  : Plot a segmented region map with dissolved Voronoi-style polygons
- [`plot_region_property_map()`](https://rafaelsdesouza.com.br/sagui/reference/plot_region_property_map.md)
  : Plot a region property map in Voronoi-style polygons
- [`plot_region_property_raster()`](https://rafaelsdesouza.com.br/sagui/reference/plot_region_property_raster.md)
  : Plot a region property map on a pixel grid
- [`reindex_regions()`](https://rafaelsdesouza.com.br/sagui/reference/reindex_regions.md)
  : Reindex positive region labels to a compact sequence
- [`fill_region_holes()`](https://rafaelsdesouza.com.br/sagui/reference/fill_region_holes.md)
  : Fill zero-valued holes inside each labeled region
- [`filter_regions_by_size()`](https://rafaelsdesouza.com.br/sagui/reference/filter_regions_by_size.md)
  : Filter regions by pixel count
- [`reject_background_regions()`](https://rafaelsdesouza.com.br/sagui/reference/reject_background_regions.md)
  : Reject background-like regions after segmentation
- [`smooth_region_field_laplacian()`](https://rafaelsdesouza.com.br/sagui/reference/smooth_region_field_laplacian.md)
  : Smooth region values over an adjacency graph

## Clump candidates

Optional candidate detection and diagnostic helpers.

- [`detect_sed_clumps()`](https://rafaelsdesouza.com.br/sagui/reference/detect_sed_clumps.md)
  : Detect compact SED-clump seeds
- [`segment_clumps()`](https://rafaelsdesouza.com.br/sagui/reference/segment_clumps.md)
  : Segment compact SED clumps and the diffuse galaxy body
- [`segment_clump_profiles()`](https://rafaelsdesouza.com.br/sagui/reference/segment_clump_profiles.md)
  : Split clump footprints into internal profile regions
- [`compute_clump_score()`](https://rafaelsdesouza.com.br/sagui/reference/compute_clump_score.md)
  : Compute a compact SED-clump score
- [`clump_evidence_table()`](https://rafaelsdesouza.com.br/sagui/reference/clump_evidence_table.md)
  : Summarise evidence for candidate clumps
- [`grow_clump_footprints()`](https://rafaelsdesouza.com.br/sagui/reference/grow_clump_footprints.md)
  : Grow clump seeds into compact footprints
- [`plot_clump_diagnostics()`](https://rafaelsdesouza.com.br/sagui/reference/plot_clump_diagnostics.md)
  : Plot clump-segmentation diagnostics

## Package

- [`sagui`](https://rafaelsdesouza.com.br/sagui/reference/sagui-package.md)
  [`sagui-package`](https://rafaelsdesouza.com.br/sagui/reference/sagui-package.md)
  : sagui
