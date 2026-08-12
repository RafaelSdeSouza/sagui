# Sagui website content map

This map records the Phase B migration from the former Quarto site to one
pkgdown documentation system. No scientific claim is declared obsolete.

| Source content | Destination | Treatment | Parity notes |
|---|---|---|---|
| `index.qmd` purpose and hero | `README.md` / homepage | Retained and tightened | Semantic H1, compact local logo, verified paper links, and real manuscript mosaic. |
| Three feature blocks | Homepage workflow | Recast | Preparation, support, segmentation, and regional measurement remain separate. |
| Installation and first segmentation | Homepage and Getting Started | Replaced by role | Fixed-seed nine-band galaxy simulation replaces the local-file-only first run. |
| `index.qmd` manuscript gallery | Paper Examples | Moved | Real assets retained with file-level provenance status and categorical-label caveats. |
| Exact and sparse runs | SED Segmentation | Moved | Uses `cluster_pretransform`; deprecated `pretransform` is not taught. |
| Design choices | Preparation, Support Masks, SED Segmentation | Expanded | PSF matching, support, transforms, and flux measurement are distinguished. |
| `examples.qmd` common setup and data products | Getting Started and Paper Examples | Redistributed | Local user paths remain explicitly non-executable templates. |
| Sagui-10 workflow | Paper Examples | Retained | Intended settings remain; unavailable source cube and tables are disclosed. |
| Target-S/N selection | Choose Region Count and Paper Examples | Retained | Producer script and machine-readable tables are linked. |
| Bagpipes result | Paper Examples and Python hand-off | Retained | Optional downstream modelling is not presented as part of segmentation. |
| `python.qmd` CSV/FITS examples | Python and Astropy | Condensed | Region IDs, zero-background convention, and WCS caveat remain explicit. |
| Function list | Generated Reference index | Replaced | Every documented public export is grouped by scientific role. |
| Citation text | `inst/CITATION`, sidebar, Paper nav | Canonicalised | Fourteen named authors, collaboration note, journal record, DOI, and arXiv ID. |
| Quarto styling | `pkgdown/extra.css` | Reimplemented | Astronomical blue/coral identity, responsive layout, focus states, and native themes. |
| Remote COIN brand dependency | Text affiliation in footer and About | Removed | The local site no longer depends on a remote logo; COIN affiliation is retained as text. |

## Redirect parity

- `/examples.html` → `/articles/paper-examples-reproduction.html`
- `/python.html` → `/articles/python-astropy-handoff.html`
- Legacy paper-example fragments retained: `#sagui-10-main-example`,
  `#sagui-11-faint-bridge`, and `#target-sn-selection`.
- Homepage anchors retained: `#installation`, `#quick-start`, `#core-api`, and
  `#about`.

The former Quarto files may be retired only after pkgdown build, redirect,
fragment, link, and rendered-content checks pass.
