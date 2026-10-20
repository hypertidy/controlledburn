# Package index

## Overview

- [`controlledburn`](https://hypertidy.github.io/controlledburn/reference/controlledburn-package.md)
  [`controlledburn-package`](https://hypertidy.github.io/controlledburn/reference/controlledburn-package.md)
  : controlledburn: Sparse Geometry Rasterization with Exact Coverage
  Fractions
- [`overview_burn()`](https://hypertidy.github.io/controlledburn/reference/overview_burn.md)
  : Coarse overview of a burn

## Burn

Rasterize polygon, line, and point geometry into the sparse four-table
result.

- [`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
  [`burn_scanline()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
  : Sparse rasterization of polygon, line, and point geometry

## Consume the result

Materialize, crop, and extract from a burn. No dense raster is allocated
until you ask for one.

- [`materialise_chunk()`](https://hypertidy.github.io/controlledburn/reference/materialise_chunk.md)
  [`materialize_chunk()`](https://hypertidy.github.io/controlledburn/reference/materialise_chunk.md)
  : Materialise a controlledburn sparse result into a dense matrix
- [`crop_burn()`](https://hypertidy.github.io/controlledburn/reference/crop_burn.md)
  : Crop a controlledburn result to a target extent
- [`extract_burn()`](https://hypertidy.github.io/controlledburn/reference/extract_burn.md)
  : Extract values from a controlledburn result at point locations

## Methods

Inspect and visualise a burn.

- [`plot(`*`<controlledburn>`*`)`](https://hypertidy.github.io/controlledburn/reference/plot.controlledburn.md)
  : Plot a burn
- [`summary(`*`<controlledburn>`*`)`](https://hypertidy.github.io/controlledburn/reference/summary.controlledburn.md)
  [`print(`*`<summary.controlledburn>`*`)`](https://hypertidy.github.io/controlledburn/reference/summary.controlledburn.md)
  : Summarise a burn without materialising it

## Deprecated

- [`burn_sparse()`](https://hypertidy.github.io/controlledburn/reference/burn_sparse.md)
  : Deprecated: use burn() instead
