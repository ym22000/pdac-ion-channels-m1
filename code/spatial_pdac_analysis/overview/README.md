# Overview

This folder contains the general overview of the spatial dataset.

It gives a quick way to inspect the object and the main metadata before looking at specific KCN or niche analyses.

Main script:
- `run_overview.R`

Main outputs:
- object summaries
- section counts
- metadata tables
- UMAP figures

Use this folder first if you want a quick check of what is inside the Seurat object.

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `dplyr`, `tidyr`, `tibble`, `ggplot2`, `forcats`, `patchwork`
- Main input:
  `../inputs/PDAC_Updated_ST.rds`

