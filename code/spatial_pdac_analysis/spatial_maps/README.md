# Spatial Maps

This folder contains the basic spatial reference maps used in the project.

These maps provide the visual tissue context needed before interpreting KCN patterns or niche-level statistics.

Subfolders:
- `ecotypes/`: maps colored by `cc_ischia_10`
- `cell_types/`: maps colored by `first_type`
- `marker_projections/`: projections of a few reference markers

These maps are the starting point for reading the tissue before looking at KCN-specific analyses.

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `dplyr`, `tibble`, `ggplot2`, `patchwork`
- Main input:
  `../inputs/PDAC_Updated_ST.rds`

