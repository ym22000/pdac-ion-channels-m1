# Cell Type Maps

This folder contains the spatial maps colored by the dominant spot identity stored in `first_type`.

These maps give a second view of the tissue that is often easier to interpret in cell-state terms than ecotypes alone.

Main script:
- `run_cell_type_maps.R`

Main output:
- section-wise maps for the selected normal and primary slides

These maps are useful next to the ecotype maps because they give a second label system for the same spots.

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `dplyr`, `tibble`, `ggplot2`, `patchwork`
- Main input:
  `../../inputs/PDAC_Updated_ST.rds`

