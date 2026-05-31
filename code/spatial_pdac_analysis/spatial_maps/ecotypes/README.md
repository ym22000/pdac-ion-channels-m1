# Ecotype Maps

This folder contains the spatial maps colored by the published ecotype label `cc_ischia_10`.

The ecotype maps are the main reference for reading the tissue with the same label system used in the paper.

Main script:
- `run_ecotype_maps.R`

Main output:
- `figures/normal_and_primary_ecotype_maps.pdf`

This is the main spatial reference used later to interpret KCN patterns.

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `dplyr`, `tibble`, `ggplot2`, `patchwork`
- Main input:
  `../../inputs/PDAC_Updated_ST.rds`

