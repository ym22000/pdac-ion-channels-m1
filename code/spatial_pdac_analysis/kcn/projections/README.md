# KCN Projections

This folder contains the KCN spatial expression maps.

These maps are the direct visual layer of the KCN analysis and make it easy to compare genes across the same tissue sections.

Subfolders:
- `all_sections/`: base KCN maps
- `cc5/`: KCN maps with CC5 contour
- `cc1_cc5/`: KCN maps with combined CC1+CC5 contour
- `cc2_cc3/`: KCN maps with combined CC2+CC3 contour

## Requirements
- `R 4.5.x`
- Main R packages across these projection folders:
  `Seurat`, `dplyr`, `tibble`, `ggplot2`
- The contour-based versions also use:
  `sf`, `RColorBrewer`
- Main input:
  `../../inputs/PDAC_Updated_ST.rds`

