# Combined CC2 and CC3 Contours

This folder contains the contour maps where CC2 and CC3 are treated as one merged niche.

CC2 and CC3 often appeared as a related tumor-associated territory, so this folder keeps the merged contour version.

Main script:
- `run_cc2_cc3_combined_contours_on_grey_spots.R`

Main output:
- grey spot maps with one contour around the combined CC2+CC3 islands

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `dplyr`, `tibble`, `ggplot2`, `sf`
- Main input:
  `../../inputs/PDAC_Updated_ST.rds`

