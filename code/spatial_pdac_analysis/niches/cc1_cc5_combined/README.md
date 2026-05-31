# Combined CC1 and CC5 Contours

This folder contains the contour maps where CC1 and CC5 are treated as one merged niche.

CC1 and CC5 often looked close enough spatially that it was useful to test them as one broader niche.

Main script:
- `run_cc1_cc5_combined_contours_on_grey_spots.R`

Main output:
- grey spot maps with one contour around the combined CC1+CC5 islands

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `dplyr`, `tibble`, `ggplot2`, `sf`
- Main input:
  `../../inputs/PDAC_Updated_ST.rds`

