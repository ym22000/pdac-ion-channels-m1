# Niches

This folder groups the contour-based niche plots.

These plots make the niche boundaries explicit, which is useful when a KCN seems to sit inside or next to a specific territory.

Subfolders:
- `cc5/`: CC5 contours
- `cc1_cc5_combined/`: combined CC1+CC5 contours
- `cc2_cc3_combined/`: combined CC2+CC3 contours

These plots are useful when we want to compare KCN expression with a niche boundary instead of only with spot colors.

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `dplyr`, `tibble`, `ggplot2`, `sf`, `patchwork`
- Main input:
  `../inputs/PDAC_Updated_ST.rds`
- These scripts also rely on valid spatial coordinates stored in the Seurat image slots.

