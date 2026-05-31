# KCN Mean Positive Expression by Cell Type

This folder contains the mean expression among positive spots for each KCN by `first_type`.

This gives the same type of summary as the ecotype version, but with dominant spot identity, which is sometimes easier biologically.

Main script:
- `run_kcn_mean_positive_expression_by_cell_type.R`

This is useful when the interpretation is easier with dominant spot identity than with ecotypes.

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `dplyr`, `tidyr`, `tibble`, `ggplot2`
- Main input:
  `../../inputs/PDAC_Updated_ST.rds`

