# KCN Mean Positive Expression by Ecotype

This folder contains the mean expression among positive spots for each KCN by ecotype.

Percent-positive alone can miss intensity differences, so this folder adds the strength of expression where the gene is present.

Main script:
- `run_kcn_mean_positive_expression_by_ecotype.R`

This is complementary to the percent-positive summaries because it focuses on signal strength where the gene is detected.

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `dplyr`, `tidyr`, `tibble`, `ggplot2`
- Main input:
  `../../inputs/PDAC_Updated_ST.rds`

