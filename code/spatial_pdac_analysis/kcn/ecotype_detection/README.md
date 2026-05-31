# KCN Detection by Ecotype

This folder contains the percent-positive summaries for the 25-gene KCN list by ecotype.

This is one of the simplest ways to see where a KCN is detected most often across the published spatial ecotypes.

Main script:
- `run_kcn_detection_by_ecotype.R`

Main outputs:
- one barplot PDF per KCN
- summary tables by ecotype

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `dplyr`, `tidyr`, `tibble`, `ggplot2`
- Main input:
  `../../inputs/PDAC_Updated_ST.rds`

