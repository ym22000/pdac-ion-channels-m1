# KCN Fisher Niche Enrichment In Primary PDAC

This folder contains the Fisher tests used to quantify whether a KCN is enriched in CC1+CC5 or CC2+CC3.

The goal here is to move beyond visual impression and quantify whether a KCN is overrepresented in one of the main niche groups.

Main script:
- `run_kcn_niche_fisher_tests.R`

Main outputs:
- odds ratios
- p-values
- FDR tables
- summary tables by gene and niche

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `dplyr`, `tidyr`, `tibble`
- Main input:
  `../../inputs/PDAC_Updated_ST.rds`

