# KCN Analyses

This folder contains the KCN-focused analyses of the spatial dataset.

This is the core KCN part of the ST project, from simple detection summaries to niche-level tests and projections.

Subfolders:
- `ecotype_detection/`
- `ecotype_mean_positive/`
- `wilcoxon_ecotype_expression_contrasts_primary/`
- `cell_type_mean_positive/`
- `projections/`
- `fisher_niche_enrichment_primary/`
- `kcnma1_examples/`

Together, these folders describe where the KCN genes are detected, how strong the signal is, and how it relates to the main niche groups.

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `dplyr`, `tidyr`, `tibble`, `ggplot2`
- Some KCN projection scripts also use:
  `sf`, `patchwork`, `RColorBrewer`
- Main input:
  `../inputs/PDAC_Updated_ST.rds`

