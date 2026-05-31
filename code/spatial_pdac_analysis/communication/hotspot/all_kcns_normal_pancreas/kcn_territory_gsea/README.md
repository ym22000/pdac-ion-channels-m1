# KCN-High Territory GSEA in Normal Pancreas

This folder contains the Hallmark GSEA run on KCN-high territories in the
three normal pancreas sections.

This analysis was added to describe the broad biological programs linked to the
tissue territory where a KCN is strongest in healthy tissue. It mirrors the
primary-PDAC territory GSEA, but keeps the normal sections separate so the
interpretation stays biologically clean.

Main scripts:
- `run_all_kcn_territory_gsea_normal_pancreas.R`

Main outputs:
- one Hallmark GSEA PDF per KCN
- territory summary tables
- top positive and negative genes

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `dplyr`, `tibble`, `Matrix`, `ggplot2`, `fgsea`, `msigdbr`
- Main input:
  the spatial Seurat object in `../../../inputs/PDAC_Updated_ST.rds`
- This analysis defines the territory inside each section from the top fraction of positive spots for the target KCN.

Here, the territory is defined inside each section as the top part of the positive spots for the target KCN.

