# Data

This folder stores the main curated RDS objects reused across the project.

Files:
- `Stroma_Subset2021.rds`: stromal-focused scRNA-seq object
- `Stroma_Subset2021_fibroblast_focus.rds`: fibroblast-focused version of the stromal object
- `PDAC_Updated_ST.rds`: spatial transcriptomics Seurat object

These files are treated as analysis inputs and are not rebuilt here.

## Source studies

- `Stroma_Subset2021.rds`
- `Stroma_Subset2021_fibroblast_focus.rds`
  - derived from the stromal single-cell PDAC resource reported by Oh et al., *Nature Communications*.
- `PDAC_Updated_ST.rds`
  - derived from the spatial PDAC resource reported by Khaliq et al., *Nature Genetics*.

These local objects are curated working copies prepared for the project workflow. The corresponding public studies should be cited when the datasets are reused, and the original source papers remain the primary reference for data generation and study design.

## Requirements
- No extra software is needed to browse this folder.
- To use these files in the project scripts, you mainly need `R 4.5.x` with `Seurat`.
- `PDAC_Updated_ST.rds` is used by the spatial workflow.
- `Stroma_Subset2021.rds` and `Stroma_Subset2021_fibroblast_focus.rds` are used by the stromal single-cell workflows.

