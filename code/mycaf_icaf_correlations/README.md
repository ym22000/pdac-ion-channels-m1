# KCN - myCAF/iCAF Correlations

This folder contains the pseudobulk correlation analysis between the KCN gene list and human myCAF / iCAF markers in the stromal single-cell dataset.

Main script:
- `run_kcn_mycaf_icaf_correlations.R`

Main inputs:
- `kcn_union_25.tsv`
- `genes_top30_myCAF_human.txt`
- `genes_top30_iCAF_human.txt`

Main outputs:
- `results/tables/`
- `results/figures/`

This analysis was used to test whether the KCN candidates align more with myCAF-like or iCAF-like programs.

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `Matrix`, `ggplot2`, `dplyr`, `tidyr`
- Main inputs:
  `kcn_union_25.tsv`, `genes_top30_myCAF_human.txt`, `genes_top30_iCAF_human.txt`
- This analysis also expects access to the stromal Seurat object from the main project data.

