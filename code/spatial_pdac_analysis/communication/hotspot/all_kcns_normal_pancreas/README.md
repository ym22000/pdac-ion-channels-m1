# Hotspot for All KCNs in Normal Pancreas Sections

This folder contains the Hotspot analyses run on the three normal pancreas
sections kept in the spatial project.

The goal here is simple: keep the same spatial autocorrelation and local
correlation framework used for primary PDAC sections, but apply it to a
normal-only context. This makes it easier to compare which KCN patterns look
resident or homeostatic in normal tissue and which ones become more tumor-like
in PDAC.

Main scripts:
- `run_export_normal_pancreas_sections.R`
- `run_hotspot_all_kcns_normal_pancreas.py`

Subfolders:
- `local_correlations/`: local gene-gene correlations around each KCN
- `kcn_territory_gsea/`: Hallmark GSEA on normal-pancreas KCN-high territories
- `exports_normal_pancreas/`: exported section matrices and coordinates reused
  by the Python scripts

Main outputs:
- one Hotspot PDF per available KCN
- summary tables by gene and section
- one local-correlation PDF per available KCN
- one normal-pancreas territory-GSEA PDF per KCN, including placeholder pages
  when too few sections pass the territory filters

## Requirements
- `R 4.5.x`
- `Python 3.10+`
- Main R packages:
  `Seurat`, `dplyr`, `tibble`, `Matrix`, `ggplot2`, `fgsea`, `msigdbr`
- Main Python packages:
  `anndata`, `hotspot`, `numpy`, `pandas`, `scipy`, `matplotlib`
- Main reused input:
  `inputs/PDAC_Updated_ST.rds`

