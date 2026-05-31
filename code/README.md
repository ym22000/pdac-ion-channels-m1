# Code

This folder contains the main analysis branches used in the project.

## Main branches

- `moffitt_stromal_exploration_outputs/`
  Stromal exploration and reference visualizations in the Moffitt dataset.
- `moffitt_kcn_deg_3_methods/`
  Conservative KCN overlap without `MAST`.
- `moffitt_kcn_deg_4_methods/`
  Four-method KCN overlap including `MAST`.
- `moffitt_qpsc_vs_mycaf_gsea/`
  Direct pseudobulk Hallmark GSEA between qPSC and myCAF.
- `user_moffitt_validation/`
  External validation in the USER snRNA-seq atlas.
- `bulk/`
  Bulk perturbation analyses, including BKCa silencing in CAFs.
- `spatial_pdac_analysis/`
  Spatial transcriptomic analyses and follow-up interpretation.
- `bulk_survival_analysis/`
  Secondary multicohort survival analyses.
- `project_flowchart/`
  Workflow figure used in the report.

## Fast orientation by report figure

- `F1`
  - conceptual figure assembled in the report
- `F2`
  - `moffitt_stromal_exploration_outputs/`
  - `moffitt_kcn_deg_4_methods/`
  - `moffitt_qpsc_vs_mycaf_gsea/`
- `F3`
  - `user_moffitt_validation/`
- `F4`
  - `bulk/BKCA_shBKCA/`
- `F5`
  - `spatial_pdac_analysis/overview/`
  - `spatial_pdac_analysis/spatial_maps/ecotypes/`
  - `spatial_pdac_analysis/communication/hotspot/all_kcns_primary/`
- `F6`
  - `spatial_pdac_analysis/communication/hotspot/kcnma1_tumor_epithelial_overlay_primary/`
  - `spatial_pdac_analysis/communication/hotspot/kcnn4_tumor_epithelial_overlay_primary/`
  - `spatial_pdac_analysis/communication/hotspot/kcnma1_pathway_overlay_primary/`
- `F7`
  - conceptual figure assembled in the report from the interpretation of `F2` to `F6`

## Languages and software

- Main languages:
  - `R`
  - `Python`
  - shell commands for local execution

- Common R packages:
  - `Seurat`
  - `DESeq2`
  - `MAST`
  - `fgsea`
  - `msigdbr`
  - `ggplot2`
  - `patchwork`
  - `dplyr`
  - `tidyr`
  - `tibble`
  - `Matrix`

- Additional R packages used in specific branches:
  - `escape`
  - `GSVA`
  - `AUCell`
  - `VennDiagram`
  - `survival`
  - `sf`
  - `readxl`
  - `mistyR`

- Common Python packages:
  - `scanpy`
  - `anndata`
  - `hotspot`
  - `commot`
  - `numpy`
  - `pandas`
  - `scipy`
  - `matplotlib`
  - `seaborn`
  - `scikit-learn`

## Suggested reading order

1. Read the root `README.md` for the figure-by-figure guide.
2. Open the branch-specific `README.md`.
3. Rerun the main `run_*.R` or `run_*.py` script from that folder.
4. Inspect the local `figures/`, `tables/`, and `rds/` outputs.
