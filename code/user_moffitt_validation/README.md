# USER validation against Moffitt stromal states

This folder contains the USER validation branch of the project. It is used to place stromal KCN candidates back into a broader PDAC microenvironment that includes tumor, fibroblast, immune, endothelial, and other compartments.

## What this folder contains

- a fibroblast-focused validation workflow against Moffitt stromal states
- a global USER UMAP and patient-composition summary
- cell-level and pseudobulk KCN mapping across USER cell subsets
- a direct patient-aware `Tumor vs Fibroblast` pseudobulk comparison
- a transcriptome-wide `Tumor vs Fibroblast` pseudobulk volcano with KCN highlights

## Main scripts

- `run_user_moffitt_validation.R`
- `kcn25_mast_cell_subsets/run_user_kcn25_mast_cell_subsets.R`
- `kcn25_pseudobulk_cell_subsets/run_user_kcn25_pseudobulk_cell_subsets.R`
- `kcn25_pseudobulk_tumor_vs_fibroblast/run_user_kcn25_pseudobulk_tumor_vs_fibroblast.R`
- `all_genes_pseudobulk_tumor_vs_fibroblast/run_user_all_genes_pseudobulk_tumor_vs_fibroblast.R`

## Main report-related outputs

- `figures/Figure_01_USER_global_umap_cell_subsets.pdf`
- `figures/Figure_06_USER_patient_cell_subset_composition.pdf`
- `all_genes_pseudobulk_tumor_vs_fibroblast/figures/USER_all_genes_pseudobulk_tumor_vs_fibroblast_volcano.pdf`

## Core tables

- `tables/USER_full_metadata.tsv`
- `tables/USER_fibro_cluster_annotation_summary.tsv`
- `tables/USER_fibro_cleaning_rules.tsv`
- `all_genes_pseudobulk_tumor_vs_fibroblast/tables/USER_all_genes_pseudobulk_tumor_vs_fibroblast.tsv`
- `all_genes_pseudobulk_tumor_vs_fibroblast/tables/USER_kcn25_within_all_genes_tumor_vs_fibroblast.tsv`
- `kcn25_pseudobulk_tumor_vs_fibroblast/tables/USER_kcn25_pseudobulk_tumor_vs_fibroblast.tsv`

## Inputs

- `inputs/gene_sorted-naivedata_scp.mtx`
- `inputs/naivedata_scp.genes.csv`
- `inputs/naivedata_scp.barcodes.csv`
- `inputs/combinenaivedata-reprocessed-clean-detailed-annotations.tsv`
- `inputs/combinenaivedata-reprocessed-clean-detailed-UMAP.tsv`
- `resources/Moffitt_SuppData1_CellType2_Markers.xlsx`

## Methods in brief

1. Rebuild the global USER metadata table and reuse the published UMAP coordinates.
2. Restrict one branch to fibroblasts and project these cells onto Moffitt stromal programs.
3. Recluster fibroblasts and build a cleaned stromal view.
4. Run cell-level `MAST` for sensitive KCN localization across cell subsets.
5. Run patient-aware pseudobulk `DESeq2` for more conservative cell-subset and `Tumor vs Fibroblast` contrasts.

## Practical interpretation

- `MAST` is the most sensitive screening layer.
- pseudobulk `DESeq2` is the patient-aware validation layer.
- the direct `Tumor vs Fibroblast` contrast is the clearest USER readout for comparing stromal and tumor tendencies of the selected KCN genes.
