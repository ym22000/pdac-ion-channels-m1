# KCNMA1 cell-type-context MISTy in primary PDAC

This folder contains a focused `mistyR` run for one target: `KCNMA1`.

## Goal
- model `KCNMA1` as a spatial target
- use the `rctd_multi1` cell-type decomposition as predictors
- separate same-spot (`intra`), immediate-neighbor (`juxta_ct`), and broader-context (`para_ct`) effects

## Inputs
- spatial object: `../../../../inputs/PDAC_Updated_ST.rds`
- target assay: `SCT`, layer `data`
- predictor assay: `rctd_multi1`, layer `data`
- sections kept: primary PDAC only (`Origin == "Pancreas"`)

## Main script
- `run_kcnma1_cell_type_misty_primary.R`

## How to rerun
From this folder:

```bash
Rscript run_kcnma1_cell_type_misty_primary.R
```

Required R packages:
- `Seurat`
- `dplyr`
- `tibble`
- `tidyr`
- `ggplot2`
- `mistyR`

## Main outputs
- `tables/kcnma1_misty_section_info.tsv`
- `tables/kcnma1_misty_section_status.tsv`
- `tables/kcnma1_misty_contributions_summary.tsv`
- `tables/kcnma1_misty_importances_summary.tsv`
- `tables/kcnma1_misty_top_predictors_by_view.tsv`
- `figures/kcnma1_misty_view_contributions.pdf`
- `figures/kcnma1_misty_predictor_importance_heatmap.pdf`

## Quick reading
- `intra` means that `KCNMA1` is mainly explained by the cell-type composition of the same spot
- `juxta_ct` means that the immediate neighborhood adds explanatory power
- `para_ct` means that a broader tissue context contributes

In this run, `KCNMA1` is mainly explained by a local fibroblast-rich context, with additional contribution from the surrounding tumor/TAM environment.

## Notes
- predictors with zero variance are removed section by section before fitting
- MISTy importances are not directional: they indicate which predictors matter most, not whether the association is positive or negative

