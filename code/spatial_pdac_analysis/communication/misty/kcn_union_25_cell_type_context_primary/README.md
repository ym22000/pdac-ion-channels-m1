# 25-KCN cell-type-context MISTy in primary PDAC

This folder extends the single-target `KCNMA1` run to the full 25-KCN union.

## Goal
- treat each KCN as a spatial target
- use the `rctd_multi1` cell-type decomposition as predictors
- compare same-spot (`intra`), immediate-neighbor (`juxta_ct`), and broader-context (`para_ct`) effects

## Inputs
- spatial object: `../../../../inputs/PDAC_Updated_ST.rds`
- KCN panel: `../../../../../mycaf_icaf_correlations/kcn_union_25.tsv`
- target assay: `SCT`, layer `data`
- predictor assay: `rctd_multi1`, layer `data`
- sections kept: primary PDAC only (`Origin == "Pancreas"`)

## Main scripts
- `run_kcn_union_25_cell_type_misty_primary.R`
- `make_kcn_union_25_bipartite_network_per_kcn_pdfs.R`
- `make_kcn_union_25_bipartite_network_per_kcn_visual_pdfs.R`

## How to rerun
From this folder:

```bash
Rscript run_kcn_union_25_cell_type_misty_primary.R
Rscript make_kcn_union_25_bipartite_network_per_kcn_pdfs.R
Rscript make_kcn_union_25_bipartite_network_per_kcn_visual_pdfs.R
```

Required R packages:
- `Seurat`
- `dplyr`
- `tibble`
- `tidyr`
- `ggplot2`
- `mistyR`

## Main outputs
- `tables/kcn_union_25_misty_contributions_summary.tsv`
- `tables/kcn_union_25_misty_importances_summary.tsv`
- `tables/kcn_union_25_misty_top_predictors_by_view.tsv`
- `tables/kcn_union_25_misty_network_edges.tsv`
- `tables/kcn_union_25_misty_network_edges_top3.tsv`
- `tables/kcn_union_25_misty_schematic_summary.tsv`
- `figures/kcn_union_25_misty_view_fraction_heatmap.pdf`
- `figures/kcn_union_25_misty_top3_network_heatmap.pdf`
- `figures/by_kcn/`
- `figures/by_kcn_visual/`

## Output structure
- `tables/`: raw summaries used for interpretation and plotting
- `results/`: one MISTy run per `KCN x section`
- `figures/by_kcn/`: simple one-PDF-per-KCN bipartite summaries
- `figures/by_kcn_visual/`: cleaner presentation-style one-PDF-per-KCN summaries

## Quick reading
- `intra` means that the KCN is mainly explained by the composition of the same spot
- `juxta_ct` means that the immediate neighborhood matters most
- `para_ct` means that a broader tissue field contributes most

In practice:
- `KCNMA1`, `KCNMB1`, `KCNMB4`, `KCNE4`, `KCND2`, and `KCNJ8` fall in a mainly fibroblast-associated context
- `KCNN4` is the clearest tumor-epithelial case
- `KCNK3` is the strongest endothelial / TAM-context example

## Notes
- `KCNA7` is absent from the `SCT` assay and is not modeled here
- predictors with zero variance are removed section by section before fitting
- view comparison is based on mean absolute view fractions, because raw view coefficients can be signed
- MISTy importances are not directional

