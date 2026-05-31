# SIGMAR1 Context In Primary PDAC

This folder extends the focused `SIGMAR1` Hotspot run with the same three
context analyses already used for the KCN panel and for the `KCNMA1` subunit
modules:

- Hotspot local correlations
- `first_type` KNN neighborhood around `SIGMAR1-high` spots
- territory Hallmark GSEA for `SIGMAR1-high` spots

The goal is to place `SIGMAR1` on the same interpretative scale as the KCN
spatial analyses:

- does `SIGMAR1` live with a recognizable local gene module?
- what cell-type neighborhood surrounds `SIGMAR1-high` territories?
- what biological program characterizes those territories?

Scripts:
- `run_sigmar1_local_correlations_primary.py`
- `run_sigmar1_first_type_neighborhood.py`
- `run_sigmar1_territory_gsea_primary.R`

Outputs:
- `local_correlations/tables/`
- `local_correlations/figures/sigmar1_local_correlations_primary_sections.pdf`
- `first_type_neighborhood/tables/`
- `first_type_neighborhood/figures/sigmar1_first_type_knn_neighborhood.pdf`
- `territory_gsea/tables/`
- `territory_gsea/figures/sigmar1_territory_hallmark_gsea.pdf`

How to rerun:
1. run the parent Hotspot screen in `../`
2. run `run_sigmar1_local_correlations_primary.py`
3. run `run_sigmar1_first_type_neighborhood.py`
4. run `run_sigmar1_territory_gsea_primary.R`

The `local_correlations/tables/` folder also keeps `_checkpoint_*` files so a
long rerun can still be inspected before the script finishes writing the final
summary tables.

Method settings:
- Hotspot model: `bernoulli`
- neighbors: `30` for Hotspot local correlations
- KNN neighborhood: `15` neighbors
- `SIGMAR1-high`: top `10%` of positive spots, minimum `20` high spots and `30`
  positive spots
- territory GSEA assay: `SCT/data`

