# Spatial PDAC analysis

This folder contains the spatial transcriptomic analyses used to place the selected KCN genes back into tissue context.

The workflow combines broad dataset summaries, ecotype maps, KCN-centered views, and communication-oriented spatial follow-up analyses.

## Main input

- `inputs/PDAC_Updated_ST.rds`

## Main subfolders

- `overview/`
  Dataset-level summaries such as section counts and tissue composition.
- `spatial_maps/`
  Ecotype maps, cell-type maps, and marker projections.
- `kcn/`
  KCN-centered descriptive analyses across ecotypes, cell types, and niche contexts.
- `communication/`
  Hotspot, local correlation, pathway overlay, and MISTy analyses.
- `sparkx/`
  Additional spatially variable gene validation.

## Report links

- `F5A`
  - `overview/`
- `F5B`
  - `spatial_maps/ecotypes/`
- `F5C`
  - `communication/hotspot/all_kcns_primary/`
- `F6A`
  - `communication/hotspot/kcnma1_tumor_epithelial_overlay_primary/`
  - `communication/hotspot/kcnn4_tumor_epithelial_overlay_primary/`
- `F6B`
  - `communication/hotspot/kcnma1_pathway_overlay_primary/`

## Recommended entry points

- `overview/run_overview.R`
- `spatial_maps/ecotypes/run_ecotype_maps.R`
- `communication/hotspot/README.md`
- `communication/misty/README.md`

## Software

- main R packages:
  - `Seurat`
  - `Matrix`
  - `dplyr`
  - `tidyr`
  - `tibble`
  - `ggplot2`
  - `patchwork`
  - `sf`
  - `fgsea`
  - `msigdbr`
  - `escape`
- main Python packages:
  - `scanpy`
  - `anndata`
  - `hotspot`
  - `commot`
  - `numpy`
  - `pandas`
  - `scipy`
  - `matplotlib`
  - `seaborn`

## Practical interpretation

This folder is the main place to answer tissue-level questions such as:
- is a KCN more stromal or more tumor-associated in space?
- is a KCN pattern diffuse or patch-like?
- does a KCN signal follow CAF-related or tumor epithelial programs?
