# USER KCN pseudobulk cell subsets

This folder contains the patient-aware pseudobulk one-vs-rest analysis of the 25 KCN union genes across the main USER cell subsets.

## Goal

Estimate where each KCN is enriched across the global USER cell landscape with a more conservative design than cell-level `MAST`.

## Design

- keep only the 25 KCN union genes
- aggregate raw counts by `patient x cell_subset`
- for each target cell subset, compare:
  - `target`
  - `rest`
- run `DESeq2` with `~ patient + group`

## Main script

- `run_user_kcn25_pseudobulk_cell_subsets.R`

## Main outputs

- `tables/USER_kcn25_pseudobulk_by_cell_subset.tsv`
- `tables/USER_kcn25_pseudobulk_by_cell_subset_significant.tsv`
- `tables/USER_kcn25_pseudobulk_cell_subset_counts.tsv`
- `tables/USER_kcn25_pseudobulk_samples.tsv`
- `tables/USER_kcn25_pseudobulk_comparison_info.tsv`
