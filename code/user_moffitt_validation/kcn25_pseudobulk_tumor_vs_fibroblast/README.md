# USER KCN pseudobulk Tumor vs Fibroblast

This folder contains the patient-aware pseudobulk `Tumor vs Fibroblast` comparison restricted to the 25 KCN union genes in the USER dataset.

## Design

- keep only `Tumor` and `Fibroblast`
- keep only patients with both compartments
- aggregate raw counts by `patient x cell_subset`
- run `DESeq2` with `~ patient + group`

Positive `log2FC` means `Tumor up`. Negative `log2FC` means `Fibroblast up`.

## Main script

- `run_user_kcn25_pseudobulk_tumor_vs_fibroblast.R`

## Main outputs

- `tables/USER_kcn25_pseudobulk_tumor_vs_fibroblast.tsv`
- `tables/USER_kcn25_pseudobulk_tumor_vs_fibroblast_significant.tsv`
- `tables/USER_kcn25_tumor_vs_fibroblast_samples.tsv`
