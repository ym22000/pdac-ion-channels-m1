# USER all-gene Tumor vs Fibroblast pseudobulk

## Summary
This folder contains the patient-aware pseudobulk `DESeq2` analysis for the direct `Tumor vs Fibroblast` comparison across all USER genes.

The design matches the focused KCN-only run:
- same USER input files
- same patient pairing strategy
- same model: `~ patient + group`

The main addition is a transcriptome-wide volcano plot with:
- `Tumor up` genes in `#7400b8`
- `Fibroblast up` genes in `#72efdd`
- significant union KCN genes highlighted in black

## Main outputs
- `tables/USER_all_genes_pseudobulk_tumor_vs_fibroblast.tsv`
- `tables/USER_kcn25_within_all_genes_tumor_vs_fibroblast.tsv`
- `tables/USER_all_genes_pseudobulk_tumor_vs_fibroblast_info.tsv`
- `figures/USER_all_genes_pseudobulk_tumor_vs_fibroblast_volcano.pdf`
- `rds/USER_all_genes_pseudobulk_tumor_vs_fibroblast_results.rds`

## How to run

```r
Rscript run_user_all_genes_pseudobulk_tumor_vs_fibroblast.R
```

## Notes
- Only patients with both `Tumor` and `Fibroblast` cells are retained.
- Genes with fewer than `10` total counts across pseudobulk samples are removed before `DESeq2`.
- The volcano labels the top significant union KCN genes.

