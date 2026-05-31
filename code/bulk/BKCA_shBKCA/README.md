# BKCa / shBKCa bulk RNA-seq analysis

This folder contains the bulk RNA-seq analysis of BKCa silencing in CAFs. It is the main perturbation branch used in the report to test whether the BKCa axis is associated with coherent transcriptomic changes.

## Design

- samples: `8`
  - `4 CTRL`
  - `4 shBKCa`
- gene-level matrix after collapsing: `17,889` genes
- main contrast: `shBKCa vs CTRL`

## Main script

- `run_bkca_shbkca_bulk_rnaseq_analysis.R`

## Main outputs

- `figures/Volcano_plot.pdf`
- `figures/GSEA_Hallmark_bubbleplot.pdf`
- `figures/CAF_signature_heatmap_myCAF_iCAF_ifCAF_landscape.pdf`
- `tables/DESeq2_results_gene_level.csv`
- `tables/CAF_signature_statistics.csv`
- `tables/CAF_signature_reference.csv`
- `tables/Hallmark_fgsea_results.csv`
- `rds/bkca_shbkca_bulk_analysis_bundle.rds`

## Report links

- `F4B`
  - `figures/Volcano_plot.pdf`
- `F4C`
  - `figures/GSEA_Hallmark_bubbleplot.pdf`
- `F4D`
  - `figures/CAF_signature_heatmap_myCAF_iCAF_ifCAF_landscape.pdf`

## Methods in brief

1. Import the local bulk RNA-seq counts generated in the laboratory.
2. Run `DESeq2` on raw counts for the `shBKCa vs CTRL` contrast.
3. Rank genes and run Hallmark `fgsea`.
4. Score selected CAF signatures and summarize their condition-level shifts.

## Interpretation in one line

This branch is used to test whether BKCa silencing in CAFs is linked to changes in CAF-related programs, especially myCAF-like and ifCAF-like directions.
