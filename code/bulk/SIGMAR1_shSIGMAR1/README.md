# SigmaR1 / shSigmaR1 bulk RNA-seq analysis

## Summary
This folder contains a cleaned and reproducible downstream analysis of the SigmaR1 perturbation bulk RNA-seq dataset, centered on essential QC, CAF-oriented transcriptional programs, and the combined CAF heatmap without apCAF.

## Requirements
- `R 4.5.x`
- Main R packages:
  `DESeq2`, `SummarizedExperiment`, `AnnotationDbi`, `org.Hs.eg.db`, `apeglm`, `ggplot2`, `dplyr`, `fgsea`, `msigdbr`, `pheatmap`, `viridis`, `viridisLite`
- Local files expected in this folder:
  `inputs/dds.RData`, `inputs/vsd.RData`, and the marker files stored in `resources/`

## Key points
- collapsed genes: 14729
- samples: 6 (3 CTRL, 3 shSigmaR1)
- DEG `padj < 0.05`: 1546
- top Wilcoxon signature: `myCAF_top30` (p=0.0809, FDR=0.4043)
- strongest positive `SIGMAR1` correlation: `EZR` (rho=0.943, FDR=0.2714)
- strongest negative `SIGMAR1` correlation: `CXADR` (rho=-1, FDR=0.1583)
- top Hallmark pathway: `HALLMARK_G2M_CHECKPOINT` (NES=-2.047, FDR=1.8e-10)

## Main figures
- `figures/PCA_samples.pdf`
- `figures/PCA_SIGMAR1_module.pdf`
- `figures/Volcano_plot.pdf`
- `figures/CAF_marker_violin_panels.pdf`
- `figures/CAF_gene_contribution_barplot.pdf`
- `figures/CAF_signature_heatmap.pdf`
- `figures/CAF_signature_score_panels.pdf`
- `figures/GSEA_Hallmark_bubbleplot.pdf`

## Main tables
- `tables/expression_summary.csv`
- `tables/DESeq2_results_gene_level.csv`
- `tables/CAF_marker_statistics.csv`
- `tables/CAF_signature_statistics.csv`
- `tables/CAF_signature_reference.csv`
- `tables/CAF_gene_contribution.csv`
- `tables/SIGMAR1_marker_correlations.csv`
- `tables/Hallmark_fgsea_results.csv`
- `rds/sigmar1_shsigmar1_bulk_analysis_bundle.rds`

