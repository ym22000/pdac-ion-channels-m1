# Moffitt qPSC vs myCAF GSEA

This folder contains the direct pseudobulk Hallmark comparison between qPSC and myCAF in the stromal Moffitt dataset.

## Goal

Identify broad transcriptional programs that distinguish a more resident-like qPSC state from a more activated, matrix-rich myCAF state.

## Main script

- `run_moffitt_qpsc_vs_mycaf_gsea.R`

## Input

- `../../data/Stroma_Subset2021_fibroblast_focus.rds`

## Main outputs

- `tables/moffitt_qpsc_vs_mycaf_pseudobulk_samples.tsv`
- `tables/moffitt_qpsc_vs_mycaf_deseq2_results.tsv`
- `tables/moffitt_qpsc_vs_mycaf_ranked_stats.tsv`
- `tables/moffitt_qpsc_vs_mycaf_hallmark_fgsea.tsv`
- `tables/moffitt_qpsc_vs_mycaf_hallmark_fgsea_significant.tsv`
- `figures/moffitt_qpsc_vs_mycaf_hallmark_gsea_summary.pdf`
- `figures/moffitt_qpsc_vs_mycaf_hallmark_top_curves.pdf`

## Methods in brief

1. Keep only qPSC and myCAF cells.
2. Build patient-aware pseudobulk samples by `patient x subtype`.
3. Keep paired patients with enough cells in both groups.
4. Run `DESeq2` with `~ patient + subtype`.
5. Rank genes with the Wald statistic and run Hallmark `fgsea`.

## Report link

- `F2D`
