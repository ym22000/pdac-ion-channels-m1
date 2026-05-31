# KCN Wilcoxon Ecotype Expression Contrasts In Primary PDAC

This folder contains a section-wise alternative to the pooled Fisher niche test.

For each primary PDAC section, each KCN is compared between one ecotype group
and the rest of the section. The main contrasts are:
- `CC1+CC5` vs rest
- `CC2+CC3` vs rest

The goal is to keep the analysis closer to the actual structure of the dataset:
- one section at a time
- expression values kept continuous
- one FDR correction per section and per contrast

Main script:
- `run_kcn_ecotype_expression_contrasts_primary.R`
- `make_kcn_ecotype_summary_pdf.py`
- `make_kcn_ecotype_by_section_pdf.py`

Main outputs:
- `tables/section_contrast_availability.tsv`
- `tables/kcn_ecotype_vs_rest_by_section.tsv`
- `tables/kcn_ecotype_vs_rest_summary.tsv`
- `tables/kcn_preferred_ecotype_contrast_primary.tsv`
- `reports/kcn_ecotype_summary_primary.pdf`
- `reports/kcn_ecotype_by_section_primary.pdf`

Metrics reported:
- numbers of inside/outside spots
- percentages of positive spots
- mean expression and mean-positive expression
- `log2FC` based on mean expression
- `AUC` computed directly from the Wilcoxon rank-sum / Mann-Whitney statistic
- Wilcoxon p-value and FDR

This block was added because the pooled Fisher test is easy to compute but too
coarse for sparse Visium signals. The section-wise contrast keeps more
information and lets the project assess whether a KCN effect is repeated across
multiple primary sections instead of being driven by a single pooled table.

## Requirements
- `R 4.5.x`
- Main R packages:
  `Seurat`, `dplyr`, `tidyr`, `tibble`
- Main input:
  `../../inputs/PDAC_Updated_ST.rds`

