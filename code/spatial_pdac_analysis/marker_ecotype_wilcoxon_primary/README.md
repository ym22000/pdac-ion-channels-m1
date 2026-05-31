# Marker Wilcoxon Ecotype Expression Contrasts In Primary PDAC

This folder contains the same section-wise ecotype-vs-rest contrast framework
used for KCN genes, but applied to a small set of stromal / myCAF-like markers:

- `ACTA2`
- `CCN2`
- `POSTN`
- `TAGLN`
- `COL12A1`

For each primary PDAC section, each marker is compared between:
- `CC1+CC5` vs rest
- `CC2+CC3` vs rest

Main script:
- `run_marker_ecotype_expression_contrasts_primary.R`
- `make_marker_ecotype_by_section_pdf.py`

Main outputs:
- `tables/marker_gene_availability.tsv`
- `tables/section_contrast_availability.tsv`
- `tables/marker_ecotype_vs_rest_by_section.tsv`
- `tables/marker_ecotype_vs_rest_summary.tsv`
- `reports/marker_ecotype_by_section_primary.pdf`

Metrics reported:
- numbers of inside/outside spots
- percentages of positive spots
- mean expression and mean-positive expression
- `log2FC` based on mean expression
- `AUC` computed directly from the Wilcoxon rank-sum / Mann-Whitney statistic
- Wilcoxon p-value and FDR

