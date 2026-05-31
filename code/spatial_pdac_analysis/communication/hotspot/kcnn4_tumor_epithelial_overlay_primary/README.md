# KCNN4 tumor epithelial signature overlay in primary PDAC

This folder scores a small tumor epithelial signature and compares it to the spatial pattern of `KCNN4`.

## Signature
- `EPCAM`
- `KRT19` (`CK19`)
- `MUC1`
- `KRT7`
- `KRT8`
- `KRT18`

`PanCK` is not a gene symbol, so it is approximated here with epithelial keratins.

## Main script
- `run_kcnn4_tumor_epithelial_overlay_primary.R`

## Main outputs
- `tables/pathway_gene_set_summary.tsv`
- `tables/kcnn4_pathway_correlations_by_section.tsv`
- `tables/kcnn4_pathway_correlations_summary.tsv`
- `figures/tumor_epithelial_kcnn4_hotspot_vs_pathway_primary.pdf`

