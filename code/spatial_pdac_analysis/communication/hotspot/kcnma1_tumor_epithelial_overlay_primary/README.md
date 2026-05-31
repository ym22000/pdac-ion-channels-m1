# KCNMA1 tumor epithelial signature overlay in primary PDAC

This folder scores a small tumor epithelial signature and compares it to the spatial pattern of `KCNMA1`.

## Signature
- `EPCAM`
- `KRT19` (`CK19`)
- `MUC1`
- `KRT7`
- `KRT8`
- `KRT18`

`PanCK` is not a gene symbol, so it is approximated here with epithelial keratins.

## Main script
- `run_kcnma1_tumor_epithelial_overlay_primary.R`

## Main outputs
- `tables/pathway_gene_set_summary.tsv`
- `tables/kcnma1_pathway_correlations_by_section.tsv`
- `tables/kcnma1_pathway_correlations_summary.tsv`
- `figures/tumor_epithelial_kcnma1_hotspot_vs_pathway_primary.pdf`

