# KCNN4 Pathway Overlay Primary

This folder applies the same pathway-activity overlay logic as the `KCNMA1`
comparison, but focuses on `KCNN4` and three pro-tumoral epithelial programs in
primary PDAC sections.

The goal is to compare:
- the left panel: the `KCNN4` Hotspot-style spatial map based on raw exported counts
- the right panel: a per-spot pathway activity score computed with `escape`

The pathway score is not a single-gene readout. It is a continuous program
score calculated from a manual gene set and then standardized across spots
before plotting.

In this local object, the package `escape` expects an assay named `RNA`, while
the spatial data are stored in `Spatial`. The script therefore creates a simple
`RNA` alias from `Spatial` before scoring so the package can run on the
available expression matrix.

Main script:
- `run_kcnn4_pathway_overlay_primary.R`

Pathways/modules tested:
- `ductal_tumoral_epithelial`
- `tumor_stress_secretory`
- `Tumor_proliferation_rate`

Why these pathways:
- `ductal_tumoral_epithelial` tests whether `KCNN4` follows a broad ductal
  epithelial tumor program
- `tumor_stress_secretory` tests whether `KCNN4` follows a more stressed,
  secretory, injury-like tumor layer
- `Tumor_proliferation_rate` tests whether `KCNN4` also follows the cell-cycle
  and proliferative tumor program used in the reference PDAC paper

Outputs:
- three double-panel PDFs in `figures/`
- `tables/pathway_gene_set_members.tsv`
- `tables/pathway_gene_set_summary.tsv`
- `tables/pathway_scores_by_spot.tsv`
- `tables/pathway_score_summary_by_section.tsv`
- `tables/kcnn4_hotspot_stats_by_section.tsv`
- `tables/kcnn4_pathway_correlations_by_section.tsv`
- `tables/kcnn4_pathway_correlations_summary.tsv`

Method summary:
- Seurat object: `inputs/PDAC_Updated_ST.rds`
- sections shown: primary PDAC only
- scoring method: `escape::runEscape(..., method = "ssGSEA", groups = 1000)`
- pathway scaling: z-score per pathway across scored spots
- left panel colormap: viridis
- right panel colormap: magma

Requirements:
- `R 4.5.x`
- R packages:
  `Seurat`, `Matrix`, `dplyr`, `tibble`, `ggplot2`, `patchwork`, `escape`, `viridisLite`

