# KCND3 Pathway Overlay Primary

This folder applies the pathway-activity overlay framework to `KCND3` in
primary PDAC sections, with one custom immune/iCAF-like stromal module.

The goal is to compare:
- the left panel: the `KCND3` Hotspot-style spatial map based on raw exported counts
- the right panel: a per-spot pathway activity score computed with `escape`

The pathway score is a continuous program score calculated from a manual gene
set and then standardized across spots before plotting.

In this local object, the package `escape` expects an assay named `RNA`, while
the spatial data are stored in `Spatial`. The script therefore creates a simple
`RNA` alias from `Spatial` before scoring so the package can run on the
available expression matrix.

Main script:
- `run_kcnd3_pathway_overlay_primary.R`

Pathway/module tested:
- `iCAF_immune_border`

Gene set:
- `C7`
- `C1R`
- `CCL19`
- `CCL21`
- `MFAP4`
- `CLU`
- `DCN`
- `CXCL12`
- `SOD3`
- `OGN`

Why this pathway:
- it captures the immune-border / resident-stroma axis that repeatedly
  co-varies with `KCND3` in the local-correlation analyses
- it provides a focused contrast to the ductal-tumoral and contractile-stromal
  pathway overlays run for `KCNN4` and `KCNMA1`

Outputs:
- one double-panel PDF in `figures/`
- `tables/pathway_gene_set_members.tsv`
- `tables/pathway_gene_set_summary.tsv`
- `tables/pathway_scores_by_spot.tsv`
- `tables/pathway_score_summary_by_section.tsv`
- `tables/kcnd3_hotspot_stats_by_section.tsv`
- `tables/kcnd3_pathway_correlations_by_section.tsv`
- `tables/kcnd3_pathway_correlations_summary.tsv`

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

