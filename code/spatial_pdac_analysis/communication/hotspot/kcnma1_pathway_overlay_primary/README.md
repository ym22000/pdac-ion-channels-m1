# KCNMA1 Pathway Overlay Primary

This folder reproduces the pathway-activity logic used in the reference PDAC
paper, but focuses on `KCNMA1` in primary PDAC sections.

The goal is to compare:
- the left panel: `KCNMA1` Hotspot-style spatial map based on raw exported counts
- the right panel: a per-spot pathway activity score computed with `escape`

The pathway score is not a single-gene expression value. It is a continuous
program score calculated from a manual gene set and then standardized across
primary PDAC spots before plotting.
In this local object, the package `escape` expects an assay named `RNA`, while
the spatial data are stored in `Spatial`. The script therefore creates a simple
`RNA` alias from `Spatial` before scoring so the package can run on the
available expression matrix.

Main script:
- `run_kcnma1_pathway_overlay_primary.R`

Pathways/modules tested:
- `Cancer_associated_fibroblasts`
- `Matrix`
- `Matrix_remodeling`
- `myCAF_contractile`
- `NK_cells` (immune-oriented contrast control)

Outputs:
- four double-panel PDFs in `figures/`
- `tables/pathway_gene_set_members.tsv`
- `tables/pathway_gene_set_summary.tsv`
- `tables/pathway_scores_by_spot.tsv`
- `tables/pathway_score_summary_by_section.tsv`
- `tables/kcnma1_hotspot_stats_by_section.tsv`
- `tables/kcnma1_pathway_correlations_by_section.tsv`
- `tables/kcnma1_pathway_correlations_summary.tsv`

Method summary:
- Seurat object: `inputs/PDAC_Updated_ST.rds`
- spots used for scoring: primary PDAC spots only
- scoring method: `escape::runEscape(..., method = "ssGSEA", groups = 1000)`
- pathway scaling: z-score per pathway across primary PDAC spots
- left panel colormap: viridis
- right panel colormap: magma

Why these four pathways:
- `Cancer_associated_fibroblasts` tests broad fibroblast overlap
- `Matrix` tests ECM-rich territories
- `Matrix_remodeling` tests active stromal remodeling
- `myCAF_contractile` tests the most specific contractile/myofibroblastic
  hypothesis for `KCNMA1`
- `NK_cells` acts as an immune contrast control to check that `KCNMA1`
  is not simply following any structured pathway signal

Requirements:
- `R 4.5.x`
- R packages:
  `Seurat`, `Matrix`, `dplyr`, `tibble`, `ggplot2`, `patchwork`, `escape`

