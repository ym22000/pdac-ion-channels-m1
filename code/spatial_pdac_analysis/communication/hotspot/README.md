# Hotspot

This folder contains the Hotspot-based analyses.

Hotspot was used to test whether KCN signals form real spatial patches and to identify genes or programs linked to those patches.
The folder now contains both the main primary-PDAC run and a normal-pancreas comparison run built with the same framework.

Subfolders:
- `all_kcns_primary/`: main Hotspot results for the KCN set
- `all_kcns_normal_pancreas/`: matching Hotspot results for the three normal pancreas sections
- `kcnd3_pathway_overlay_primary/`: pathway-activity overlay for KCND3 in primary PDAC
- `kcnn4_pathway_overlay_primary/`: pathway-activity overlays for KCNN4 in primary PDAC
- `kcnma1_subunits/`: focused Hotspot run for KCNMA1 combined with KCNMB1 or KCNMB4
- `kcnma1_pathway_overlay_primary/`: pathway-activity overlays for KCNMA1 in primary PDAC
- `marker_panel_primary/`: Hotspot run for the spatial reference marker panel plus POSTN
- `sigmar1_primary/`: focused Hotspot run for SIGMAR1 across primary PDAC sections
- `marker_local_correlations_primary/`: tables-only local-correlation run for C7/C1R/CCL19/CCL21/MFAP4/CLU/DCN

This is where the spatial autocorrelation, local correlations, and KCN-high territory GSEA are grouped together.
The normal-pancreas branch was added to compare resident KCN spatial modules against PDAC-specific modules and to keep the healthy-tissue territory GSEA separate from the tumor-focused runs.
The `kcnd3_pathway_overlay_primary/` branch was added to compare the KCND3 Hotspot map against a focused immune/iCAF-like stromal pathway score, matching the resident-border interpretation suggested by the local-correlation analyses.
The `kcnn4_pathway_overlay_primary/` branch was added to compare the KCNN4 Hotspot map against continuous pro-tumoral epithelial pathway scores, using the same spot-level pathway framework as the KCNMA1 overlay branch.
The `kcnma1_subunits/` branch was added to test whether BK-related spatial modules become clearer when the alpha subunit KCNMA1 is combined with a beta subunit in the same pseudo-gene signal.
The `kcnma1_pathway_overlay_primary/` branch was added to compare the KCNMA1 Hotspot map directly against continuous stromal pathway scores computed spot by spot with `escape`, following the logic used in the reference PDAC paper, and now also includes an immune-oriented contrast control.
The `sigmar1_primary/` branch was added to place `SIGMAR1` on the same spatial scale as the KCN analyses without mixing it into the KCN panel itself.

Recommended order when rerunning this folder:
1. `all_kcns_primary/` for the reference KCN panel
2. `all_kcns_normal_pancreas/` for the healthy-tissue comparison
3. `kcnma1_subunits/` for BK-oriented follow-up analyses
4. `kcnma1_pathway_overlay_primary/` for KCNMA1-versus-pathway overlay maps
5. `kcnn4_pathway_overlay_primary/` for KCNN4-versus-tumor-program overlay maps
6. `kcnd3_pathway_overlay_primary/` for KCND3-versus-immune/iCAF overlay maps
7. `sigmar1_primary/` for the focused `SIGMAR1` comparison

## Requirements
- `Python 3.10+` for spatial autocorrelation and local correlations
- `R 4.5.x` for the territory GSEA step
- Main Python packages:
  `anndata`, `hotspot`, `numpy`, `pandas`, `scipy`, `matplotlib`
- Main R packages:
  `Seurat`, `Matrix`, `dplyr`, `tibble`, `ggplot2`, `fgsea`, `msigdbr`, `patchwork`, `escape`
- These analyses reuse the section exports generated from the spatial Seurat object.

