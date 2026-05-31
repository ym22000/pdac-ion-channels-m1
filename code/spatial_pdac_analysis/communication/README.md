# Communication

This folder contains the spatial communication-style analyses used in the project.

These analyses were added to test whether KCN-rich territories also have specific spatial signaling or gene-coordination patterns.
They now include both the main primary-PDAC analyses and a normal-pancreas comparison block for Hotspot and territory-level GSEA.
They also include two focused follow-up branches:
- a `SIGMAR1` branch, kept separate from the KCN panel
- a `KCNMA1` BK-subunit branch, used to compare `KCNMA1+KCNMB1` and `KCNMA1+KCNMB4`
- a `KCNMA1` pathway-overlay branch, used to compare the Hotspot map of `KCNMA1`
  against continuous stromal pathway activity scores
- a `KCNN4` pathway-overlay branch, used to compare the Hotspot map of `KCNN4`
  against continuous pro-tumoral epithelial pathway activity scores
- a `KCND3` pathway-overlay branch, used to compare the Hotspot map of `KCND3`
  against a continuous immune/iCAF-like stromal pathway activity score

Subfolders:
- `commot/`: ligand-receptor communication inference
- `hotspot/`: spatial autocorrelation, local correlations, and territory GSEA
  across primary PDAC sections, plus mirrored autocorrelation, local-correlation,
  and territory-GSEA runs in normal pancreas sections, along with focused
  `SIGMAR1`, BK-subunit, and pathway-overlay context analyses

## Requirements
- This part of the project needs both `R 4.5.x` and `Python 3.10+`.
- The R side is mainly used for exporting data from the Seurat object and for the territory GSEA step.
- The Python side is used for `COMMOT`, `Hotspot`, and local correlation runs.
- Main R packages:
  `Seurat`, `Matrix`, `dplyr`, `tibble`, `ggplot2`, `fgsea`, `msigdbr`, `patchwork`, `escape`
- Main Python packages:
  `anndata`, `hotspot`, `commot`, `scanpy`, `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`, `scikit-learn`

