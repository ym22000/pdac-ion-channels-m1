# Bulk RNA-seq Analyses

This folder groups the local bulk perturbation analyses used in the project.

Subfolders:
- `BKCA_shBKCA/`: BKCa / KCNMA1 perturbation dataset
- `SIGMAR1_shSIGMAR1/`: SigmaR1 perturbation dataset

These analyses were used to see how ion-channel perturbation is linked to CAF-like programs and pathway changes in bulk RNA-seq.

## Requirements
- `R 4.5.x`
- `BKCA_shBKCA/` uses a lighter R stack:
  `ggplot2`, `dplyr`, `fgsea`, `msigdbr`, `pheatmap`
- `SIGMAR1_shSIGMAR1/` uses a fuller DESeq2 workflow:
  `DESeq2`, `SummarizedExperiment`, `AnnotationDbi`, `org.Hs.eg.db`, `apeglm`, `ggplot2`, `dplyr`, `fgsea`, `msigdbr`, `pheatmap`, `viridis`, `viridisLite`
- Each dataset folder also expects its own local `inputs/` and `resources/` files to already be present.

