# MISTy spatial context analyses

This folder contains `mistyR` analyses used to place selected KCN genes in a cell-type spatial context.

## Current subfolders
- `kcnma1_cell_type_context_primary/`: focused MISTy run for `KCNMA1`
- `kcn_union_25_cell_type_context_primary/`: same analysis extended to the 25-KCN union

## General idea
For each target gene, MISTy separates:
- `intra`: effect of the same spot
- `juxta_ct`: effect of the immediate neighbors
- `para_ct`: effect of the broader surrounding tissue

These analyses do not assign a gene to one cell directly. They estimate which cell-type context best explains the spatial signal.

## Typical interpretation
- fibroblast-dominant `intra`: local stromal territory
- tumor-epithelial `intra`: local tumoral territory
- strong `juxta_ct`: border or interface effect
- strong `para_ct`: broader tissue-field effect

