# USER KCN25 MAST by Cell Subset

This folder contains a compact MAST analysis of the 25 KCN union genes across
the major USER cell subsets.

Goal:
- ask where each KCN is enriched in the global USER atlas
- keep the analysis restricted to the 25 KCN genes only
- use a simple one-vs-rest MAST design for each `cell_subsets` label

Main script:
- `run_user_kcn25_mast_cell_subsets.R`

Inputs:
- `../inputs/gene_sorted-naivedata_scp.mtx`
- `../inputs/naivedata_scp.genes.csv`
- `../inputs/naivedata_scp.barcodes.csv`
- `../inputs/combinenaivedata-reprocessed-clean-detailed-annotations.tsv`
- `../../mycaf_icaf_correlations/kcn_union_25.tsv`

Main outputs:
- `tables/USER_kcn25_mast_by_cell_subset.tsv`
- `tables/USER_kcn25_mast_by_cell_subset_significant.tsv`
- `tables/USER_kcn25_cell_subset_counts.tsv`
- `tables/USER_kcn25_missing_in_user.tsv`

Minimal output columns:
- `cell_subset`
- `KCN`
- `log2FC`
- `FDR`
- `Lecture`

Method:
- build a compact Seurat object restricted to the 25 KCN genes
- normalize with `LogNormalize`
- run `FindAllMarkers(..., test.use = "MAST")`
- compare each USER `cell_subsets` group against the rest
- adjust p-values with the default Seurat/MAST workflow

