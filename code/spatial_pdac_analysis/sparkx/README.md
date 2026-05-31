# SPARK-X analysis for the KCN panel in primary PDAC sections

## Summary
This folder contains a SPARK-X analysis of the KCN panel in the primary PDAC Visium sections.

The main idea here is to complement the Hotspot results with another published spatial method. Hotspot already showed that some KCN genes form clear spatial patches, but it is still useful to test the same question with SPARK-X because it is a standard tool for detecting spatially variable genes in larger spatial transcriptomics datasets.

The analysis is done section by section on the exported primary Visium data. The script runs SPARK-X on all genes in each section, then keeps the results for the KCN union panel and saves one summary PDF per gene.

## Requirements
- `R 4.5.x`
- Main R packages:
  `Matrix`, `ggplot2`, `dplyr`, `tidyr`, `readr`, `purrr`, `tibble`, `SPARK`
- The script expects the primary section exports already created in:
  `communication/commot/shared_exports/exports/`
- It also expects the KCN union table in:
  `../mycaf_icaf_correlations/kcn_union_25.tsv`

## Main script
- `run_sparkx_all_kcns_primary.R`

## Main moffitt_stromal_exploration_outputs
- `tables/sparkx_by_gene_and_section.tsv`
- `tables/sparkx_gene_summary.tsv`
- `tables/sparkx_gene_availability_by_section.tsv`
- `tables/sparkx_gene_availability_summary.tsv`
- `figures/*.pdf`

## How to run
From this folder:

```r
Rscript run_sparkx_all_kcns_primary.R
```

If `Rscript` is not in the system `PATH`, run it with the full local path to the R installation.

## Method notes
- `SPARK-X` is used instead of the older `SPARK` model because these sections have more than 3,000 spots and the official SPARK documentation recommends `SPARK-X` for larger datasets.
- The analysis uses the default SPARK-X multi-kernel setting (`option = "mixture"`).
- Results are summarized with the combined p-value and the adjusted p-value returned by SPARK-X.

