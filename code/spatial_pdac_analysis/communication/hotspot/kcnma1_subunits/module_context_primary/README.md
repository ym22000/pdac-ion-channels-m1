# KCNMA1 Subunit Module Context in Primary PDAC

This folder extends the KCNMA1 BK-subunit spatial analysis beyond Hotspot
autocorrelation alone. Two summed raw-count modules are used throughout:

- `KCNMA1_KCNMB1 = KCNMA1 + KCNMB1`
- `KCNMA1_KCNMB4 = KCNMA1 + KCNMB4`

The goal is to compare their local tissue context across primary PDAC sections:

1. `local_correlations/`
   - Hotspot local correlations with other genes
   - identifies which genes share the same local spatial pattern
   - moffitt_stromal_exploration_outputs both summary tables and per-module PDF reports

2. `first_type_neighborhood/`
   - KNN neighborhood analysis around module-high spots
   - compares local `first_type` composition against the rest of the section
   - moffitt_stromal_exploration_outputs both summary tables and per-module neighborhood PDFs

3. `territory_gsea/`
   - Hallmark GSEA for module-high territories
   - compares high-module spots to the rest of each section
   - moffitt_stromal_exploration_outputs summary tables and one Hallmark GSEA PDF per module

How to rerun:
1. run the parent module Hotspot screen in `../`
2. run `run_module_local_correlations_primary.py`
3. run `run_module_first_type_neighborhood.py`
4. run `run_module_territory_gsea_primary.R`

Read this folder as a context-comparison layer:
- `KCNMA1_KCNMB1` asks whether the BK alpha subunit aligns with a stronger
  `KCNMB1`-linked stromal module
- `KCNMA1_KCNMB4` tests the same logic with `KCNMB4`
