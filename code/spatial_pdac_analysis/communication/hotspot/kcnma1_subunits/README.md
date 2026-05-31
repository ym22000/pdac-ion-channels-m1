# KCNMA1 Subunits Hotspot

This folder contains a focused Hotspot run for two KCNMA1-based BK subunit modules in primary PDAC:

- `KCNMA1_KCNMB1`
- `KCNMA1_KCNMB4`

The goal is to test whether spatial autocorrelation becomes clearer when the BK alpha subunit `KCNMA1` is analyzed together with a beta subunit instead of alone.

To keep the analysis directly comparable with the main primary-PDAC Hotspot run, each module is defined as a pseudo-gene built from the **sum of raw counts per spot**:

- `KCNMA1_KCNMB1 = KCNMA1 + KCNMB1`
- `KCNMA1_KCNMB4 = KCNMA1 + KCNMB4`

The workflow is the same as for the main KCN run:

1. load exported primary PDAC sections
2. create the two pseudo-genes
3. build a spatial KNN graph for each section
4. run Hotspot autocorrelation section by section
5. save per-module PDFs and summary tables

Main script:
- `run_hotspot_kcnma1_subunits_primary.py`

Outputs:
- `figures/`: one PDF per combined module
- `tables/hotspot_by_gene_and_section.tsv`
- `tables/hotspot_gene_summary.tsv`
- `tables/hotspot_gene_availability.tsv`
- `module_context_primary/`: local correlations, first-type neighborhood, and territory GSEA for the two modules

Interpretation note:
- these pseudo-genes are designed to compare tissue-level BK-related modules
  rather than to prove strict same-cell co-expression
- the context subfolder is where the two modules are compared on the same scale
  as the rest of the spatial KCN analyses

Method settings:
- Hotspot model: `bernoulli`
- neighbors: `30`
- minimum detected spots: `10`
- input counts: raw exported counts

This run is intentionally separate from `all_kcns_primary/` because it tests a specific biological hypothesis about possible BK channel spatial modules.

