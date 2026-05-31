# SIGMAR1 Hotspot Primary

This folder contains a focused Hotspot spatial autocorrelation run for `SIGMAR1`
across the primary PDAC sections used in the main spatial analyses.

The goal is to test whether `SIGMAR1` forms reproducible local spatial patches
in tumor tissue, using the same framework as the main KCN Hotspot pipeline so
the result stays directly comparable.

Main script:
- `run_hotspot_sigmar1_primary.py`

Outputs:
- `figures/sigmar1_hotspot_primary_sections_viridis.pdf`
- `tables/hotspot_by_gene_and_section.tsv`
- `tables/hotspot_gene_summary.tsv`
- `tables/hotspot_gene_availability.tsv`
- `sigmar1_context_primary/` for local correlations, `first_type` neighborhood,
  and territory GSEA around `SIGMAR1-high` spots

Interpretation note:
- this folder tests whether `SIGMAR1` behaves like a robust spatial patch-forming
  gene in primary PDAC
- the context folder is intentionally separate so the Hotspot screen can stay
  readable while still allowing a second-pass biological interpretation

Method settings:
- Hotspot model: `bernoulli`
- neighbors: `30`
- minimum detected spots: `10`
- input counts: raw exported counts

