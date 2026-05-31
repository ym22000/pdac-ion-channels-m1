# Marker Local Correlations Primary

This folder contains a tables-only Hotspot local-correlation run for selected
stromal inflammatory-border markers across primary PDAC sections:

- `C7`
- `C1R`
- `CCL19`
- `CCL21`
- `MFAP4`
- `CLU`
- `DCN`

The goal is to identify the genes that tend to follow the same local spatial
pattern as each marker, with a dedicated KCN-filtered summary added on top of
the complete partner table.

Main script:
- `run_marker_local_correlations_primary.py`

Outputs:
- `tables/local_correlations_by_target_and_section.tsv`
- `tables/local_correlations_partner_summary.tsv`
- `tables/marker_kcn_local_correlations_by_target_and_section.tsv`
- `tables/marker_kcn_local_correlations_summary.tsv`
- `tables/target_gene_section_summary.tsv`
- `tables/target_gene_availability.tsv`

Method settings:
- Hotspot model: `bernoulli`
- neighbors: `30`
- minimum detected spots: `10`
- local correlations computed from a section-level spatial graph
- genes retained for local-correlation matrices:
  top spatial genes plus forced inclusion of the marker set and the KCN list

