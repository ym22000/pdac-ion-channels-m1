# Marker Panel Hotspot Primary

This folder contains a focused Hotspot spatial autocorrelation run for the
reference marker panel already used in the spatial projection figures, with
`POSTN` added as an extra stromal ECM marker.

The goal is to test whether these reference markers form reproducible spatial
patches in primary PDAC sections using the same framework as the main KCN
Hotspot pipeline.

Main script:
- `run_hotspot_marker_panel_primary.py`

Marker panel:
- `COL1A1`
- `FAP`
- `EPCAM`
- `KRT20` (`CK20` in figures)
- `C7`
- `C1R`
- `CCL19`
- `CCL21`
- `MFAP4`
- `CLU`
- `DCN`
- `CCN2`
- `ACTA2`
- `THBS1`
- `COL12A1`
- `POSTN`

Outputs:
- one Hotspot PDF per marker in `figures/`
- `tables/hotspot_by_marker_and_section.tsv`
- `tables/hotspot_marker_summary.tsv`
- `tables/hotspot_marker_availability.tsv`
- `tables/marker_panel.tsv`

Method settings:
- Hotspot model: `bernoulli`
- neighbors: `30`
- minimum detected spots: `10`
- input counts: raw exported counts

