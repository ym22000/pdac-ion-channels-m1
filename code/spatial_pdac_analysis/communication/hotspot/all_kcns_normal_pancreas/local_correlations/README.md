# Local Correlations in Normal Pancreas Sections

This folder contains the Hotspot local correlation run for KCN genes in the
three normal pancreas sections.

This step goes one level deeper than gene-wise autocorrelation. It asks which
other genes follow the same local spatial pattern as each KCN inside normal
pancreas tissue. In practice, this helps separate resident modules from
PDAC-specific modules when the same KCN is studied in both contexts.

Main script:
- `run_kcn_local_correlations_normal_pancreas.py`

Main outputs:
- `local_correlations_by_target_and_section.tsv`
- `local_correlations_partner_summary.tsv`
- one PDF per KCN

## Requirements
- `Python 3.10+`
- Main Python packages:
  `anndata`, `hotspot`, `numpy`, `pandas`, `scipy`, `matplotlib`
- Reused inputs:
  section exports produced by `../run_export_normal_pancreas_sections.R`

