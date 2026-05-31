# Hotspot Local Correlations

This folder contains the Hotspot local-correlation analysis run for all KCN targets.

After testing whether a KCN is spatially structured, the next step is to see which genes tend to follow the same local spatial pattern.

Main script:
- `run_kcn_local_correlations_primary.py`

Main outputs:
- one PDF per target KCN
- partner summary tables
- section summary tables

## Requirements
- `Python 3.10+`
- Main Python packages:
  `anndata`, `hotspot`, `numpy`, `pandas`, `scipy`, `matplotlib`
- Main reused inputs:
  the section exports already prepared for the Hotspot run
- This step assumes the primary-section exports have already been created.

These results are useful to identify the genes that tend to co-vary locally with each KCN in space.

