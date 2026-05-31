# Hotspot for All KCNs in Primary Sections

This folder contains the main Hotspot run for all KCN genes available in the spatial object.

This is the main folder for spatial structure in the KCN project because it puts all KCN Hotspot results in one place.

Main script:
- `run_hotspot_all_kcns_primary.py`

Subfolders:
- `local_correlations/`: local gene-gene correlations around each KCN
- `kcnn4_territory_gsea/`: territory GSEA around KCN-high spots

Related folder:
- `../all_kcns_normal_pancreas/`: matching Hotspot, local-correlation, and normal-pancreas territory-GSEA runs in the three normal pancreas sections

Main outputs:
- one Hotspot PDF per KCN
- summary tables by gene and section

## Requirements
- `Python 3.10+`
- Main Python packages:
  `anndata`, `hotspot`, `numpy`, `pandas`, `scipy`, `matplotlib`
- Main reused inputs:
  the exported primary-section matrices and coordinates from `communication/commot/shared_exports/exports/`
- This folder is the main entry point for KCN-level Hotspot analyses in the spatial project.

