# KCN-high neighborhood analysis against dominant spot cell types

## Summary
This folder contains a simple K-nearest-neighbor spatial analysis to ask which `first_type` labels tend to sit near the spots where a KCN gene is high.

The main idea is not to draw raw spot-level KNN graphs, because those quickly become hard to read for many genes and many sections. Instead, the script summarizes the local neighborhood composition of `KCN-high` spots and compares it with the rest of the section.

This gives an easy answer to questions like:
- are `KCNN4-high` spots close to tumor epithelial spots?
- are `KCNMA1-high` spots close to `myCAF` spots?
- are `KCNJ8-high` spots close to `PVL` or immune-border territories?

## Requirements
- `Python 3.10+`
- Main Python packages:
  `numpy`, `pandas`, `scipy`, `matplotlib`
- The script expects the exported primary sections already present in:
  `communication/commot/shared_exports/exports/`
- It also expects the 25-gene KCN table in:
  `../mycaf_icaf_correlations/kcn_union_25.tsv`

## Main script
- `run_kcn_high_first_type_knn_neighborhood.py`

## Main moffitt_stromal_exploration_outputs
- `tables/kcn_high_first_type_knn_by_section.tsv`
- `tables/kcn_high_first_type_knn_summary.tsv`
- `tables/kcn_high_first_type_knn_availability.tsv`
- `tables/kcn_high_first_type_knn_availability_summary.tsv`
- `figures/kcn_high_first_type_knn_global_dotplot.pdf`
- `figures/*_first_type_knn_neighborhood.pdf`

## Method notes
- The analysis is run section by section on the primary PDAC exports.
- `KCN-high` spots are defined as the top 10% of positive spots for a gene in a section, with at least 20 spots, and only when the gene is detected in at least 30 spots.
- A `K = 15` nearest-neighbor graph is used to summarize the local neighborhood.
- For each `first_type`, the script compares the mean neighbor fraction around `KCN-high` spots with the mean neighbor fraction around the rest of the section.
- The per-gene PDF contains:
  - a barplot of median neighbor enrichment by `first_type`
  - a small network with one central `KCN-high` node and surrounding `first_type` nodes, where edge color shows enrichment or depletion and edge width shows effect size

## How to run
From this folder:

```bash
python run_kcn_high_first_type_knn_neighborhood.py
```

