# Kaplan-Meier within-cohort median split

This folder contains a simple Kaplan-Meier visualization layer for the 25 KCN genes.

## Principle

For each gene:
- split samples at the median within each cohort
- draw one Kaplan-Meier curve per cohort
- arrange the four cohorts in a fixed `2 x 2` layout

This branch is used for visualization only.

## Main script

- `../run_kcn25_multicohort_km_within_cohort_median_survival.R`

## Main outputs

- `tables/kcn25_multicohort_km_within_cohort_median.tsv`
- `tables/kcn25_multicohort_km_within_cohort_median_gene_summary.tsv`
- `figures/<GENE>_km_within_cohort_median.pdf`
