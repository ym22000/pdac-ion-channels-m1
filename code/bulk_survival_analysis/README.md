# Multicohort Survival Analyses

This folder contains the current multicohort survival analyses for the 25 KCN genes across several PDAC cohorts.

Main scripts:
- `build_fused_survival_bundle.R`
- `run_kcn25_multicohort_two_stage_meta_survival.R`
- `run_kcn25_multicohort_km_within_cohort_median_survival.R`

Main subfolders:
- `integrated_survival_bundle/`: fused cohort-level inputs
- `survival_kcn25_two_stage_meta/`: cohort-wise adjusted Cox models pooled by random-effects meta-analysis
- `survival_kcn25_km_within_cohort_median/`: one 2x2 Kaplan-Meier panel per gene with a median split inside each cohort

These analyses were used to test whether the KCN candidates also have prognostic value in bulk patient cohorts.

## Recommended interpretation

The most defensible survival layer in this folder is now:
- `survival_kcn25_two_stage_meta/`

Why:
- each cohort is modeled separately first
- cohort-level adjusted log(HR) estimates are then pooled with a random-effects model
- heterogeneity and prediction intervals are explicit

The Kaplan-Meier branch is kept as a visual support layer only:
- `survival_kcn25_km_within_cohort_median/`
- one within-cohort median split per gene
- one `2 x 2` panel per gene across the four cohorts

## Requirements
- `R 4.5.x`
- Main R packages across this folder:
  `survival`, `metafor`, `dplyr`, `tidyr`, `ggplot2`
- Main input bundles:
  `PH_GSE183795_allKCN_bundle.rds`, `Hussain_GSE62452_bundle.rds`,
  `Zhang_GSE28735_bundle.rds`, `Bailey_QCMG_UQ_2016_allKCN_bundle.rds`
- The integrated run also depends on the fused bundle created by `build_fused_survival_bundle.R`.

