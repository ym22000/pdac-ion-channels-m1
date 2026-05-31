# Two-stage multicohort survival meta-analysis

This folder contains the main survival branch used as a secondary clinical opening in the project.

## Principle

The analysis is done in two stages:

1. fit one adjusted Cox model per cohort
2. pool the cohort-specific `log(HR)` estimates with a random-effects meta-analysis

This makes cohort-to-cohort variability explicit and is more interpretable than a single pooled model.

## Model

- stage 1:
  - `Surv(os_months, event) ~ gene_value + stage_group + grade_group`
- stage 2:
  - random-effects meta-analysis on `log(HR)`
  - `REML`
  - `Knapp-Hartung`

## Main script

- `../run_kcn25_multicohort_two_stage_meta_survival.R`

## Main outputs

- `tables/kcn25_multicohort_two_stage_stage1_cox.tsv`
- `tables/kcn25_multicohort_two_stage_meta.tsv`
- `tables/kcn25_multicohort_two_stage_meta_significant.tsv`
- `tables/kcn25_multicohort_two_stage_leave_one_out.tsv`
- `figures/Figure_01_KCN25_two_stage_meta_summary.pdf`
- `figures/by_gene/`
