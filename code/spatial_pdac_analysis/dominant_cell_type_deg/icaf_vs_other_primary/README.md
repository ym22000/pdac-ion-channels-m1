# iCAF vs Other Spot-Level Differential Expression in Primary PDAC

This folder contains a spot-level differential expression analysis comparing
spots annotated as `first_type == iCAF` against all other spots, section by
section across the primary PDAC Visium sections.

## Inputs

- Shared primary exports from:
  `code/spatial_pdac_analysis/communication/commot/shared_exports/exports/`
- For each section:
  - `counts.mtx`
  - `genes.tsv`
  - `spots.tsv`
  - `metadata.tsv`

## Method

For each primary section:

1. Define the target group as spots with `first_type == iCAF`.
2. Compare iCAF spots versus all non-iCAF spots.
3. For each gene, compute:
   - `mean_normexpr_icaf`
   - `mean_normexpr_other`
   - `log2fc_icaf_vs_other`
   - `pct_positive_icaf`
   - `pct_positive_other`
   - `delta_pct_positive`
   - `auc_icaf_vs_other`
   - `mannwhitney_pvalue`
   - `mannwhitney_fdr`

Normalized expression uses library-size scaling
`1e4 * counts / total_counts_per_spot`, with `log1p` values used for the
Mann-Whitney/AUC step.

## Availability rules

A section is considered available if it contains at least:

- `5` iCAF spots
- `20` non-iCAF spots

## Outputs

- `tables/icaf_vs_other_section_availability.tsv`
- `tables/icaf_vs_other_by_section.tsv`
- `tables/icaf_vs_other_gene_summary.tsv`
- `tables/icaf_vs_other_kcn_by_section.tsv`
- `tables/icaf_vs_other_kcn_summary.tsv`

## Run

```bash
python run_icaf_vs_other_primary.py
```
