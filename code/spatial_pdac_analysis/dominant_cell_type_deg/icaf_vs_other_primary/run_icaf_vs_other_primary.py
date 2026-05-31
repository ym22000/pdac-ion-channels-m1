#!/usr/bin/env python3
"""
Spot-level differential expression for iCAF spots versus all other spots
across primary PDAC sections.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import io, sparse
from scipy.stats import mannwhitneyu


SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]
EXPORT_ROOT = PROJECT_ROOT / "communication" / "commot" / "shared_exports" / "exports"
OUTPUT_DIR = SCRIPT_DIR / "tables"
KCN_LIST_PATH = PROJECT_ROOT.parent / "mycaf_icaf_correlations" / "kcn_union_25.tsv"

TARGET_FIRST_TYPE = "iCAF"
MIN_ICAF_SPOTS = 5
MIN_OTHER_SPOTS = 20
NORMALIZATION_SCALE = 1e4
LOG2FC_PSEUDOCOUNT = 1e-6


def bh_adjust(pvalues: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(pvalues), dtype=float)
    out = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    if not np.any(valid):
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranks = np.arange(1, pv.size + 1)
    adjusted = pv[order] * pv.size / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    tmp = np.empty_like(adjusted)
    tmp[order] = adjusted
    out[valid] = tmp
    return out


def load_kcn_genes(path: Path) -> list[str]:
    if not path.exists():
        return []
    df = pd.read_csv(path, sep="\t")
    if "gene" not in df.columns:
        return []
    return sorted(df["gene"].dropna().astype(str).unique().tolist())


def read_export(section_dir: Path) -> tuple[sparse.csr_matrix, list[str], list[str], pd.DataFrame]:
    counts = io.mmread(section_dir / "counts.mtx")
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    else:
        counts = counts.tocsr()

    genes = pd.read_csv(section_dir / "genes.tsv", sep="\t", header=None).iloc[:, 0].astype(str).tolist()
    spots = pd.read_csv(section_dir / "spots.tsv", sep="\t", header=None).iloc[:, 0].astype(str).tolist()

    if counts.shape == (len(genes), len(spots)):
        matrix = counts
    elif counts.shape == (len(spots), len(genes)):
        matrix = counts.T.tocsr()
    else:
        raise ValueError(
            f"Unexpected matrix shape {counts.shape} for {section_dir.name}; "
            f"expected ({len(genes)}, {len(spots)}) or ({len(spots)}, {len(genes)})."
        )

    metadata = pd.read_csv(section_dir / "metadata.tsv", sep="\t")
    metadata = metadata.set_index("spot").reindex(spots).reset_index()
    if metadata["first_type"].isna().any():
        raise ValueError(f"Missing metadata alignment for {section_dir.name}.")
    return matrix, genes, spots, metadata


def compute_section_de(section_dir: Path) -> tuple[pd.DataFrame, dict]:
    matrix, genes, spots, metadata = read_export(section_dir)
    section_name = section_dir.name

    is_icaf = metadata["first_type"].astype(str).eq(TARGET_FIRST_TYPE).to_numpy()
    n_icaf = int(is_icaf.sum())
    n_other = int((~is_icaf).sum())

    availability = {
        "orig.ident": section_name,
        "n_spots": len(spots),
        "n_icaf_spots": n_icaf,
        "n_other_spots": n_other,
        "available_for_deg": bool(n_icaf >= MIN_ICAF_SPOTS and n_other >= MIN_OTHER_SPOTS),
    }

    if n_icaf < MIN_ICAF_SPOTS or n_other < MIN_OTHER_SPOTS:
        return pd.DataFrame(), availability

    libsize = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    libsize[libsize == 0] = 1.0

    dense_counts = matrix.toarray().astype(np.float32, copy=False)
    norm = (dense_counts / libsize[None, :]) * NORMALIZATION_SCALE
    log_norm = np.log1p(norm).astype(np.float32, copy=False)

    icaf_idx = np.where(is_icaf)[0]
    other_idx = np.where(~is_icaf)[0]

    mean_icaf = norm[:, icaf_idx].mean(axis=1)
    mean_other = norm[:, other_idx].mean(axis=1)
    log2fc = np.log2((mean_icaf + LOG2FC_PSEUDOCOUNT) / (mean_other + LOG2FC_PSEUDOCOUNT))

    positive = dense_counts > 0
    pct_icaf = positive[:, icaf_idx].mean(axis=1)
    pct_other = positive[:, other_idx].mean(axis=1)
    delta_pct = pct_icaf - pct_other

    auc = np.full(len(genes), np.nan, dtype=np.float64)
    pvals = np.ones(len(genes), dtype=np.float64)

    for idx in range(len(genes)):
        x = log_norm[idx, icaf_idx]
        y = log_norm[idx, other_idx]
        if np.allclose(x, x[0]) and np.allclose(y, y[0]) and np.isclose(x[0], y[0]):
            auc[idx] = 0.5
            pvals[idx] = 1.0
            continue
        try:
            res = mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
            u_stat = float(res.statistic)
            pvals[idx] = float(res.pvalue)
            auc[idx] = u_stat / (len(x) * len(y))
        except ValueError:
            auc[idx] = 0.5
            pvals[idx] = 1.0

    fdr = bh_adjust(pvals)

    section_df = pd.DataFrame(
        {
            "orig.ident": section_name,
            "gene": genes,
            "n_icaf_spots": n_icaf,
            "n_other_spots": n_other,
            "mean_normexpr_icaf": mean_icaf,
            "mean_normexpr_other": mean_other,
            "log2fc_icaf_vs_other": log2fc,
            "pct_positive_icaf": pct_icaf,
            "pct_positive_other": pct_other,
            "delta_pct_positive": delta_pct,
            "auc_icaf_vs_other": auc,
            "mannwhitney_pvalue": pvals,
            "mannwhitney_fdr": fdr,
        }
    )
    return section_df, availability


def summarize(by_section: pd.DataFrame) -> pd.DataFrame:
    top_idx = by_section.groupby("gene")["mannwhitney_fdr"].idxmin()
    top_sections = (
        by_section.loc[top_idx, ["gene", "orig.ident", "mannwhitney_fdr"]]
        .rename(columns={"orig.ident": "top_section", "mannwhitney_fdr": "top_section_fdr"})
    )

    summary = (
        by_section.groupby("gene", as_index=False)
        .agg(
            n_sections_tested=("orig.ident", "nunique"),
            n_sections_fdr_lt_0_05=("mannwhitney_fdr", lambda s: int((s < 0.05).sum())),
            n_sections_auc_gt_0_60=("auc_icaf_vs_other", lambda s: int((s > 0.60).sum())),
            median_log2fc=("log2fc_icaf_vs_other", "median"),
            max_log2fc=("log2fc_icaf_vs_other", "max"),
            median_auc=("auc_icaf_vs_other", "median"),
            max_auc=("auc_icaf_vs_other", "max"),
            median_pct_positive_icaf=("pct_positive_icaf", "median"),
            median_pct_positive_other=("pct_positive_other", "median"),
            min_fdr=("mannwhitney_fdr", "min"),
        )
        .merge(top_sections, on="gene", how="left")
        .sort_values(["n_sections_fdr_lt_0_05", "median_log2fc", "median_auc"], ascending=[False, False, False])
    )
    return summary


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    kcn_genes = set(load_kcn_genes(KCN_LIST_PATH))

    all_sections: list[pd.DataFrame] = []
    availability_rows: list[dict] = []

    for section_dir in sorted([p for p in EXPORT_ROOT.iterdir() if p.is_dir()]):
        print(f"[iCAF DEG] Processing {section_dir.name}...")
        section_df, availability = compute_section_de(section_dir)
        availability_rows.append(availability)
        if not section_df.empty:
            all_sections.append(section_df)

    availability_df = pd.DataFrame(availability_rows)
    availability_df.to_csv(OUTPUT_DIR / "icaf_vs_other_section_availability.tsv", sep="\t", index=False)

    if not all_sections:
        raise RuntimeError("No primary sections passed the iCAF availability threshold.")

    by_section = pd.concat(all_sections, ignore_index=True)
    by_section.to_csv(OUTPUT_DIR / "icaf_vs_other_by_section.tsv", sep="\t", index=False)

    summary = summarize(by_section)
    summary.to_csv(OUTPUT_DIR / "icaf_vs_other_gene_summary.tsv", sep="\t", index=False)

    kcn_by_section = by_section[by_section["gene"].isin(kcn_genes)].copy()
    kcn_by_section.to_csv(OUTPUT_DIR / "icaf_vs_other_kcn_by_section.tsv", sep="\t", index=False)

    kcn_summary = summary[summary["gene"].isin(kcn_genes)].copy()
    kcn_summary.to_csv(OUTPUT_DIR / "icaf_vs_other_kcn_summary.tsv", sep="\t", index=False)

    print("iCAF vs other primary spot-level DEG complete.")


if __name__ == "__main__":
    main()
