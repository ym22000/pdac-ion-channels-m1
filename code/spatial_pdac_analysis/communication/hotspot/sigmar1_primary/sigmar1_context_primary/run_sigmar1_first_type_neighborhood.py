#!/usr/bin/env python3
"""
Neighborhood analysis for SIGMAR1-high spots in primary PDAC.
"""

from __future__ import annotations

from pathlib import Path
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import io as spio
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR
while PROJECT_DIR.name != "spatial_pdac_analysis":
    if PROJECT_DIR.parent == PROJECT_DIR:
        raise RuntimeError("Could not find spatial_pdac_analysis root.")
    PROJECT_DIR = PROJECT_DIR.parent

EXPORTS_DIR = PROJECT_DIR / "communication" / "commot" / "shared_exports" / "exports"
OUT_DIR = SCRIPT_DIR / "first_type_neighborhood"
TABLES_DIR = OUT_DIR / "tables"
FIGURES_DIR = OUT_DIR / "figures"

TARGET_GENE = "SIGMAR1"
K_NEIGHBORS = 15
TOP_FRACTION = 0.10
MIN_POSITIVE_SPOTS = 30
MIN_HIGH_SPOTS = 20
MIN_CELLTYPE_TOTAL_SPOTS = 20
MIN_SECTIONS_FOR_TEST = 3
MAX_CELLTYPES_IN_PLOT = 12
EPSILON = 1e-4
NETWORK_RADIUS = 3.2


def bh_adjust(pvalues: pd.Series) -> pd.Series:
    p = pvalues.to_numpy(dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty(n, dtype=float)
    running_min = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        value = ranked[i] * n / rank
        running_min = min(running_min, value)
        adjusted[i] = running_min
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0, 1)
    return pd.Series(out, index=pvalues.index)


def load_section(section_dir: Path) -> dict:
    counts = spio.mmread(section_dir / "counts.mtx").tocsr()
    genes = pd.read_csv(section_dir / "genes.tsv", sep="\t", header=None)[0].tolist()
    spots = pd.read_csv(section_dir / "spots.tsv", sep="\t", header=None)[0].tolist()
    coords = pd.read_csv(section_dir / "coords.tsv", sep="\t")
    metadata = pd.read_csv(section_dir / "metadata.tsv", sep="\t")
    row_index = {gene: idx for idx, gene in enumerate(genes)}
    coords = coords.set_index("spot").loc[spots].reset_index()
    metadata = metadata.set_index("spot").loc[spots].reset_index()
    total_counts = np.asarray(counts.sum(axis=0)).reshape(-1)
    return {
        "counts": counts,
        "row_index": row_index,
        "coords": coords,
        "metadata": metadata,
        "total_counts": total_counts,
    }


def collect_global_celltypes(section_dirs: list[Path]) -> list[str]:
    total_counts = {}
    for section_dir in section_dirs:
        meta = pd.read_csv(section_dir / "metadata.tsv", sep="\t")
        for label, count in meta["first_type"].value_counts().items():
            total_counts[label] = total_counts.get(label, 0) + int(count)
    keep = [label for label, count in total_counts.items() if count >= MIN_CELLTYPE_TOTAL_SPOTS]
    keep.sort(key=lambda x: (-total_counts[x], x))
    return keep


def collect_global_celltype_counts(section_dirs: list[Path]) -> dict[str, int]:
    total_counts: dict[str, int] = {}
    for section_dir in section_dirs:
        meta = pd.read_csv(section_dir / "metadata.tsv", sep="\t")
        for label, count in meta["first_type"].value_counts().items():
            total_counts[label] = total_counts.get(label, 0) + int(count)
    return total_counts


def compute_neighbor_fractions(coords: np.ndarray, labels: pd.Series, selected_celltypes: list[str]) -> np.ndarray:
    n_spots = coords.shape[0]
    k_eff = min(K_NEIGHBORS + 1, n_spots)
    tree = cKDTree(coords)
    _, idx = tree.query(coords, k=k_eff)
    if k_eff == 1:
        raise RuntimeError("A section needs at least two spots for KNN analysis.")
    neighbor_idx = idx[:, 1:]
    onehot = pd.get_dummies(labels).reindex(columns=selected_celltypes, fill_value=0).to_numpy(dtype=float)
    return onehot[neighbor_idx].mean(axis=1)


def define_high_spots(section: dict) -> tuple[np.ndarray | None, int]:
    if TARGET_GENE not in section["row_index"]:
        return None, 0
    counts = section["counts"].getrow(section["row_index"][TARGET_GENE]).toarray().reshape(-1)
    positive_mask = counts > 0
    n_positive = int(positive_mask.sum())
    if n_positive < MIN_POSITIVE_SPOTS:
        return None, n_positive
    expr = np.log1p(1e4 * counts / np.maximum(section["total_counts"], 1))
    positive_idx = np.where(positive_mask)[0]
    n_high = max(MIN_HIGH_SPOTS, int(math.ceil(TOP_FRACTION * n_positive)))
    n_high = min(n_high, len(positive_idx))
    if n_high < MIN_HIGH_SPOTS:
        return None, n_positive
    ordered_idx = positive_idx[np.argsort(expr[positive_idx])[::-1]]
    high_idx = ordered_idx[:n_high]
    high_mask = np.zeros_like(positive_mask, dtype=bool)
    high_mask[high_idx] = True
    if (~high_mask).sum() == 0:
        return None, n_positive
    return high_mask, n_positive


def run_section(section_name: str, section_dir: Path, selected_celltypes: list[str]) -> tuple[pd.DataFrame, dict]:
    print(f"[SIGMAR1 neighborhood] Processing {section_name}...")
    section = load_section(section_dir)
    coords = section["coords"][["imagecol", "imagerow"]].to_numpy(dtype=float)
    labels = section["metadata"]["first_type"]
    neighbor_fractions = compute_neighbor_fractions(coords, labels, selected_celltypes)

    high_mask, n_positive = define_high_spots(section)
    availability = {
        "target_gene": TARGET_GENE,
        "orig.ident": section_name,
        "available_for_neighborhood": high_mask is not None,
        "n_positive_spots": n_positive,
    }
    if high_mask is None:
        return pd.DataFrame(), availability

    other_mask = ~high_mask
    high_mean = neighbor_fractions[high_mask].mean(axis=0)
    other_mean = neighbor_fractions[other_mask].mean(axis=0)
    delta = high_mean - other_mean
    ratio = (high_mean + EPSILON) / (other_mean + EPSILON)

    rows = []
    for i, celltype in enumerate(selected_celltypes):
        rows.append(
            {
                "target_gene": TARGET_GENE,
                "orig.ident": section_name,
                "first_type": celltype,
                "n_high_spots": int(high_mask.sum()),
                "n_other_spots": int(other_mask.sum()),
                "high_neighbor_fraction": float(high_mean[i]),
                "other_neighbor_fraction": float(other_mean[i]),
                "delta_neighbor_fraction": float(delta[i]),
                "enrichment_ratio": float(ratio[i]),
            }
        )
    return pd.DataFrame(rows), availability


def summarize_results(section_results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for celltype, group in section_results.groupby("first_type", sort=True):
        deltas = group["delta_neighbor_fraction"].to_numpy(dtype=float)
        if group.shape[0] >= MIN_SECTIONS_FOR_TEST and not np.allclose(deltas, 0):
            try:
                test = wilcoxon(deltas, alternative="two-sided", zero_method="wilcox")
                pval = float(test.pvalue)
                stat = float(test.statistic)
            except ValueError:
                pval = 1.0
                stat = np.nan
        else:
            pval = np.nan
            stat = np.nan
        rows.append(
            {
                "target_gene": TARGET_GENE,
                "first_type": celltype,
                "n_sections_used": int(group.shape[0]),
                "n_sections_positive_delta": int((group["delta_neighbor_fraction"] > 0).sum()),
                "median_high_neighbor_fraction": float(group["high_neighbor_fraction"].median()),
                "median_other_neighbor_fraction": float(group["other_neighbor_fraction"].median()),
                "median_delta_neighbor_fraction": float(group["delta_neighbor_fraction"].median()),
                "median_enrichment_ratio": float(group["enrichment_ratio"].median()),
                "wilcoxon_statistic": stat,
                "wilcoxon_pvalue": pval,
            }
        )
    summary = pd.DataFrame(rows)
    if summary.empty:
        return summary
    summary["wilcoxon_fdr_global"] = np.nan
    valid = summary["wilcoxon_pvalue"].notna()
    if valid.any():
        summary.loc[valid, "wilcoxon_fdr_global"] = bh_adjust(summary.loc[valid, "wilcoxon_pvalue"])
    return summary.sort_values("median_delta_neighbor_fraction", ascending=False)


def draw_neighborhood_network(
    ax: plt.Axes,
    ordered: pd.DataFrame,
    celltype_sizes: dict[str, int],
) -> None:
    ordered = ordered.copy()
    ordered["abs_delta"] = ordered["median_delta_neighbor_fraction"].abs()
    max_abs_delta = max(float(ordered["abs_delta"].max()), 1e-6)
    max_count = max(float(max(celltype_sizes.values())), 1.0)

    ax.scatter([0], [0], s=2600, color="#2b2b2b", edgecolor="black", linewidth=1.0, zorder=3)
    ax.text(0, 0, "SIGMAR1-high", ha="center", va="center", color="white", fontsize=10, fontweight="bold", zorder=4)

    angles = np.linspace(0, 2 * np.pi, num=len(ordered), endpoint=False)
    for angle, (_, row) in zip(angles, ordered.iterrows()):
        label = row["first_type"]
        delta = float(row["median_delta_neighbor_fraction"])
        x = NETWORK_RADIUS * np.cos(angle)
        y = NETWORK_RADIUS * np.sin(angle)
        edge_color = "#d62728" if delta >= 0 else "#1f77b4"
        line_width = 1.2 + 9.0 * abs(delta) / max_abs_delta
        node_size = 400 + 2600 * celltype_sizes.get(label, 0) / max_count

        ax.plot([0, x], [0, y], color=edge_color, linewidth=line_width, alpha=0.85, zorder=1)
        ax.scatter([x], [y], s=node_size, color="#f5f5f5", edgecolor=edge_color, linewidth=2.0, zorder=2)
        ax.text(x, y, label, ha="center", va="center", fontsize=8.5, zorder=3)

    ax.set_title("SIGMAR1: high to first_type neighborhood network")
    ax.text(
        0,
        -NETWORK_RADIUS - 1.15,
        "Edge width = |median delta neighbor fraction|   |   Red = enriched   |   Blue = depleted",
        ha="center",
        va="center",
        fontsize=8.5,
    )
    ax.set_xlim(-NETWORK_RADIUS - 2.0, NETWORK_RADIUS + 2.0)
    ax.set_ylim(-NETWORK_RADIUS - 1.6, NETWORK_RADIUS + 1.6)
    ax.axis("off")


def save_target_pdf(summary: pd.DataFrame, pdf_path: Path, celltype_sizes: dict[str, int]) -> None:
    from matplotlib.backends.backend_pdf import PdfPages

    ordered = summary.sort_values("median_delta_neighbor_fraction", ascending=False).head(MAX_CELLTYPES_IN_PLOT)

    with PdfPages(pdf_path) as pdf:
        fig, ax = plt.subplots(figsize=(9.5, 6.0))
        colors = ["#d62728" if x >= 0 else "#1f77b4" for x in ordered["median_delta_neighbor_fraction"]]
        ax.barh(ordered["first_type"][::-1], ordered["median_delta_neighbor_fraction"][::-1], color=colors[::-1])
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_title("SIGMAR1: KNN neighborhood enrichment by cell type")
        ax.set_xlabel("Median delta neighbor fraction (SIGMAR1-high - other spots)")
        ax.set_ylabel("first_type")
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(9.5, 8.0))
        draw_neighborhood_network(ax=ax, ordered=ordered, celltype_sizes=celltype_sizes)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    section_dirs = sorted([p for p in EXPORTS_DIR.iterdir() if p.is_dir()])
    selected_celltypes = collect_global_celltypes(section_dirs)
    celltype_sizes = collect_global_celltype_counts(section_dirs)

    section_frames = []
    availability_rows = []
    for section_dir in section_dirs:
        section_df, availability = run_section(section_dir.name, section_dir, selected_celltypes)
        section_frames.append(section_df)
        availability_rows.append(availability)

    by_section = pd.concat(section_frames, ignore_index=True) if section_frames else pd.DataFrame()
    by_section.to_csv(TABLES_DIR / "sigmar1_first_type_knn_by_section.tsv", sep="\t", index=False)

    summary = summarize_results(by_section) if not by_section.empty else pd.DataFrame()
    summary.to_csv(TABLES_DIR / "sigmar1_first_type_knn_summary.tsv", sep="\t", index=False)

    availability_df = pd.DataFrame(availability_rows)
    availability_df.to_csv(TABLES_DIR / "sigmar1_first_type_knn_availability.tsv", sep="\t", index=False)

    if not summary.empty:
        save_target_pdf(
            summary=summary,
            pdf_path=FIGURES_DIR / "sigmar1_first_type_knn_neighborhood.pdf",
            celltype_sizes=celltype_sizes,
        )

    print("SIGMAR1 first_type neighborhood complete.")


if __name__ == "__main__":
    main()
