#!/usr/bin/env python
"""
Hotspot spatial autocorrelation analysis for all KCN genes in normal pancreas.

This script mirrors the main PDAC Hotspot run, but it is restricted to the
three normal pancreas sections. The goal is to check which KCN genes already
show structured spatial patterns in healthy tissue and which patterns look more
tumor-specific in PDAC.

The analysis is done section by section because spatial distances only make
sense within one Visium slide. Raw counts are kept, a spatial KNN graph is
built, and then Hotspot computes gene-wise local spatial autocorrelation.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import warnings

import anndata as ad
import hotspot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[2]
EXPORTS_DIR_DEFAULT = SCRIPT_DIR / "exports_normal_pancreas" / "exports"
OUT_DIR_DEFAULT = SCRIPT_DIR
KCN_LIST_PATH_DEFAULT = PROJECT_DIR.parent / "mycaf_icaf_correlations" / "kcn_union_25.tsv"

MODEL = "bernoulli"
N_NEIGHBORS = 30
MIN_GENE_SPOTS = 10
COLOR_UPPER_QUANTILE = 0.95


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run section-wise Hotspot spatial autocorrelation for all KCN genes in normal pancreas."
    )
    parser.add_argument(
        "--exports-dir",
        default=str(EXPORTS_DIR_DEFAULT),
        help="Directory containing exported normal pancreas ST sections.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(OUT_DIR_DEFAULT),
        help="Output directory for Hotspot results.",
    )
    parser.add_argument(
        "--kcn-list",
        default=str(KCN_LIST_PATH_DEFAULT),
        help="Path to the 25-gene KCN union table.",
    )
    parser.add_argument("--max-images", type=int, default=None, help="Optional debug limit on the number of sections.")
    return parser.parse_args()


def load_target_genes(kcn_list_path: Path) -> list[str]:
    kcn_df = pd.read_csv(kcn_list_path, sep="\t")
    genes = kcn_df["gene"].dropna().astype(str).tolist()
    return list(dict.fromkeys(genes))


def load_section(section_dir: Path) -> ad.AnnData:
    counts = spio.mmread(section_dir / "counts.mtx").tocsr()
    genes = pd.read_csv(section_dir / "genes.tsv", sep="\t", header=None)[0].tolist()
    spots = pd.read_csv(section_dir / "spots.tsv", sep="\t", header=None)[0].tolist()
    coords = pd.read_csv(section_dir / "coords.tsv", sep="\t")
    metadata = pd.read_csv(section_dir / "metadata.tsv", sep="\t")

    coords = coords.set_index("spot").loc[spots].reset_index()
    metadata = metadata.set_index("spot").loc[spots].reset_index()

    adata = ad.AnnData(X=counts.T)
    adata.var_names = genes
    adata.obs_names = spots
    adata.obs = metadata.set_index("spot")
    adata.obsm["spatial"] = coords[["imagecol", "imagerow"]].to_numpy(dtype=float)
    adata.layers["csc_counts"] = adata.X.tocsc()
    adata.obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).reshape(-1)
    return adata


def keep_genes_for_hotspot(adata: ad.AnnData, target_genes: list[str]) -> ad.AnnData:
    detected = np.asarray((adata.X > 0).sum(axis=0)).reshape(-1)
    keep = detected >= MIN_GENE_SPOTS

    for gene in target_genes:
        if gene in adata.var_names:
            keep[np.where(adata.var_names == gene)[0][0]] = True

    adata = adata[:, keep].copy()

    x = adata.X
    if sparse.issparse(x):
        x = x.toarray()
    gene_var = np.var(np.asarray(x), axis=0)
    adata = adata[:, gene_var > 0].copy()
    return adata


def compute_hotspot_autocorrelations(adata: ad.AnnData) -> pd.DataFrame:
    hs = hotspot.Hotspot(
        adata,
        layer_key="csc_counts",
        model=MODEL,
        latent_obsm_key="spatial",
        umi_counts_obs_key="total_counts",
    )
    hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS, approx_neighbors=False)
    return hs.compute_autocorrelations(jobs=1)


def dense_gene_vector(adata: ad.AnnData, gene: str) -> np.ndarray:
    x = adata[:, gene].X
    if sparse.issparse(x):
        return np.asarray(x.toarray()).reshape(-1)
    return np.asarray(x).reshape(-1)


def plot_gene_section(
    image_name: str,
    gene: str,
    coords: np.ndarray,
    expr: np.ndarray,
    hotspot_stats: pd.Series,
    ax: plt.Axes,
) -> plt.Collection:
    vmax = float(np.quantile(expr, COLOR_UPPER_QUANTILE))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = float(expr.max()) if expr.size else 1.0
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0

    points = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=expr,
        s=9,
        cmap="viridis",
        vmin=0,
        vmax=vmax,
        linewidths=0,
    )
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(
        f"{image_name}\n{gene}: C={hotspot_stats['C']:.3f}, Z={hotspot_stats['Z']:.2f}, FDR={hotspot_stats['FDR']:.2e}",
        fontsize=11,
    )
    for spine in ax.spines.values():
        spine.set_visible(False)
    return points


def save_gene_pdf(
    gene: str,
    plot_payloads: list[tuple[str, np.ndarray, np.ndarray, pd.Series]],
    pdf_path: Path,
) -> None:
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(pdf_path) as pdf:
        summary_df = pd.DataFrame(
            [
                {
                    "orig.ident": image_name,
                    "C": float(stats["C"]),
                    "Z": float(stats["Z"]),
                    "FDR": float(stats["FDR"]),
                }
                for image_name, _, _, stats in plot_payloads
            ]
        ).sort_values("Z", ascending=False)

        fig, ax = plt.subplots(figsize=(8.5, 5.5))
        ax.bar(summary_df["orig.ident"], summary_df["Z"], color="#1f77b4")
        ax.set_title(f"{gene} Hotspot Z-score across normal pancreas sections")
        ax.set_ylabel("Hotspot Z-score")
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=45)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        for image_name, coords, expr, stats in plot_payloads:
            fig, ax = plt.subplots(figsize=(5.2, 4.8))
            points = plot_gene_section(image_name, gene, coords, expr, stats, ax)
            cbar = fig.colorbar(points, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(f"{gene} raw counts")
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def build_gene_summary(summary_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene, group in summary_df.dropna(subset=["Z"]).groupby("gene", sort=True):
        rows.append(
            {
                "gene": gene,
                "n_sections_tested": int(group.shape[0]),
                "n_fdr_lt_0_05": int((group["FDR"] < 0.05).sum()),
                "median_Z": float(group["Z"].median()),
                "max_Z": float(group["Z"].max()),
                "median_C": float(group["C"].median()),
            }
        )
    return pd.DataFrame(rows).sort_values("median_Z", ascending=False)


def main() -> None:
    args = parse_args()
    exports_dir = Path(args.exports_dir)
    out_dir = Path(args.out_dir)
    kcn_list_path = Path(args.kcn_list)
    fig_dir = out_dir / "figures"
    table_dir = out_dir / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    target_genes = load_target_genes(kcn_list_path)
    section_dirs = sorted([p for p in exports_dir.iterdir() if p.is_dir()])
    if args.max_images is not None:
        section_dirs = section_dirs[: args.max_images]

    availability_rows = []
    summary_rows = []
    payload_by_gene: dict[str, list[tuple[str, np.ndarray, np.ndarray, pd.Series]]] = {gene: [] for gene in target_genes}

    for section_dir in section_dirs:
        image_name = section_dir.name
        print(f"[Hotspot normal] Processing {image_name}...")
        try:
            adata = load_section(section_dir)
            adata = keep_genes_for_hotspot(adata, target_genes)
            results = compute_hotspot_autocorrelations(adata)

            for gene in target_genes:
                available = gene in adata.var_names and gene in results.index
                expr = dense_gene_vector(adata, gene) if gene in adata.var_names else np.array([])
                availability_rows.append(
                    {
                        "orig.ident": image_name,
                        "gene": gene,
                        "available_for_hotspot": available,
                        "positive_spots": int((expr > 0).sum()) if expr.size else 0,
                    }
                )

                if not available:
                    summary_rows.append(
                        {
                            "orig.ident": image_name,
                            "gene": gene,
                            "n_spots": int(adata.n_obs),
                            "n_genes_tested": int(adata.n_vars),
                            "positive_spots": int((expr > 0).sum()) if expr.size else np.nan,
                            "C": np.nan,
                            "Z": np.nan,
                            "Pval": np.nan,
                            "FDR": np.nan,
                            "model": MODEL,
                            "n_neighbors": N_NEIGHBORS,
                            "error": "Gene unavailable in filtered Hotspot result",
                        }
                    )
                    continue

                stats = results.loc[gene]
                summary_rows.append(
                    {
                        "orig.ident": image_name,
                        "gene": gene,
                        "n_spots": int(adata.n_obs),
                        "n_genes_tested": int(adata.n_vars),
                        "positive_spots": int((expr > 0).sum()),
                        "C": float(stats["C"]),
                        "Z": float(stats["Z"]),
                        "Pval": float(stats["Pval"]),
                        "FDR": float(stats["FDR"]),
                        "model": MODEL,
                        "n_neighbors": N_NEIGHBORS,
                    }
                )
                payload_by_gene[gene].append(
                    (
                        image_name,
                        adata.obsm["spatial"].copy(),
                        expr.copy(),
                        stats.copy(),
                    )
                )
        except Exception as exc:
            warnings.warn(f"Skipping {image_name}: {exc}")
            for gene in target_genes:
                availability_rows.append(
                    {
                        "orig.ident": image_name,
                        "gene": gene,
                        "available_for_hotspot": False,
                        "positive_spots": np.nan,
                    }
                )
                summary_rows.append(
                    {
                        "orig.ident": image_name,
                        "gene": gene,
                        "n_spots": np.nan,
                        "n_genes_tested": np.nan,
                        "positive_spots": np.nan,
                        "C": np.nan,
                        "Z": np.nan,
                        "Pval": np.nan,
                        "FDR": np.nan,
                        "model": MODEL,
                        "n_neighbors": N_NEIGHBORS,
                        "error": str(exc),
                    }
                )

    availability_df = pd.DataFrame(availability_rows)
    availability_df.to_csv(table_dir / "hotspot_gene_availability.tsv", sep="\t", index=False)

    summary_df = pd.DataFrame(summary_rows).sort_values(["gene", "Z"], ascending=[True, False], na_position="last")
    summary_df.to_csv(table_dir / "hotspot_by_gene_and_section.tsv", sep="\t", index=False)

    gene_summary_df = build_gene_summary(summary_df)
    gene_summary_df.to_csv(table_dir / "hotspot_gene_summary.tsv", sep="\t", index=False)

    for gene, payloads in payload_by_gene.items():
        if not payloads:
            continue
        save_gene_pdf(
            gene=gene,
            plot_payloads=payloads,
            pdf_path=fig_dir / f"{gene.lower()}_hotspot_normal_pancreas_sections_viridis.pdf",
        )


if __name__ == "__main__":
    main()
