#!/usr/bin/env python
"""
Hotspot local correlation analysis for KCN genes in primary PDAC sections.

This script extends the Hotspot autocorrelation step by looking for genes that
follow a similar local spatial pattern as each KCN gene of interest inside
primary PDAC sections. The goal is to identify candidate niche markers and
pathway-related genes that repeatedly co-vary with each KCN in space.

The analysis is done section by section because spatial distances only make
sense within the same Visium slide. The same Hotspot settings are used across
sections so the moffitt_stromal_exploration_outputs stay easy to compare.

Method resources
    - Hotspot documentation:
      https://hotspot.readthedocs.io/en/latest/
    - Spatial tutorial:
      https://hotspot.readthedocs.io/en/latest/Spatial_Tutorial.html
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
PROJECT_DIR = SCRIPT_DIR.parents[3]
EXPORTS_DIR_DEFAULT = PROJECT_DIR / "communication" / "commot" / "shared_exports" / "exports"
KCN_LIST_PATH_DEFAULT = PROJECT_DIR.parent / "mycaf_icaf_correlations" / "kcn_union_25.tsv"
OUT_DIR_DEFAULT = SCRIPT_DIR

MODEL = "bernoulli"
N_NEIGHBORS = 30
MIN_GENE_SPOTS = 10
AUTOCORR_FDR_THRESHOLD = 0.05
MAX_SPATIAL_GENES = 400
TOP_N_GENES_IN_FIGURE = 20


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run section-wise Hotspot local correlations for all KCN genes."
    )
    parser.add_argument(
        "--exports-dir",
        default=str(EXPORTS_DIR_DEFAULT),
        help="Directory containing exported primary ST sections.",
    )
    parser.add_argument(
        "--kcn-list",
        default=str(KCN_LIST_PATH_DEFAULT),
        help="Path to the 25-gene KCN union table.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(OUT_DIR_DEFAULT),
        help="Output directory for local correlation results.",
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


def build_hotspot(adata: ad.AnnData) -> hotspot.Hotspot:
    hs = hotspot.Hotspot(
        adata,
        layer_key="csc_counts",
        model=MODEL,
        latent_obsm_key="spatial",
        umi_counts_obs_key="total_counts",
    )
    hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS, approx_neighbors=False)
    return hs


def select_genes_for_local_correlations(autocorr: pd.DataFrame, target_genes: list[str]) -> list[str]:
    spatial_genes = autocorr.loc[
        (autocorr["FDR"] < AUTOCORR_FDR_THRESHOLD) & (autocorr["C"] > 0)
    ].sort_values("Z", ascending=False)

    genes = spatial_genes.index.tolist()[:MAX_SPATIAL_GENES]
    genes.extend([gene for gene in target_genes if gene in autocorr.index])
    return list(dict.fromkeys(genes))


def extract_target_profile(local_corr_z: pd.DataFrame, image_name: str, target_gene: str) -> pd.DataFrame:
    profile = (
        local_corr_z.loc[target_gene]
        .drop(labels=[target_gene], errors="ignore")
        .rename("local_corr_z")
        .reset_index()
        .rename(columns={"index": "partner_gene"})
        .sort_values("local_corr_z", ascending=False)
    )
    profile["orig.ident"] = image_name
    profile["target_gene"] = target_gene
    return profile[["target_gene", "orig.ident", "partner_gene", "local_corr_z"]]


def build_partner_summary(profile_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (target_gene, partner_gene), group in profile_df.groupby(["target_gene", "partner_gene"], sort=True):
        rows.append(
            {
                "target_gene": target_gene,
                "partner_gene": partner_gene,
                "n_sections_tested": int(group.shape[0]),
                "n_sections_positive_z": int((group["local_corr_z"] > 0).sum()),
                "median_local_corr_z": float(group["local_corr_z"].median()),
                "max_local_corr_z": float(group["local_corr_z"].max()),
                "min_local_corr_z": float(group["local_corr_z"].min()),
            }
        )

    return pd.DataFrame(rows).sort_values(
        ["target_gene", "n_sections_positive_z", "median_local_corr_z"],
        ascending=[True, False, False],
    )


def save_target_pdf(
    target_gene: str,
    target_section_df: pd.DataFrame,
    target_summary_df: pd.DataFrame,
    pdf_path: Path,
) -> None:
    from matplotlib.backends.backend_pdf import PdfPages

    top_positive = target_summary_df.head(TOP_N_GENES_IN_FIGURE).iloc[::-1]
    shared_genes = top_positive["partner_gene"].tolist()

    heatmap_df = (
        target_section_df[target_section_df["partner_gene"].isin(shared_genes)]
        .pivot(index="partner_gene", columns="orig.ident", values="local_corr_z")
        .reindex(shared_genes)
    )

    with PdfPages(pdf_path) as pdf:
        fig, ax = plt.subplots(figsize=(8.5, 5.5))
        ax.barh(
            top_positive["partner_gene"],
            top_positive["median_local_corr_z"],
            color=plt.cm.viridis(np.linspace(0.2, 0.9, top_positive.shape[0])),
        )
        ax.set_title(f"Top genes locally correlated with {target_gene}")
        ax.set_xlabel("Median local correlation Z-score")
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8.5, 6.5))
        im = ax.imshow(heatmap_df.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-8, vmax=8)
        ax.set_title(f"{target_gene} local correlation Z-scores across sections")
        ax.set_xticks(np.arange(heatmap_df.shape[1]))
        ax.set_xticklabels(heatmap_df.columns.tolist(), rotation=45, ha="right")
        ax.set_yticks(np.arange(heatmap_df.shape[0]))
        ax.set_yticklabels(heatmap_df.index.tolist())
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Local correlation Z-score")
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        for image_name, profile in target_section_df.groupby("orig.ident", sort=True):
            profile = profile.sort_values("local_corr_z", ascending=False)
            top_pos = profile.head(TOP_N_GENES_IN_FIGURE).iloc[::-1]
            top_neg = profile.tail(TOP_N_GENES_IN_FIGURE)

            fig, axes = plt.subplots(1, 2, figsize=(11, 7))
            axes[0].barh(
                top_pos["partner_gene"],
                top_pos["local_corr_z"],
                color=plt.cm.viridis(np.linspace(0.2, 0.9, top_pos.shape[0])),
            )
            axes[0].set_title(f"{image_name}: top positive")
            axes[0].set_xlabel("Local correlation Z-score")

            axes[1].barh(
                top_neg["partner_gene"],
                top_neg["local_corr_z"],
                color=plt.cm.magma(np.linspace(0.2, 0.9, top_neg.shape[0])),
            )
            axes[1].set_title(f"{image_name}: top negative")
            axes[1].set_xlabel("Local correlation Z-score")

            fig.suptitle(f"{target_gene} local spatial correlations in {image_name}")
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def main() -> None:
    args = parse_args()
    exports_dir = Path(args.exports_dir)
    kcn_list_path = Path(args.kcn_list)
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    table_dir = out_dir / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    target_genes = load_target_genes(kcn_list_path)
    section_dirs = sorted([p for p in exports_dir.iterdir() if p.is_dir()])
    if args.max_images is not None:
        section_dirs = section_dirs[: args.max_images]

    profile_rows = []
    target_meta_rows = []
    availability_rows = []

    for section_dir in section_dirs:
        image_name = section_dir.name
        print(f"[Hotspot local corr] Processing {image_name}...")
        try:
            adata = keep_genes_for_hotspot(load_section(section_dir), target_genes)
            hs = build_hotspot(adata)
            autocorr = hs.compute_autocorrelations(jobs=1)
            genes_for_corr = select_genes_for_local_correlations(autocorr, target_genes)
            local_corr_z = hs.compute_local_correlations(genes_for_corr, jobs=1)

            n_spatial_genes = int(((autocorr["FDR"] < AUTOCORR_FDR_THRESHOLD) & (autocorr["C"] > 0)).sum())

            for target_gene in target_genes:
                available = target_gene in local_corr_z.index
                availability_rows.append(
                    {
                        "target_gene": target_gene,
                        "orig.ident": image_name,
                        "available_for_local_corr": available,
                    }
                )

                if not available:
                    target_meta_rows.append(
                        {
                            "target_gene": target_gene,
                            "orig.ident": image_name,
                            "n_spots": int(adata.n_obs),
                            "n_genes_after_filtering": int(adata.n_vars),
                            "n_spatial_genes_fdr_lt_0_05": n_spatial_genes,
                            "n_genes_used_for_local_corr": int(len(genes_for_corr)),
                            "hotspot_C": np.nan,
                            "hotspot_Z": np.nan,
                            "hotspot_FDR": np.nan,
                            "error": "Target gene unavailable in local correlation matrix",
                        }
                    )
                    continue

                profile_rows.append(extract_target_profile(local_corr_z, image_name, target_gene))
                target_meta_rows.append(
                    {
                        "target_gene": target_gene,
                        "orig.ident": image_name,
                        "n_spots": int(adata.n_obs),
                        "n_genes_after_filtering": int(adata.n_vars),
                        "n_spatial_genes_fdr_lt_0_05": n_spatial_genes,
                        "n_genes_used_for_local_corr": int(len(genes_for_corr)),
                        "hotspot_C": float(autocorr.loc[target_gene, "C"]) if target_gene in autocorr.index else np.nan,
                        "hotspot_Z": float(autocorr.loc[target_gene, "Z"]) if target_gene in autocorr.index else np.nan,
                        "hotspot_FDR": float(autocorr.loc[target_gene, "FDR"]) if target_gene in autocorr.index else np.nan,
                    }
                )
        except Exception as exc:
            warnings.warn(f"Skipping {image_name}: {exc}")
            for target_gene in target_genes:
                availability_rows.append(
                    {
                        "target_gene": target_gene,
                        "orig.ident": image_name,
                        "available_for_local_corr": False,
                    }
                )
                target_meta_rows.append(
                    {
                        "target_gene": target_gene,
                        "orig.ident": image_name,
                        "n_spots": np.nan,
                        "n_genes_after_filtering": np.nan,
                        "n_spatial_genes_fdr_lt_0_05": np.nan,
                        "n_genes_used_for_local_corr": np.nan,
                        "hotspot_C": np.nan,
                        "hotspot_Z": np.nan,
                        "hotspot_FDR": np.nan,
                        "error": str(exc),
                    }
                )

    if not profile_rows:
        raise RuntimeError("No valid local correlation results were produced.")

    profile_df = pd.concat(profile_rows, ignore_index=True)
    partner_summary_df = build_partner_summary(profile_df)
    target_meta_df = pd.DataFrame(target_meta_rows)
    availability_df = pd.DataFrame(availability_rows)

    availability_df.to_csv(table_dir / "target_gene_availability.tsv", sep="\t", index=False)
    profile_df.to_csv(table_dir / "local_correlations_by_target_and_section.tsv", sep="\t", index=False)
    partner_summary_df.to_csv(table_dir / "local_correlations_partner_summary.tsv", sep="\t", index=False)
    target_meta_df.to_csv(table_dir / "target_gene_section_summary.tsv", sep="\t", index=False)

    for target_gene in target_genes:
        target_section_df = profile_df.loc[profile_df["target_gene"] == target_gene].copy()
        target_summary_df = partner_summary_df.loc[partner_summary_df["target_gene"] == target_gene].copy()
        if target_section_df.empty or target_summary_df.empty:
            continue
        save_target_pdf(
            target_gene=target_gene,
            target_section_df=target_section_df,
            target_summary_df=target_summary_df,
            pdf_path=fig_dir / f"{target_gene.lower()}_local_correlations_primary_sections.pdf",
        )


if __name__ == "__main__":
    main()
