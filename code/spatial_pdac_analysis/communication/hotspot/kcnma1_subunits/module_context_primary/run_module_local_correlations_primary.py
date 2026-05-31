#!/usr/bin/env python3
"""
Hotspot local correlations for KCNMA1-based subunit modules in primary PDAC.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import hotspot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR
while PROJECT_DIR.name != "spatial_pdac_analysis":
    if PROJECT_DIR.parent == PROJECT_DIR:
        raise RuntimeError("Could not find spatial_pdac_analysis root.")
    PROJECT_DIR = PROJECT_DIR.parent

EXPORTS_DIR = PROJECT_DIR / "communication" / "commot" / "shared_exports" / "exports"
KCN_LIST_PATH = PROJECT_DIR.parent / "mycaf_icaf_correlations" / "kcn_union_25.tsv"
OUT_DIR = SCRIPT_DIR / "local_correlations"
TABLES_DIR = OUT_DIR / "tables"
FIGURES_DIR = OUT_DIR / "figures"

MODEL = "bernoulli"
N_NEIGHBORS = 30
MIN_GENE_SPOTS = 10
AUTOCORR_FDR_THRESHOLD = 0.05
MAX_SPATIAL_GENES = 400
TOP_N_GENES_IN_FIGURE = 15

TARGET_MODULES = {
    "KCNMA1_KCNMB1": ("KCNMA1", "KCNMB1"),
    "KCNMA1_KCNMB4": ("KCNMA1", "KCNMB4"),
}


def load_kcn_genes() -> list[str]:
    tbl = pd.read_csv(KCN_LIST_PATH, sep="\t")
    return tbl["gene"].dropna().astype(str).drop_duplicates().tolist()


def append_module_rows(counts: sparse.csr_matrix, genes: list[str]) -> tuple[sparse.csr_matrix, list[str]]:
    gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}
    extra_rows = []
    extra_genes = []
    for module_name, components in TARGET_MODULES.items():
        if not all(gene in gene_to_idx for gene in components):
            continue
        row = counts.getrow(gene_to_idx[components[0]]) + counts.getrow(gene_to_idx[components[1]])
        extra_rows.append(row.tocsr())
        extra_genes.append(module_name)
    if not extra_rows:
        return counts, genes
    combined = sparse.vstack([counts] + extra_rows, format="csr")
    return combined, genes + extra_genes


def load_section(section_dir: Path) -> ad.AnnData:
    counts = spio.mmread(section_dir / "counts.mtx").tocsr()
    genes = pd.read_csv(section_dir / "genes.tsv", sep="\t", header=None)[0].tolist()
    spots = pd.read_csv(section_dir / "spots.tsv", sep="\t", header=None)[0].tolist()
    coords = pd.read_csv(section_dir / "coords.tsv", sep="\t")
    metadata = pd.read_csv(section_dir / "metadata.tsv", sep="\t")

    counts, genes = append_module_rows(counts, genes)
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
    x = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
    gene_var = np.var(x, axis=0)
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


def select_genes_for_local_correlations(autocorr: pd.DataFrame, target_genes: list[str], forced_genes: list[str]) -> list[str]:
    spatial_genes = autocorr.loc[
        (autocorr["FDR"] < AUTOCORR_FDR_THRESHOLD) & (autocorr["C"] > 0)
    ].sort_values("Z", ascending=False)
    genes = spatial_genes.index.tolist()[:MAX_SPATIAL_GENES]
    genes.extend([gene for gene in target_genes if gene in autocorr.index])
    genes.extend([gene for gene in forced_genes if gene in autocorr.index])
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
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    target_genes = list(TARGET_MODULES)
    kcn_genes = load_kcn_genes()

    profile_rows = []
    availability_rows = []
    section_summary_rows = []

    for section_dir in sorted([p for p in EXPORTS_DIR.iterdir() if p.is_dir()]):
        image_name = section_dir.name
        print(f"[Module local corr] Processing {image_name}...")
        adata = load_section(section_dir)
        adata = keep_genes_for_hotspot(adata, target_genes)

        present_targets = [gene for gene in target_genes if gene in adata.var_names]
        for gene in target_genes:
            availability_rows.append(
                {"target_gene": gene, "orig.ident": image_name, "available_for_local_corr": gene in adata.var_names}
            )
        if not present_targets:
            continue

        hs = build_hotspot(adata)
        autocorr = hs.compute_autocorrelations(jobs=1)
        genes_for_corr = select_genes_for_local_correlations(autocorr, present_targets, kcn_genes)
        local_corr_z = hs.compute_local_correlations(genes_for_corr, jobs=1)

        for target_gene in present_targets:
            profile_rows.append(extract_target_profile(local_corr_z, image_name, target_gene))
            stats = autocorr.loc[target_gene]
            section_summary_rows.append(
                {
                    "target_gene": target_gene,
                    "orig.ident": image_name,
                    "n_spots": int(adata.n_obs),
                    "n_genes_after_filtering": int(adata.n_vars),
                    "n_spatial_genes_fdr_lt_0_05": int(((autocorr["FDR"] < AUTOCORR_FDR_THRESHOLD) & (autocorr["C"] > 0)).sum()),
                    "n_genes_used_for_local_corr": int(len(genes_for_corr)),
                    "hotspot_C": float(stats["C"]),
                    "hotspot_Z": float(stats["Z"]),
                    "hotspot_FDR": float(stats["FDR"]),
                }
            )

    availability_df = pd.DataFrame(availability_rows)
    availability_df.to_csv(TABLES_DIR / "target_gene_availability.tsv", sep="\t", index=False)

    profile_df = pd.concat(profile_rows, ignore_index=True) if profile_rows else pd.DataFrame(
        columns=["target_gene", "orig.ident", "partner_gene", "local_corr_z"]
    )
    profile_df.to_csv(TABLES_DIR / "local_correlations_by_target_and_section.tsv", sep="\t", index=False)

    summary_df = build_partner_summary(profile_df) if not profile_df.empty else pd.DataFrame()
    summary_df.to_csv(TABLES_DIR / "local_correlations_partner_summary.tsv", sep="\t", index=False)

    kcn_by_section = profile_df[profile_df["partner_gene"].isin(kcn_genes)].copy() if not profile_df.empty else pd.DataFrame()
    kcn_summary = summary_df[summary_df["partner_gene"].isin(kcn_genes)].copy() if not summary_df.empty else pd.DataFrame()
    kcn_by_section.to_csv(TABLES_DIR / "module_kcn_local_correlations_by_target_and_section.tsv", sep="\t", index=False)
    kcn_summary.to_csv(TABLES_DIR / "module_kcn_local_correlations_summary.tsv", sep="\t", index=False)

    section_summary_df = pd.DataFrame(section_summary_rows)
    section_summary_df.to_csv(TABLES_DIR / "target_gene_section_summary.tsv", sep="\t", index=False)

    if not profile_df.empty and not summary_df.empty:
        for target_gene in target_genes:
            target_section_df = profile_df.loc[profile_df["target_gene"] == target_gene].copy()
            target_summary_df = summary_df.loc[summary_df["target_gene"] == target_gene].copy()
            if target_section_df.empty or target_summary_df.empty:
                continue
            save_target_pdf(
                target_gene=target_gene,
                target_section_df=target_section_df,
                target_summary_df=target_summary_df,
                pdf_path=FIGURES_DIR / f"{target_gene.lower()}_local_correlations_primary_sections.pdf",
            )

    print("Module local correlations complete.")


if __name__ == "__main__":
    main()
