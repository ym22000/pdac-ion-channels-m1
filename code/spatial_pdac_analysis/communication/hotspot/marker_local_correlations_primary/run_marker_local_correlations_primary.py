#!/usr/bin/env python
"""
Hotspot local correlation analysis for selected inflammatory-border stromal markers.

This run focuses on a small marker set used to describe resident/perivascular
or iCAF-like spatial territories in PDAC:
    - C7
    - C1R
    - CCL19
    - CCL21
    - MFAP4
    - CLU
    - DCN

The goal is to identify which genes, and especially which KCN genes, tend to
follow the same local spatial pattern as each marker across primary PDAC
sections. The workflow matches the main Hotspot local-correlation analysis, but
it writes tables only and does not generate PDFs.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import warnings

import anndata as ad
import hotspot
import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[2]
EXPORTS_DIR_DEFAULT = PROJECT_DIR / "communication" / "commot" / "shared_exports" / "exports"
KCN_LIST_PATH_DEFAULT = PROJECT_DIR.parent / "mycaf_icaf_correlations" / "kcn_union_25.tsv"
OUT_DIR_DEFAULT = SCRIPT_DIR

MODEL = "bernoulli"
N_NEIGHBORS = 30
MIN_GENE_SPOTS = 10
AUTOCORR_FDR_THRESHOLD = 0.05
MAX_SPATIAL_GENES = 400

TARGET_MARKERS = [
    "C7",
    "C1R",
    "CCL19",
    "CCL21",
    "MFAP4",
    "CLU",
    "DCN",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run section-wise Hotspot local correlations for selected marker genes."
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


def load_kcn_genes(kcn_list_path: Path) -> list[str]:
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


def keep_genes_for_hotspot(adata: ad.AnnData, forced_genes: list[str]) -> ad.AnnData:
    detected = np.asarray((adata.X > 0).sum(axis=0)).reshape(-1)
    keep = detected >= MIN_GENE_SPOTS

    for gene in forced_genes:
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


def select_genes_for_local_correlations(autocorr: pd.DataFrame, forced_genes: list[str]) -> list[str]:
    spatial_genes = autocorr.loc[
        (autocorr["FDR"] < AUTOCORR_FDR_THRESHOLD) & (autocorr["C"] > 0)
    ].sort_values("Z", ascending=False)

    genes = spatial_genes.index.tolist()[:MAX_SPATIAL_GENES]
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
    profile["target_gene"] = target_gene
    profile["orig.ident"] = image_name
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


def main() -> None:
    args = parse_args()
    exports_dir = Path(args.exports_dir)
    kcn_list_path = Path(args.kcn_list)
    out_dir = Path(args.out_dir)
    table_dir = out_dir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    kcn_genes = load_kcn_genes(kcn_list_path)
    target_genes = TARGET_MARKERS
    forced_genes = list(dict.fromkeys(target_genes + kcn_genes))

    section_dirs = sorted([p for p in exports_dir.iterdir() if p.is_dir()])
    if args.max_images is not None:
        section_dirs = section_dirs[: args.max_images]

    profile_rows = []
    target_meta_rows = []
    availability_rows = []

    for section_dir in section_dirs:
        image_name = section_dir.name
        print(f"[Hotspot marker local corr] Processing {image_name}...")
        try:
            adata = keep_genes_for_hotspot(load_section(section_dir), forced_genes)
            hs = build_hotspot(adata)
            autocorr = hs.compute_autocorrelations(jobs=1)
            genes_for_corr = select_genes_for_local_correlations(autocorr, forced_genes)
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

    kcn_set = set(kcn_genes)
    kcn_profile_df = profile_df.loc[profile_df["partner_gene"].isin(kcn_set)].copy()
    kcn_partner_summary_df = build_partner_summary(kcn_profile_df) if not kcn_profile_df.empty else pd.DataFrame(
        columns=[
            "target_gene",
            "partner_gene",
            "n_sections_tested",
            "n_sections_positive_z",
            "median_local_corr_z",
            "max_local_corr_z",
            "min_local_corr_z",
        ]
    )

    availability_df.to_csv(table_dir / "target_gene_availability.tsv", sep="\t", index=False)
    profile_df.to_csv(table_dir / "local_correlations_by_target_and_section.tsv", sep="\t", index=False)
    partner_summary_df.to_csv(table_dir / "local_correlations_partner_summary.tsv", sep="\t", index=False)
    target_meta_df.to_csv(table_dir / "target_gene_section_summary.tsv", sep="\t", index=False)
    kcn_profile_df.to_csv(table_dir / "marker_kcn_local_correlations_by_target_and_section.tsv", sep="\t", index=False)
    kcn_partner_summary_df.to_csv(table_dir / "marker_kcn_local_correlations_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
