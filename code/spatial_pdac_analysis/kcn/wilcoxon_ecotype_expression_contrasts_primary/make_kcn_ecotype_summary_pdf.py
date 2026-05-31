from pathlib import Path

import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


BASE_DIR = Path(
    r"C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/code/spatial_pdac_analysis/kcn/wilcoxon_ecotype_expression_contrasts_primary"
)
TABLE_DIR = BASE_DIR / "tables"
REPORT_DIR = BASE_DIR / "reports"
REPORT_DIR.mkdir(parents=True, exist_ok=True)


def fmt_decimal(value, digits=3):
    if pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}f}"


def fmt_scientific(value):
    if pd.isna(value):
        return "NA"
    return f"{float(value):.2e}"


summary_df = pd.read_csv(TABLE_DIR / "kcn_ecotype_vs_rest_summary.tsv", sep="\t")
by_section_df = pd.read_csv(TABLE_DIR / "kcn_ecotype_vs_rest_by_section.tsv", sep="\t")

section_stats = (
    by_section_df.groupby(["gene", "contrast"], as_index=False)
    .agg(
        median_p=("wilcoxon_p", "median"),
        median_fdr=("wilcoxon_fdr", "median"),
        best_fdr=("wilcoxon_fdr", "min"),
    )
)

report_df = summary_df.merge(section_stats, on=["gene", "contrast"], how="left")

report_df = report_df[
    [
        "gene",
        "contrast",
        "n_sections_tested",
        "n_sections_fdr_lt_0_05",
        "median_log2fc_mean_expression",
        "median_auc",
        "median_pct_inside_positive",
        "median_pct_outside_positive",
        "median_delta_pct_positive",
        "median_p",
        "median_fdr",
        "best_fdr",
    ]
].copy()

report_df = report_df.rename(
    columns={
        "gene": "KCN",
        "contrast": "Contrast",
        "n_sections_tested": "nSec",
        "n_sections_fdr_lt_0_05": "nFDR<0.05",
        "median_log2fc_mean_expression": "medLog2FC",
        "median_auc": "medAUC",
        "median_pct_inside_positive": "%Pos in",
        "median_pct_outside_positive": "%Pos out",
        "median_delta_pct_positive": "dPct",
        "median_p": "medP",
        "median_fdr": "medFDR",
        "best_fdr": "bestFDR",
    }
)

display_df = report_df.copy()
for col in ["medLog2FC", "medAUC", "%Pos in", "%Pos out", "dPct"]:
    display_df[col] = display_df[col].map(fmt_decimal)
for col in ["medP", "medFDR", "bestFDR"]:
    display_df[col] = display_df[col].map(fmt_scientific)

contrast_order = {"CC1_CC5": 0, "CC2_CC3": 1}
display_df["contrast_order"] = display_df["Contrast"].map(contrast_order)

# Rank KCN globally by strength of FDR support across both contrasts.
kcn_rank_df = report_df.copy()
kcn_rank_df["abs_log2fc"] = kcn_rank_df["medLog2FC"].abs()
kcn_rank_df = (
    kcn_rank_df.groupby("KCN", as_index=False)
    .agg(
        max_n_fdr=("nFDR<0.05", "max"),
        min_best_fdr=("bestFDR", "min"),
        max_abs_log2fc=("abs_log2fc", "max"),
    )
    .sort_values(
        ["max_n_fdr", "min_best_fdr", "max_abs_log2fc", "KCN"],
        ascending=[False, True, False, True],
    )
    .reset_index(drop=True)
)
kcn_order = {kcn: idx for idx, kcn in enumerate(kcn_rank_df["KCN"])}
display_df["kcn_order"] = display_df["KCN"].map(kcn_order)
display_df = display_df.sort_values(["kcn_order", "contrast_order"]).reset_index(drop=True)
display_df = display_df.drop(columns=["kcn_order"])
display_df = display_df.drop(columns=["contrast_order"])

# Keep the KCN label only on the first of the two rows to make the grouping clearer.
previous_kcn = None
for idx in display_df.index:
    current_kcn = display_df.at[idx, "KCN"]
    if current_kcn == previous_kcn:
        display_df.at[idx, "KCN"] = ""
    previous_kcn = current_kcn

pdf_path = REPORT_DIR / "kcn_ecotype_summary_primary.pdf"

rows_per_page = 24
n_pages = math.ceil(len(display_df) / rows_per_page)

legend_lines = [
    "Chaque KCN est montre sur deux lignes: d'abord CC1+CC5, puis CC2+CC3.",
    "nSec: nombre de coupes testees; nFDR<0.05: coupes significatives.",
    "medLog2FC: mediane du log2FC (ecotype vs reste); medAUC: mediane de l'AUC (>0.5 = plus haut dans l'ecotype).",
    "%Pos in / %Pos out: pourcentage median de spots positifs dans l'ecotype / hors ecotype; dPct: difference des deux.",
    "medP / medFDR: medianes des p-values et FDR section par section; bestFDR: meilleure FDR observee.",
]

with PdfPages(pdf_path) as pdf:
    for page_idx in range(n_pages):
        start = page_idx * rows_per_page
        end = min(start + rows_per_page, len(display_df))
        page_df = display_df.iloc[start:end]

        fig = plt.figure(figsize=(16.5, 11.7))
        ax = fig.add_axes([0.03, 0.22, 0.94, 0.70])
        ax.axis("off")

        title = "KCN ecotype-vs-rest summary in primary PDAC (all KCN, both contrasts)"
        if n_pages > 1:
            title = f"{title} (page {page_idx + 1}/{n_pages})"
        fig.text(0.03, 0.95, title, fontsize=18, fontweight="bold", ha="left")

        table = ax.table(
            cellText=page_df.values,
            colLabels=page_df.columns,
            cellLoc="center",
            loc="upper left",
            bbox=[0, 0, 1, 1],
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8.7)
        table.scale(1, 1.45)

        for (row, col), cell in table.get_celld().items():
            if row == 0:
                cell.set_text_props(weight="bold", color="black")
                cell.set_facecolor("#d9ead3")
            else:
                data_idx = start + row - 1
                if data_idx % 2 == 0:
                    cell.set_facecolor("#f7f7f7")
                else:
                    cell.set_facecolor("#ffffff")

        legend_y = 0.14
        fig.text(0.03, legend_y + 0.04, "Legend", fontsize=11, fontweight="bold", ha="left")
        for i, line in enumerate(legend_lines):
            fig.text(0.03, legend_y - i * 0.025, line, fontsize=9, ha="left")

        pdf.savefig(fig)
        plt.close(fig)

print(f"Saved PDF to: {pdf_path}")
