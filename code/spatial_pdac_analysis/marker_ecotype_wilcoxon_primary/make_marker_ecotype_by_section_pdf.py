from pathlib import Path

import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


BASE_DIR = Path(
    r"C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/code/spatial_pdac_analysis/marker_ecotype_wilcoxon_primary"
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


by_section_df = pd.read_csv(TABLE_DIR / "marker_ecotype_vs_rest_by_section.tsv", sep="\t")

report_df = by_section_df[
    [
        "orig.ident",
        "gene",
        "contrast",
        "n_inside",
        "n_outside",
        "log2fc_mean_expression",
        "auc",
        "pct_inside_positive",
        "pct_outside_positive",
        "delta_pct_positive",
        "wilcoxon_p",
        "wilcoxon_fdr",
    ]
].copy()

report_df = report_df.rename(
    columns={
        "orig.ident": "Section",
        "gene": "Marker",
        "contrast": "Contrast",
        "n_inside": "nIn",
        "n_outside": "nOut",
        "log2fc_mean_expression": "log2FC",
        "auc": "AUC",
        "pct_inside_positive": "%Pos in",
        "pct_outside_positive": "%Pos out",
        "delta_pct_positive": "dPct",
        "wilcoxon_p": "P",
        "wilcoxon_fdr": "FDR",
    }
)

contrast_order = {"CC1_CC5": 0, "CC2_CC3": 1}
report_df["contrast_order"] = report_df["Contrast"].map(contrast_order)
report_df["abs_log2fc"] = report_df["log2FC"].abs()

display_df = report_df.copy()
for col in ["log2FC", "AUC", "%Pos in", "%Pos out", "dPct"]:
    display_df[col] = display_df[col].map(fmt_decimal)
for col in ["P", "FDR"]:
    display_df[col] = display_df[col].map(fmt_scientific)

pdf_path = REPORT_DIR / "marker_ecotype_by_section_primary.pdf"

rows_per_page = 20
legend_lines = [
    "Par marqueur: ligne 1 = CC1+CC5, ligne 2 = CC2+CC3.",
    "nIn / nOut: nombre de spots dans l'ecotype / hors ecotype pour la coupe.",
    "log2FC: enrichissement moyen dans l'ecotype; AUC > 0.5 = plus haut dans l'ecotype.",
    "%Pos in / %Pos out: pourcentage de spots positifs; dPct: difference des deux; P et FDR: Wilcoxon.",
]

with PdfPages(pdf_path) as pdf:
    for section_name in sorted(display_df["Section"].unique()):
        section_df = display_df[display_df["Section"] == section_name].copy().reset_index(drop=True)

        section_rank_df = (
            report_df[report_df["Section"] == section_name]
            .groupby("Marker", as_index=False)
            .agg(
                min_fdr=("FDR", "min"),
                max_abs_log2fc=("abs_log2fc", "max"),
            )
            .sort_values(
                ["min_fdr", "max_abs_log2fc", "Marker"],
                ascending=[True, False, True],
            )
            .reset_index(drop=True)
        )
        marker_order = {marker: idx for idx, marker in enumerate(section_rank_df["Marker"])}
        section_df["marker_order"] = section_df["Marker"].map(marker_order)
        section_df = section_df.sort_values(["marker_order", "contrast_order"]).reset_index(drop=True)

        previous_marker = None
        for idx in section_df.index:
            current_marker = section_df.at[idx, "Marker"]
            if current_marker == previous_marker:
                section_df.at[idx, "Marker"] = ""
            previous_marker = current_marker

        n_pages = math.ceil(len(section_df) / rows_per_page)

        for page_idx in range(n_pages):
            start = page_idx * rows_per_page
            end = min(start + rows_per_page, len(section_df))
            page_df = section_df.iloc[start:end].drop(columns=["Section", "contrast_order", "abs_log2fc", "marker_order"])

            fig = plt.figure(figsize=(16.5, 11.7))
            ax = fig.add_axes([0.03, 0.20, 0.94, 0.72])
            ax.axis("off")

            title = f"Marker ecotype summary by section: {section_name}"
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
            table.set_fontsize(9)
            table.scale(1, 1.40)

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

            legend_y = 0.12
            fig.text(0.03, legend_y + 0.05, "Legend", fontsize=11, fontweight="bold", ha="left")
            for i, line in enumerate(legend_lines):
                fig.text(0.03, legend_y - i * 0.023, line, fontsize=9, ha="left")

            pdf.savefig(fig)
            plt.close(fig)

print(f"Saved PDF to: {pdf_path}")
