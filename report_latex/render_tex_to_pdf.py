from __future__ import annotations

import re
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer


ROOT = Path(__file__).resolve().parent
TEX_PATH = ROOT / "main.tex"
OUT_PATH = ROOT / "report_skeleton_render.pdf"
BASE_FONT = "Times-Roman"


def register_fonts() -> None:
    global BASE_FONT
    font_dir = Path("C:/Windows/Fonts")
    regular = font_dir / "TCM_____.TTF"
    italic = font_dir / "TCMI____.TTF"
    bold = font_dir / "TCCM____.TTF"

    if regular.exists():
        pdfmetrics.registerFont(TTFont("ClassicTex", str(regular)))
        if italic.exists():
            pdfmetrics.registerFont(TTFont("ClassicTex-Italic", str(italic)))
        if bold.exists():
            pdfmetrics.registerFont(TTFont("ClassicTex-Bold", str(bold)))

        pdfmetrics.registerFontFamily(
            "ClassicTex",
            normal="ClassicTex",
            bold="ClassicTex-Bold" if bold.exists() else "ClassicTex",
            italic="ClassicTex-Italic" if italic.exists() else "ClassicTex",
            boldItalic="ClassicTex-Bold" if bold.exists() else "ClassicTex",
        )
        BASE_FONT = "ClassicTex"


def build_styles():
    styles = getSampleStyleSheet()

    styles.add(
        ParagraphStyle(
            name="CoverTop",
            parent=styles["Normal"],
            fontName=BASE_FONT,
            fontSize=18,
            leading=24,
            alignment=TA_CENTER,
            spaceAfter=8,
        )
    )
    styles.add(
        ParagraphStyle(
            name="CoverSub",
            parent=styles["Normal"],
            fontName=BASE_FONT,
            fontSize=14,
            leading=18,
            alignment=TA_CENTER,
            spaceAfter=8,
        )
    )
    styles.add(
        ParagraphStyle(
            name="CoverTitle",
            parent=styles["Normal"],
            fontName=BASE_FONT,
            fontSize=20,
            leading=26,
            alignment=TA_CENTER,
            spaceAfter=10,
        )
    )
    styles.add(
        ParagraphStyle(
            name="Body",
            parent=styles["Normal"],
            fontName=BASE_FONT,
            fontSize=10.5,
            leading=15,
            spaceAfter=8,
        )
    )
    styles.add(
        ParagraphStyle(
            name="SectionCustom",
            parent=styles["Heading1"],
            fontName=BASE_FONT,
            fontSize=15,
            leading=20,
            textColor=colors.black,
            spaceBefore=12,
            spaceAfter=8,
        )
    )
    styles.add(
        ParagraphStyle(
            name="SubsectionCustom",
            parent=styles["Heading2"],
            fontName=BASE_FONT,
            fontSize=12.5,
            leading=17,
            textColor=colors.black,
            spaceBefore=10,
            spaceAfter=6,
        )
    )
    styles.add(
        ParagraphStyle(
            name="CaptionCustom",
            parent=styles["Normal"],
            fontName=BASE_FONT,
            fontSize=9.5,
            leading=13,
            spaceBefore=4,
            spaceAfter=10,
        )
    )
    return styles


def tex_to_html(text: str) -> str:
    text = text.strip()
    text = text.replace(r"\textwidth", "")
    text = re.sub(r"\\textit\{([^{}]+)\}", r"<i>\1</i>", text)
    text = re.sub(r"\\textbf\{([^{}]+)\}", r"<b>\1</b>", text)
    text = re.sub(r"\\texttt\{([^{}]+)\}", r"<font face='Courier'>\1</font>", text)
    text = text.replace(r"\textit", "")
    text = text.replace(r"\par", "")
    text = text.replace(r"\newpage", "")
    text = text.replace(r"\textasciitilde{}", "~")
    text = text.replace(r"\%", "%")
    text = text.replace(r"\&", "&")
    text = text.replace(r"---", "—")
    text = text.replace(r"--", "–")
    text = text.replace(r"$", "")
    text = escape(text, {'"': "&quot;"})
    text = text.replace("&lt;i&gt;", "<i>").replace("&lt;/i&gt;", "</i>")
    text = text.replace("&lt;b&gt;", "<b>").replace("&lt;/b&gt;", "</b>")
    text = text.replace("&lt;font face='Courier'&gt;", "<font face='Courier'>")
    text = text.replace("&lt;/font&gt;", "</font>")
    text = re.sub(r"(\[(?:\d+(?:,\s*)?)+\])", r"<link href='#references'>\1</link>", text)
    return text


def extract_braced_text(line: str) -> str:
    line = line.strip()
    if not line.startswith("{"):
        return line
    match = re.search(r"\{.*? (.+?)\\par\}", line)
    if match:
        return match.group(1).strip()
    match = re.search(r"\{.*? (.+?)\}", line)
    if match:
        return match.group(1).strip()
    return line


def parse_tex():
    text = TEX_PATH.read_text(encoding="utf-8")
    lines = text.splitlines()
    elements: list[tuple[str, dict]] = []

    in_title = False
    in_figure = False
    figure_data = {}
    paragraph_buffer: list[str] = []

    for raw_line in lines:
        line = raw_line.strip()
        if not line or line.startswith("%"):
            if paragraph_buffer:
                elements.append(("paragraph", {"text": " ".join(paragraph_buffer)}))
                paragraph_buffer = []
            continue

        if line.startswith(r"\begin{titlepage}"):
            in_title = True
            continue
        if line.startswith(r"\end{titlepage}"):
            in_title = False
            elements.append(("titlepage_end", {}))
            continue

        if in_title:
            if line.startswith(r"\centering") or line.startswith(r"\vfill") or line.startswith(r"\vspace"):
                continue
            elements.append(("cover_line", {"text": extract_braced_text(line)}))
            continue

        if line.startswith(r"\tableofcontents"):
            elements.append(("toc", {}))
            continue
        if line.startswith(r"\newpage"):
            elements.append(("pagebreak", {}))
            continue

        if line.startswith(r"\begin{figure}"):
            in_figure = True
            figure_data = {}
            continue
        if line.startswith(r"\end{figure}"):
            in_figure = False
            elements.append(("figure", figure_data.copy()))
            figure_data = {}
            continue
        if in_figure:
            inc = re.search(r"\\includegraphics\[.*?\]\{(.+?)\}", line)
            cap = re.search(r"\\caption\{(.+)\}", line)
            if inc:
                figure_data["image"] = inc.group(1)
            if cap:
                figure_data["caption"] = cap.group(1)
            continue

        sec = re.search(r"\\section\*?\{(.+?)\}", line)
        sub = re.search(r"\\subsection\{(.+?)\}", line)
        if sec:
            if paragraph_buffer:
                elements.append(("paragraph", {"text": " ".join(paragraph_buffer)}))
                paragraph_buffer = []
            elements.append(("section", {"text": sec.group(1)}))
            continue
        if sub:
            if paragraph_buffer:
                elements.append(("paragraph", {"text": " ".join(paragraph_buffer)}))
                paragraph_buffer = []
            elements.append(("subsection", {"text": sub.group(1)}))
            continue

        if line.startswith("\\documentclass") or line.startswith("\\usepackage") or line.startswith("\\graphicspath") or line.startswith("\\setlength") or line.startswith("\\onehalfspacing") or line.startswith("\\titleformat") or line.startswith("\\begin{document}") or line.startswith("\\end{document}"):
            continue

        paragraph_buffer.append(line)

    if paragraph_buffer:
        elements.append(("paragraph", {"text": " ".join(paragraph_buffer)}))

    return elements


def build_story(elements):
    styles = build_styles()
    story = []

    cover_index = 0
    cover_styles = [
        "CoverTop",
        "CoverSub",
        "CoverSub",
        "CoverTitle",
        "CoverSub",
        "Body",
        "Body",
        "Body",
        "Body",
        "Body",
    ]

    for kind, payload in elements:
        if kind == "cover_line":
            text = tex_to_html(payload["text"])
            style_name = cover_styles[min(cover_index, len(cover_styles) - 1)]
            story.append(Spacer(1, 0.15 * cm))
            story.append(Paragraph(text, styles[style_name]))
            cover_index += 1
        elif kind == "titlepage_end":
            story.append(PageBreak())
        elif kind == "toc":
            story.append(Paragraph("Table of Contents", styles["SectionCustom"]))
            toc_items = [
                "Abstract",
                "List of Abbreviations",
                "1. Introduction",
                "2. Objectives",
                "3. Materials and Methods",
                "4. Results and Discussion",
                "5. Conclusions and Perspectives",
                "Clinical Opening and Limits",
                "References",
                "M1 Skills Summary",
                "Seminar Attendance Sheet",
            ]
            for item in toc_items:
                story.append(Paragraph(tex_to_html(item), styles["Body"]))
            story.append(PageBreak())
        elif kind == "pagebreak":
            story.append(PageBreak())
        elif kind == "section":
            heading = payload["text"]
            if heading == "References":
                story.append(Paragraph("<a name='references'/>" + tex_to_html(heading), styles["SectionCustom"]))
            else:
                story.append(Paragraph(tex_to_html(heading), styles["SectionCustom"]))
        elif kind == "subsection":
            story.append(Paragraph(tex_to_html(payload["text"]), styles["SubsectionCustom"]))
        elif kind == "paragraph":
            story.append(Paragraph(tex_to_html(payload["text"]), styles["Body"]))
        elif kind == "figure":
            img_name = payload.get("image")
            caption = payload.get("caption", "")
            if img_name:
                img_path = ROOT.parent / "report_figures" / img_name
                if img_path.exists():
                    img = Image(str(img_path))
                    max_width = 16 * cm
                    max_height = 21 * cm
                    scale = min(max_width / img.drawWidth, max_height / img.drawHeight)
                    img.drawWidth *= scale
                    img.drawHeight *= scale
                    story.append(img)
            if caption:
                story.append(Paragraph(tex_to_html(caption), styles["CaptionCustom"]))

    return story


def main() -> None:
    register_fonts()
    elements = parse_tex()
    story = build_story(elements)
    doc = SimpleDocTemplate(
        str(OUT_PATH),
        pagesize=A4,
        leftMargin=2.2 * cm,
        rightMargin=2.2 * cm,
        topMargin=2.0 * cm,
        bottomMargin=2.0 * cm,
    )
    doc.build(story)
    print(f"PDF written to {OUT_PATH}")


if __name__ == "__main__":
    main()
