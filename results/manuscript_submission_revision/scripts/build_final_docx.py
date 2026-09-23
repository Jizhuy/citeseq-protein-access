#!/usr/bin/env python3
"""Build the final presubmission DOCX from Markdown and inline figures."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Iterable

REPO = Path(__file__).resolve().parents[3]
REV = Path(__file__).resolve().parents[1]
MANUSCRIPT_ROOT = REPO.parent
MD = REV / "manuscript" / "CITEseq_Protein_Access_Model_Associated_Gains.md"
FIGS = REV / "figures"
TEMPLATE = MANUSCRIPT_ROOT / "CITEseq_Information_vs_Model_Reliability_Research_Manuscript_Final_Formatted.docx"
OUT_DOCX = REV / "manuscript" / "CITEseq_Protein_Access_Model_Gains_Submission_Ready.docx"
DEPS = REPO / ".deps_docx"
LOGS = REV / "logs"

from docx import Document  # noqa: E402
from docx.enum.section import WD_SECTION  # noqa: E402
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING  # noqa: E402
from docx.oxml import OxmlElement  # noqa: E402
from docx.oxml.ns import qn  # noqa: E402
from docx.shared import Cm, Inches, Pt, RGBColor  # noqa: E402

NAVY = RGBColor(0x1B, 0x3A, 0x5F)
TITLE = "Evaluating protein-access and residual totalVI-associated gains in CITE-seq representation learning under distribution shift"
AUTHOR_HEADER = "Jizhu Yang — Protein-access vs residual totalVI-associated contrast"
FALLBACK_FIGS = (
    MANUSCRIPT_ROOT
    / "multiomics_robustness"
    / "results"
    / "final_information_vs_model_manuscript"
    / "figures"
)
PARENT_FIGS = MANUSCRIPT_ROOT / "presubmission_revision" / "figures"
FIG_MAP = {
    1: "Figure1_central_logic.png",
    2: "Figure2_protein_access_vs_model.png",
    3: "Figure3_transfer_fairness.png",
    4: "Figure4_protein_degradation.png",
    5: "Figure5_uncertainty_stochasticity.png",
    6: "Figure6_lawlor_external.png",
}
FALLBACK_FIG_MAP = {
    1: "Figure1_central_logic.png",
    2: "Figure2_information_vs_model.png",
    3: "Figure3_transfer_fairness.png",
    4: "Figure4_protein_degradation.png",
    5: "Figure5_uncertainty_stochasticity.png",
    6: "Figure6_lawlor_external.png",
}


def normalize_text(text: str) -> str:
    """Convert the manuscript's lightweight Markdown/LaTeX into Word-friendly text."""
    replacements = {
        "\\mathrm{access}": "access",
        "\\mathrm{assoc}": "assoc",
        "\\mathrm{total}": "total",
        "\\mathrm{concat}": "concat",
        "\\mathrm{RNA\\,PCA}": "RNA_PCA",
        "\\mathrm{totalVI}": "totalVI",
        "\\mathrm{pAURC}_{[0.5,1.0]}": "pAURC_[0.5,1.0]",
        "\\approx": "approx.",
        "\\le": "<=",
        "−": "-",
        "–": "-",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    text = re.sub(r"\\\((.*?)\\\)", r"\1", text)
    text = re.sub(r"\\\[(.*?)\\\]", r"\1", text, flags=re.S)
    text = re.sub(r"([A-Za-z0-9]+)_\{([^{}]+)\}", r"\1_\2", text)
    text = text.replace("\\", "")
    return re.sub(r"\s+", " ", text).strip()


def set_run_font(run, size: float = 10, bold: bool = False, italic: bool = False, color=None) -> None:
    run.font.name = "Times New Roman"
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    if color is not None:
        run.font.color.rgb = color


def clear_document_body(doc) -> None:
    body = doc._body._element
    for child in list(body):
        if child.tag != qn("w:sectPr"):
            body.remove(child)


def add_page_number(paragraph) -> None:
    run = paragraph.add_run()
    fld_begin = OxmlElement("w:fldChar")
    fld_begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = " PAGE "
    fld_end = OxmlElement("w:fldChar")
    fld_end.set(qn("w:fldCharType"), "end")
    run._r.append(fld_begin)
    run._r.append(instr)
    run._r.append(fld_end)
    set_run_font(run, size=8, color=NAVY)


def configure_section(section, two_column: bool = False) -> None:
    section.page_width = Cm(21.0)
    section.page_height = Cm(29.7)
    section.left_margin = Cm(2.0)
    section.right_margin = Cm(2.0)
    section.top_margin = Cm(2.0)
    section.bottom_margin = Cm(2.0)
    section.header.is_linked_to_previous = False
    section.footer.is_linked_to_previous = False
    header = section.header.paragraphs[0]
    header.clear()
    header.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    hrun = header.add_run(AUTHOR_HEADER)
    set_run_font(hrun, size=8, italic=True, color=NAVY)
    footer = section.footer.paragraphs[0]
    footer.clear()
    footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
    add_page_number(footer)
    cols = section._sectPr.xpath("./w:cols")
    cols_el = cols[0] if cols else OxmlElement("w:cols")
    cols_el.set(qn("w:num"), "2" if two_column else "1")
    cols_el.set(qn("w:space"), "720")
    if not cols:
        section._sectPr.append(cols_el)


def configure_styles(doc) -> None:
    for style_name in ("Normal", "Body Text"):
        if style_name in doc.styles:
            style = doc.styles[style_name]
            style.font.name = "Times New Roman"
            style._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
            style.font.size = Pt(10)
    for style_name in ("Heading 1", "Heading 2", "Heading 3"):
        if style_name in doc.styles:
            style = doc.styles[style_name]
            style.font.name = "Times New Roman"
            style._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
            style.font.color.rgb = NAVY


def parse_md(text: str):
    text = re.sub(r"<!--.*?-->", "", text, flags=re.S)
    lines = text.splitlines()
    blocks = []
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if not stripped:
            i += 1
            continue
        if stripped.startswith("# "):
            blocks.append(("h1", stripped[2:].strip()))
        elif stripped.startswith("## "):
            blocks.append(("h2", stripped[3:].strip()))
        elif stripped.startswith("### "):
            blocks.append(("h3", stripped[4:].strip()))
        elif stripped.startswith("|") and i + 1 < len(lines) and re.match(r"^\|[\s\-:|]+\|$", lines[i + 1].strip()):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                current = lines[i].strip()
                if not re.match(r"^\|[\s\-:|]+\|$", current):
                    rows.append([normalize_text(c.strip()) for c in current.strip("|").split("|")])
                i += 1
            blocks.append(("table", rows))
            continue
        elif re.match(r"^\*\*Figure [1-6]\.", stripped) or re.match(r"^Figure [1-6]\.", stripped):
            blocks.append(("legend", stripped))
        else:
            para = [stripped]
            i += 1
            while i < len(lines):
                nxt = lines[i].strip()
                if (
                    not nxt
                    or nxt.startswith("#")
                    or nxt.startswith("|")
                    or nxt.startswith("**Figure ")
                    or re.match(r"^Figure [1-6]\.", nxt)
                ):
                    break
                para.append(nxt)
                i += 1
            blocks.append(("p", " ".join(para)))
            continue
        i += 1
    return blocks


def add_paragraph(
    doc,
    text: str,
    size: float = 10,
    bold: bool = False,
    italic: bool = False,
    align=None,
    space_after: float = 6,
    color=None,
):
    p = doc.add_paragraph()
    p.alignment = align if align is not None else WD_ALIGN_PARAGRAPH.JUSTIFY
    pf = p.paragraph_format
    pf.space_before = Pt(0)
    pf.space_after = Pt(space_after)
    pf.line_spacing_rule = WD_LINE_SPACING.SINGLE
    pf.first_line_indent = Pt(0)
    text = normalize_text(text)
    parts = re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*)", text)
    for part in parts:
        if not part:
            continue
        if part.startswith("**") and part.endswith("**"):
            run = p.add_run(part[2:-2])
            set_run_font(run, size=size, bold=True, color=NAVY)
        elif part.startswith("*") and part.endswith("*"):
            run = p.add_run(part[1:-1])
            set_run_font(run, size=size, italic=True)
        else:
            run = p.add_run(part)
            set_run_font(run, size=size, bold=bold, italic=italic, color=color)
    return p


def add_heading(doc, text: str, level: int) -> None:
    size = {1: 16, 2: 12, 3: 11}[level]
    align = WD_ALIGN_PARAGRAPH.CENTER if level == 1 else WD_ALIGN_PARAGRAPH.LEFT
    add_paragraph(doc, text, size=size, bold=True, align=align, space_after=8, color=NAVY)


def figure_path(n: int, caveats: list[str]) -> Path | None:
    primary = FIGS / FIG_MAP[n]
    if primary.exists():
        return primary
    fallback = FALLBACK_FIGS / FALLBACK_FIG_MAP[n]
    if fallback.exists():
        note = f"Figure {n}: used fallback source {fallback} because {primary} was unavailable."
        if note not in caveats:
            caveats.append(note)
        return fallback
    caveats.append(f"Figure {n}: missing {primary} and fallback {fallback}.")
    return None


def insert_figure(doc, n: int, caption: str | None, caveats: list[str]) -> bool:
    path = figure_path(n, caveats)
    if path is None:
        add_paragraph(doc, f"[Missing figure file: {FIG_MAP[n]}]", size=9, bold=True, color=NAVY)
        return False
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.first_line_indent = Pt(0)
    run = p.add_run()
    run.add_picture(str(path), width=Inches(6.35))
    if caption:
        cap = add_paragraph(
            doc,
            re.sub(r"^\*\*|\*\*$", "", caption),
            size=8,
            italic=True,
            align=WD_ALIGN_PARAGRAPH.JUSTIFY,
            space_after=10,
        )
        cap.paragraph_format.first_line_indent = Pt(0)
    return True


def add_table(doc, rows: list[list[str]]) -> None:
    if not rows:
        return
    width = max(len(r) for r in rows)
    table = doc.add_table(rows=len(rows), cols=width)
    table.style = "Table Grid"
    for ri, row in enumerate(rows):
        for ci in range(width):
            cell = table.cell(ri, ci)
            cell.text = row[ci] if ci < len(row) else ""
            for paragraph in cell.paragraphs:
                paragraph.paragraph_format.first_line_indent = Pt(0)
                paragraph.paragraph_format.space_after = Pt(0)
                for run in paragraph.runs:
                    set_run_font(run, size=7.3, bold=(ri == 0), color=NAVY if ri == 0 else None)
    doc.add_paragraph()


def collect_legends(blocks) -> dict[int, str]:
    legends: dict[int, str] = {}
    for kind, payload in blocks:
        if kind != "legend":
            continue
        match = re.match(r"^(?:\*\*)?Figure (\d)\.", payload)
        if match:
            n = int(match.group(1))
            clean = re.sub(r"^\*\*|\*\*$", "", payload)
            if n not in legends or len(clean) > len(legends[n]):
                legends[n] = clean
    return legends


def make_document(caveats: list[str]):
    if TEMPLATE.exists():
        doc = Document(str(TEMPLATE))
        clear_document_body(doc)
        caveats.append(f"Template formatting cloned from {TEMPLATE}.")
    else:
        doc = Document()
        caveats.append(f"Template unavailable at {TEMPLATE}; used professional single-column A4 fallback.")
    configure_styles(doc)
    for section in doc.sections:
        configure_section(section, two_column=False)
    caveats.append(
        "Two-column layout not applied: continuous two-column/body and full-width figure section breaks are fragile in python-docx, so the presubmission DOCX uses stable single-column A4 Times New Roman with navy headings."
    )
    return doc


def build_docx() -> tuple[Path, int, list[str]]:
    text = MD.read_text(encoding="utf-8")
    blocks = parse_md(text)
    legends = collect_legends(blocks)
    caveats: list[str] = []
    doc = make_document(caveats)
    inserted: set[int] = set()
    in_legends = False
    in_refs = False
    in_methods = False
    in_tables = False

    for kind, payload in blocks:
        if kind == "h1":
            add_heading(doc, payload, 1)
            continue
        if kind == "h2":
            heading = payload
            lower = heading.lower()
            in_legends = lower.startswith("figure legend")
            in_refs = lower.startswith("reference")
            in_methods = lower.startswith("methods")
            in_tables = lower.startswith("table")
            if in_legends:
                continue
            add_heading(doc, heading, 2)
            continue
        if kind == "h3":
            if not in_legends:
                add_heading(doc, payload, 3)
            continue
        if kind == "legend":
            # Figure callouts that begin a paragraph are parsed like legends; keep
            # them inline unless they are in the dedicated Figure Legends section.
            if in_legends:
                continue
            add_paragraph(doc, payload, size=10, space_after=6)
            match = re.match(r"^(?:\*\*)?Figure ([1-6])\.", payload)
            if match:
                n = int(match.group(1))
                if n not in inserted:
                    if insert_figure(doc, n, legends.get(n, payload), caveats):
                        inserted.add(n)
            continue
        if kind == "table":
            if not in_legends:
                add_table(doc, payload)
            continue
        if kind == "p":
            if in_legends:
                continue
            center_front_matter = (
                payload.startswith("**Jizhu Yang**")
                or payload.startswith("Research Manuscript")
                or payload.startswith("Unpublished")
                or payload.startswith("September")
            )
            if center_front_matter:
                add_paragraph(
                    doc,
                    payload.replace("**", ""),
                    size=11,
                    bold=payload.startswith("**Jizhu Yang**"),
                    align=WD_ALIGN_PARAGRAPH.CENTER,
                    space_after=2,
                )
                continue
            size = 9 if (in_methods or in_refs or in_tables) else 10
            p = add_paragraph(doc, payload, size=size, space_after=4 if in_refs else 6)
            if in_refs and re.match(r"^\d+\.\s", normalize_text(payload)):
                p.paragraph_format.left_indent = Cm(0.55)
                p.paragraph_format.first_line_indent = Cm(-0.55)
            match = re.search(r"\*\*Figure ([1-6])\.", payload)
            if match:
                n = int(match.group(1))
                if n not in inserted:
                    if insert_figure(doc, n, legends.get(n, payload), caveats):
                        inserted.add(n)

    missing_insertions = [n for n in range(1, 7) if n not in inserted]
    if missing_insertions:
        caveats.append(f"No end-of-document figure dump was added; missing inline insertions: {missing_insertions}.")

    LOGS.mkdir(parents=True, exist_ok=True)
    OUT_DOCX.parent.mkdir(parents=True, exist_ok=True)
    doc.save(OUT_DOCX)
    # Also publish concise filenames under revision manuscript/ and repo root if writable
    alt = REV / "manuscript" / "CITEseq_Protein_Access_Model_Gains.docx"
    if OUT_DOCX.resolve() != alt.resolve():
        alt.write_bytes(OUT_DOCX.read_bytes())
    print(f"wrote {OUT_DOCX}")
    print("formatting_caveat: stable single-column A4 layout used instead of fragile two-column section switching.")
    for caveat in caveats:
        print(f"caveat: {caveat}")
    return OUT_DOCX, abstract_word_count(text), caveats


def section_text(markdown: str, heading: str, stop_headings: Iterable[str] = ()) -> str:
    pattern = rf"\n## {re.escape(heading)}\n(?P<body>.*?)(?=\n## (?:{'|'.join(map(re.escape, stop_headings))})\n|\Z)"
    if stop_headings:
        match = re.search(pattern, markdown, flags=re.S)
    else:
        match = re.search(rf"\n## {re.escape(heading)}\n(?P<body>.*?)(?=\n## |\Z)", markdown, flags=re.S)
    return match.group("body") if match else ""


def abstract_word_count(markdown: str) -> int:
    abstract = section_text(markdown, "Abstract")
    abstract = re.sub(r"\*\*", "", normalize_text(abstract))
    return len(re.findall(r"\b[\w.-]+\b", abstract))


def all_docx_text(doc) -> str:
    chunks = [p.text for p in doc.paragraphs]
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                chunks.extend(p.text for p in cell.paragraphs)
    for section in doc.sections:
        chunks.extend(p.text for p in section.header.paragraphs)
        chunks.extend(p.text for p in section.footer.paragraphs)
    return "\n".join(chunks)


def validation_contexts(text: str, needle: str) -> list[str]:
    lines = [line.strip() for line in text.splitlines() if needle in line]
    return lines


def validate(docx_path: Path, abstract_wc: int, caveats: list[str]) -> dict:
    doc = Document(str(docx_path))
    text = all_docx_text(doc)
    body_text = "\n".join(p.text for p in doc.paragraphs)
    front_matter = body_text.split("Abstract", 1)[0] if "Abstract" in body_text else body_text
    lower = text.lower()
    ack_context = section_text(MD.read_text(encoding="utf-8"), "Acknowledgements")
    bo_contexts = validation_contexts(text, "Bo Li")
    checks = {
        "title_correct": TITLE in text,
        "sole_author_jizhu_yang": "Jizhu Yang" in front_matter and "Bo Li" not in front_matter,
        "bo_li_only_in_acknowledgements": bool(bo_contexts) and all("valuable discussions" in c for c in bo_contexts) and "Bo Li" in ack_context,
        "n_images": len(doc.inline_shapes),
        "six_images": len(doc.inline_shapes) == 6,
        "protein_gated_present": "protein-gated" in lower,
        "g_access_present": "G_access" in text,
        "lawlor_means_0_794_0_884": ("0.794" in text and "0.884" in text),
        "old_lawlor_0_800_0_888_absent": ("0.800" not in text and "0.888" not in text),
        "sign_test_p_00195": "0.00195" in text,
        "g_assoc_crosses_zero": (
            ("-0.060" in text or "−0.060" in text)
            and ("crosses zero" in lower or "crossed zero" in lower)
        ),
        "case_c_absent": "Case C" not in text,
        "application_absent": "application" not in lower,
        "bootstrap_interval_0_042_present": "0.042" in text,
        "bootstrap_interval_minus_0_060_present": ("-0.060" in text or "−0.060" in text),
        "strict_train_only_present": "strict train-only" in lower,
        "totalvi_primary_0_741": "0.741" in text,
        "github_url_present": ("https://github.com/Jizhuy/robust-multiomics-integration" in text),
        "data_availability_before_references": text.find("Data Availability") != -1
        and text.find("References") != -1
        and text.find("Data Availability") < text.find("References"),
        "code_availability_before_references": text.find("Code Availability") != -1
        and text.find("References") != -1
        and text.find("Code Availability") < text.find("References"),
        "acknowledgements_before_references": text.find("Acknowledgements") != -1
        and text.find("References") != -1
        and text.find("Acknowledgements") < text.find("References"),
        "tables_before_references": bool(doc.tables) and text.find("Table 1.") != -1 and text.find("Table 1.") < text.find("References"),
        "abstract_word_count": abstract_wc,
        "formatting_caveats": caveats,
    }
    LOGS.mkdir(parents=True, exist_ok=True)
    report = LOGS / "docx_validation.json"
    report.write_text(json.dumps(checks, indent=2, ensure_ascii=False), encoding="utf-8")
    print("validation", json.dumps(checks, ensure_ascii=False, sort_keys=True))
    print(f"wrote {report}")
    return checks


def main() -> None:
    out, abs_wc, caveats = build_docx()
    checks = validate(out, abs_wc, caveats)
    if not checks["six_images"]:
        raise SystemExit("validation failed: expected exactly 6 embedded images")
    required_bools = [
        "title_correct",
        "sole_author_jizhu_yang",
        "bo_li_only_in_acknowledgements",
        "protein_gated_present",
        "g_access_present",
        "case_c_absent",
        "application_absent",
        "bootstrap_interval_0_042_present",
        "bootstrap_interval_minus_0_060_present",
        "strict_train_only_present",
        "totalvi_primary_0_741",
        "github_url_present",
        "lawlor_means_0_794_0_884",
        "old_lawlor_0_800_0_888_absent",
        "sign_test_p_00195",
        "g_assoc_crosses_zero",
    ]
    failed = [key for key in required_bools if not checks.get(key)]
    if failed:
        raise SystemExit(f"validation failed: {failed}")


if __name__ == "__main__":
    main()
