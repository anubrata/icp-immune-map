"""
Step 0 (atlas) - Turn the Human Immune Cell Atlas document into a marker table for step 3.

The atlas (resources/human_immune_cell_atlas.docx) is a Word document with one table of
marker genes per section ("7. Natural killer cells", ...). This script reads those tables
and writes them as a plain table that step 3 uses. Run it again whenever the document changes.

Input:  resources/human_immune_cell_atlas.docx
Output: resources/immune_cell_markers.tsv, one row per cell population or state:
          population_id  section number . row number, e.g. "7.1"
          section        e.g. "7. Natural killer cells"
          population     e.g. "General NK cells"
          kind           "cell population", or "functional state" for the cross-lineage section
                         (whose signatures, per the atlas, are not specific to immune cells)
          markers        marker genes, comma-separated (used for scoring)
          low_markers    genes the atlas lists as low/absent, e.g. "FCER1G-low" (not scored)
          source_text    the cell text exactly as written in the document

Run:  uv run analysis/step0_extract_atlas_markers.py
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import pandas as pd

from common import RESOURCES_DIR, show_path

ATLAS_DOCX = RESOURCES_DIR / "human_immune_cell_atlas.docx"
MARKERS_TSV = RESOURCES_DIR / "immune_cell_markers.tsv"
WORD = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
SYMBOL_FIXES = {"SPI-C": "SPIC"}  # spellings in the document that differ from the official symbol


def read_docx_blocks(path: Path):
    """Yield ("paragraph", text) and ("table", rows of cell texts) in document order."""
    with zipfile.ZipFile(path) as docx:
        body = ET.fromstring(docx.read("word/document.xml")).find(f"{WORD}body")

    def text(element) -> str:
        return "".join(t.text or "" for t in element.iter(f"{WORD}t")).strip()

    for element in body:
        if element.tag == f"{WORD}p" and text(element):
            yield "paragraph", text(element)
        elif element.tag == f"{WORD}tbl":
            yield "table", [[text(cell) for cell in row.findall(f"{WORD}tc")] for row in element.iter(f"{WORD}tr")]


def parse_markers(cell_text: str) -> tuple[list[str], list[str]]:
    """Split a marker cell into (markers, low markers).

    'KLRC2, FCER1G-low'                     -> (['KLRC2'], ['FCER1G'])
    'CD3D, CD3E, TRAC with low CD4/CD8A'    -> (['CD3D', 'CD3E', 'TRAC'], ['CD4', 'CD8A'])
    'CXCR4, SELL, ITGAM (context-dependent)' -> (['CXCR4', 'SELL', 'ITGAM'], [])
    """
    text = re.sub(r"\([^)]*\)", "", cell_text)  # drop notes in brackets
    listed, _, low_part = text.partition(" with low ")
    markers, low = [], [gene for gene in re.split(r"[/,\s]+", low_part) if gene]
    for gene in (g.strip() for g in listed.split(",")):
        if gene.endswith("-low"):
            low.append(gene.removesuffix("-low"))
        elif gene:
            markers.append(gene)
    return [SYMBOL_FIXES.get(g, g) for g in markers], [SYMBOL_FIXES.get(g, g) for g in low]


def extract(path: Path) -> pd.DataFrame:
    """One row per table row of the atlas, labelled with the numbered section heading above it."""
    rows, section = [], None
    for kind, content in read_docx_blocks(path):
        if kind == "paragraph" and re.match(r"\d+\.\s", content):
            section = content
        elif kind == "table" and section:
            number = section.split(".")[0]
            for row_number, cells in enumerate(content[1:], start=1):  # first row is the header
                population, genes = cells[0], cells[1]
                markers, low = parse_markers(genes)
                rows.append({
                    "population_id": f"{number}.{row_number}",
                    "section": section,
                    "population": population,
                    "kind": "functional state" if "across lineages" in section.lower() else "cell population",
                    "markers": ",".join(markers),
                    "low_markers": ",".join(low),
                    "source_text": genes,
                })
    return pd.DataFrame(rows)


def main() -> None:
    table = extract(ATLAS_DOCX)
    table.to_csv(MARKERS_TSV, sep="\t", index=False)
    print(f"{len(table)} populations/states in {table['section'].nunique()} sections -> {show_path(MARKERS_TSV)}")
    for section, rows in table.groupby("section", sort=False):
        print(f"  {section}: {len(rows)}")


if __name__ == "__main__":
    main()
