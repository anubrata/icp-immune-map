import zipfile

import pandas as pd

from common import RESOURCES_DIR
from step0_extract_atlas_markers import ATLAS_FILE, extract, parse_markers

NAMESPACE = 'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"'


def paragraph(text):
    return f"<w:p><w:r><w:t>{text}</w:t></w:r></w:p>"


def table(rows):
    return "<w:tbl>" + "".join("<w:tr>" + "".join(f"<w:tc>{paragraph(cell)}</w:tc>" for cell in row) + "</w:tr>"
                               for row in rows) + "</w:tbl>"


def test_extracts_each_table_row_under_its_section(tmp_path):
    body = (paragraph("Human Immune Cell Atlas")
            + paragraph("1. Pan-immune and major lineage markers")
            + table([["Cell population", "Representative genes"], ["NK cells", "NKG7, GNLY, KLRD1"],
                     ["Red-pulp macrophages", "CD163, SPI-C"]])
            + paragraph("13. Common immune-cell functional states across lineages")
            + table([["State", "Representative genes"], ["Cytotoxicity", "PRF1, GZMB"]])
            + paragraph("15. Recommended annotation procedure") + paragraph("Identify major cell lineages."))
    docx = tmp_path / "atlas.docx"
    with zipfile.ZipFile(docx, "w") as archive:
        archive.writestr("word/document.xml", f"<w:document {NAMESPACE}><w:body>{body}</w:body></w:document>")

    markers = extract(docx)
    assert markers["population_id"].tolist() == ["1.1", "1.2", "13.1"]
    assert markers["markers"].tolist() == ["NKG7,GNLY,KLRD1", "CD163,SPIC", "PRF1,GZMB"]  # SPI-C -> SPIC
    assert markers["kind"].tolist() == ["cell population", "cell population", "functional state"]
    assert markers.loc[2, "section"] == "13. Common immune-cell functional states across lineages"


def test_reads_the_plain_text_export(tmp_path):
    """Word's plain-text export: UTF-16, one table cell per line, wrapped paragraphs, invisible marks."""
    text = "\r\n".join([
        "Human Immune Cell Atlas: Gene Markers", "for Cell Types and Functional States",
        "‎1. Pan-immune and major lineage markers", "Cell population", "Representative genes",
        "NK cells", "NKG7, GNLY, KLRD1", "Double-negative αβ T cells", "CD3D, CD3E, TRAC with low CD4/CD8A/CD8B",
        "Note: these are candidate panels and several", "categories overlap.",
        "13. Common immune-cell functional states across lineages", "State", "Representative genes",
        "Cytotoxicity", "PRF1, GZMB", "15. Recommended annotation procedure",
        "‎1.‎\tIdentify major cell lineages using multiple gene markers.", ""])
    path = tmp_path / "atlas.txt"
    path.write_bytes(text.encode("utf-16"))

    markers = extract(path)
    assert markers["population_id"].tolist() == ["1.1", "1.2", "13.1"]
    assert markers["population"].tolist() == ["NK cells", "Double-negative αβ T cells", "Cytotoxicity"]
    assert markers.loc[1, "markers"] == "CD3D,CD3E,TRAC" and markers.loc[1, "low_markers"] == "CD4,CD8A,CD8B"
    assert markers["kind"].tolist() == ["cell population", "cell population", "functional state"]


def test_low_markers_and_notes_are_separated():
    assert parse_markers("KLRC2, FCER1G-low, ZBTB16-low") == (["KLRC2"], ["FCER1G", "ZBTB16"])
    assert parse_markers("CD3D, CD3E, TRAC with low CD4/CD8A/CD8B") == (["CD3D", "CD3E", "TRAC"], ["CD4", "CD8A", "CD8B"])
    assert parse_markers("CXCR4, SELL, ITGAM (context-dependent)") == (["CXCR4", "SELL", "ITGAM"], [])


def test_committed_marker_table_is_up_to_date_with_the_atlas():
    saved = pd.read_csv(RESOURCES_DIR / "immune_cell_markers.tsv", sep="\t", dtype=str, keep_default_na=False)
    fresh = extract(ATLAS_FILE)
    pd.testing.assert_frame_equal(saved, fresh)
