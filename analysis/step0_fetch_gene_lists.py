"""
Step 0 - Download the reference gene lists used in step 3.

Run this rarely (to refresh the lists). Its outputs are committed to the repository, so the
website build never depends on these servers being up.

  resources/immport_gene_lists.gmt  ImmPort's immune gene lists: ~150 lists built from the Gene
                                    Ontology term "immune system process" (GO:0002376) and the
                                    Reactome pathway "Immune System" (R-HSA-168256), with the
                                    gene symbols in each list.
  resources/hgnc_symbols.tsv        Official (HGNC) gene symbols with their Ensembl IDs and
                                    previous symbols, so genes that were renamed still match.

Run:  uv run analysis/step0_fetch_gene_lists.py
"""

from __future__ import annotations

import io
import urllib.request
from datetime import date

import pandas as pd

from common import RESOURCES_DIR, show_path

IMMPORT_GMT_URL = "https://s3.immport.org/release/genelists/current/all_gene_lists.gmt"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
HGNC_COLUMNS = ["symbol", "ensembl_gene_id", "prev_symbol"]


def download(url: str) -> bytes:
    print(f"Downloading {url}")
    with urllib.request.urlopen(url, timeout=300) as response:
        return response.read()


def main() -> None:
    RESOURCES_DIR.mkdir(exist_ok=True)

    # ImmPort: one gene list per line -> name <tab> link <tab> symbol <tab> symbol ...
    lines = download(IMMPORT_GMT_URL).decode("utf-8").splitlines()
    lists = [line for line in lines if line.count("\t") >= 2]
    if len(lists) < 50:
        raise SystemExit(f"The ImmPort file looks wrong: only {len(lists)} gene lists found.")
    gmt_path = RESOURCES_DIR / "immport_gene_lists.gmt"
    gmt_path.write_text("\n".join(lists) + "\n")
    n_genes = len({symbol for line in lists for symbol in line.split("\t")[2:] if symbol})
    print(f"  {len(lists)} immune gene lists covering {n_genes:,} genes -> {show_path(gmt_path)}")

    # HGNC: keep only the columns needed to translate old symbols / Ensembl IDs to current symbols.
    hgnc = pd.read_csv(io.BytesIO(download(HGNC_URL)), sep="\t", usecols=HGNC_COLUMNS, dtype=str)
    hgnc = hgnc[HGNC_COLUMNS].fillna("").sort_values("symbol")
    hgnc_path = RESOURCES_DIR / "hgnc_symbols.tsv"
    hgnc.to_csv(hgnc_path, sep="\t", index=False)
    print(f"  {len(hgnc):,} HGNC genes -> {show_path(hgnc_path)}")

    readme = f"""# Reference files

The ImmPort and HGNC files were downloaded on {date.today().isoformat()} by `analysis/step0_fetch_gene_lists.py`;
re-run it to refresh them. The atlas document is kept by hand: after changing it, run
`analysis/step0_extract_atlas_markers.py` to rebuild the marker table.

| File | Source | Contents |
| --- | --- | --- |
| `immport_gene_lists.gmt` | [ImmPort gene lists]({IMMPORT_GMT_URL}) ([docs](https://docs.immport.org/apidocumentation/genelist/overview/)) | {len(lists)} immune gene lists ({n_genes:,} genes) curated from GO "immune system process" (GO:0002376) and Reactome "Immune System" (R-HSA-168256). One list per line: name, link, gene symbols. |
| `hgnc_symbols.tsv` | [HGNC complete set]({HGNC_URL}) ([download page](https://www.genenames.org/download/archive/)) | Current official symbol, Ensembl gene ID and previous symbols for {len(hgnc):,} human genes. Used to match genes whose symbol changed. |
| `human_immune_cell_atlas.txt` | Human Immune Cell Atlas (project document, plain-text export) | Candidate marker genes for ~200 human immune cell populations and functional states, in 14 sections. |
| `immune_cell_markers.tsv` | Made from the atlas by `analysis/step0_extract_atlas_markers.py` | The atlas tables, one row per population: section, name, kind, markers, low markers, original text. |
"""
    (RESOURCES_DIR / "README.md").write_text(readme)


if __name__ == "__main__":
    main()
