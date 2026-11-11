"""
Step 0 (mouse datasets) - Mouse gene annotation and human orthologs from MGI.

Mouse data files often name genes only by their Ensembl ID (ENSMUSG...). This script builds one
table that gives every mouse gene its official symbol, name, chromosome and type, and its human
ortholog, from two reports of MGI (Mouse Genome Informatics, the official mouse gene database):

  MRK_ENSEMBL.rpt             every mouse gene (MGI ID) with its Ensembl gene ID
  HOM_MouseHumanSequence.rpt  mouse-human homology classes (from the Alliance of Genome Resources):
                              each class groups the mouse and human genes descended from one ancestral gene

Why orthologs: the Human Immune Cell Atlas lists human genes. A mouse gene can stand in for an
atlas marker only through its human ortholog, never because the names look alike (mouse Gnly
does not exist; mouse Cd8b1 is human CD8B).

One-to-one only: a mouse gene gets a `human_ortholog` only when its homology class holds exactly
one mouse and one human gene. Classes with several genes on either side (e.g. the Ccl3/Ccl4 cluster,
Fcgr genes, Ly49/KIR receptors) are ambiguous, so those genes get no ortholog and are not matched
to atlas markers; their human genes are kept in `human_genes` for reference.

Output: resources/mouse_genes.tsv.gz, one row per mouse Ensembl gene ID:
          gene_id          Ensembl gene ID (ENSMUSG...)
          gene_symbol      official MGI symbol (e.g. Nkg7)
          description      gene name
          chromosome       e.g. "7", "X", "MT"
          biotype          e.g. "protein coding", "lncRNA"
          mgi_id           MGI accession (e.g. MGI:1931250)
          human_ortholog   the human gene symbol when the ortholog is one-to-one, else ""
          ortholog_type    one-to-one, one-to-many, many-to-one, many-to-many, or none
          human_genes      all human genes in the homology class (';'-separated)

Run:  uv run analysis/step0_fetch_mouse_genes.py
      (--from-dir <folder> uses already downloaded copies of the two reports)
"""

from __future__ import annotations

import argparse
import http.client
import tempfile
import urllib.request
from pathlib import Path

import pandas as pd

from common import RESOURCES_DIR, show_path

MGI_URL = "https://www.informatics.jax.org/downloads/reports/{name}"
GENE_REPORT, HOMOLOGY_REPORT = "MRK_ENSEMBL.rpt", "HOM_MouseHumanSequence.rpt"
OUTPUT = RESOURCES_DIR / "mouse_genes.tsv.gz"

# MRK_ENSEMBL.rpt has no header line; its columns, in order:
GENE_REPORT_COLUMNS = ["mgi_id", "gene_symbol", "description", "cm_position", "chromosome", "gene_id",
                       "ensembl_transcripts", "ensembl_proteins", "feature_type", "start", "end", "strand", "biotypes"]
MOUSE_TAXON, HUMAN_TAXON = "10090", "9606"
OUTPUT_COLUMNS = ["gene_id", "gene_symbol", "description", "chromosome", "biotype", "mgi_id",
                  "human_ortholog", "ortholog_type", "human_genes"]


def download(name: str, folder: Path, attempts: int = 3) -> Path:
    """Download one MGI report into `folder` (with a few retries)."""
    url, path = MGI_URL.format(name=name), folder / name
    for attempt in range(1, attempts + 1):
        print(f"Downloading {url} (attempt {attempt})")
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "icp-immune-map (python urllib)"})
            with urllib.request.urlopen(request, timeout=600) as response:
                path.write_bytes(response.read())
            return path
        except (OSError, http.client.HTTPException) as error:
            print(f"  failed: {error!r}")
    raise SystemExit(f"Could not download {url}. Download it another way and use --from-dir.")


def read_gene_report(path: Path) -> pd.DataFrame:
    """MRK_ENSEMBL.rpt -> one row per Ensembl gene ID with symbol, name, chromosome and type."""
    genes = pd.read_csv(path, sep="\t", header=None, names=GENE_REPORT_COLUMNS, dtype=str, keep_default_na=False)
    if genes["gene_id"].duplicated().any():
        raise SystemExit(f"{path.name}: an Ensembl gene ID appears twice; the report format may have changed")
    return genes


def ortholog_table(homology: pd.DataFrame) -> pd.DataFrame:
    """One row per mouse gene (MGI ID) in HOM_MouseHumanSequence.rpt: its human ortholog (one-to-one
    only), the type of the relationship, and every human gene it is paired with.

    Pairs are counted over the whole report, because one human gene can sit in several classes
    with a different mouse gene in each. One-to-one means: this mouse gene has exactly one human
    partner, and that human gene has no other mouse partner."""
    mouse = homology[homology["NCBI Taxon ID"] == MOUSE_TAXON][["DB Class Key", "Mouse MGI ID"]]
    human = homology[homology["NCBI Taxon ID"] == HUMAN_TAXON][["DB Class Key", "Symbol"]]
    pairs = mouse.merge(human, on="DB Class Key").drop_duplicates(["Mouse MGI ID", "Symbol"])
    humans_of = pairs.groupby("Mouse MGI ID")["Symbol"].agg(lambda s: sorted(set(s)))
    mice_of = pairs.groupby("Symbol")["Mouse MGI ID"].nunique()

    rows = []
    for mgi_id in mouse["Mouse MGI ID"].unique():
        humans = humans_of.get(mgi_id, [])
        shared = any(mice_of[h] > 1 for h in humans)  # a partner that is also paired with another mouse gene
        if not humans:
            kind = "none"
        elif len(humans) == 1:
            kind = "many-to-one" if shared else "one-to-one"
        else:
            kind = "many-to-many" if shared else "one-to-many"
        rows.append({"mgi_id": mgi_id, "human_ortholog": humans[0] if kind == "one-to-one" else "",
                     "ortholog_type": kind, "human_genes": ";".join(humans)})
    return pd.DataFrame(rows).set_index("mgi_id")


def build(genes: pd.DataFrame, homology: pd.DataFrame) -> pd.DataFrame:
    """The output table: every mouse Ensembl gene with its annotation and human ortholog."""
    orthologs = ortholog_table(homology)
    table = genes.join(orthologs, on="mgi_id")
    # "protein coding gene" -> "protein coding", "lncRNA gene" -> "lncRNA"
    table["biotype"] = table["feature_type"].str.replace(r"\s+gene$", "", regex=True)
    table["ortholog_type"] = table["ortholog_type"].fillna("none")
    table[["human_ortholog", "human_genes"]] = table[["human_ortholog", "human_genes"]].fillna("")
    return table[OUTPUT_COLUMNS].sort_values("gene_id").reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build resources/mouse_genes.tsv.gz from two MGI reports.")
    parser.add_argument("--from-dir", type=Path, help=f"folder with already downloaded {GENE_REPORT} and {HOMOLOGY_REPORT}")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as scratch:  # downloaded reports are deleted afterwards
        folder = args.from_dir or Path(scratch)
        gene_path = folder / GENE_REPORT if args.from_dir else download(GENE_REPORT, folder)
        homology_path = folder / HOMOLOGY_REPORT if args.from_dir else download(HOMOLOGY_REPORT, folder)
        genes = read_gene_report(gene_path)
        homology = pd.read_csv(homology_path, sep="\t", dtype=str, keep_default_na=False)
    table = build(genes, homology)
    # mtime 0: the same table gives the same file, so git only sees real changes
    table.to_csv(OUTPUT, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})

    kinds = table["ortholog_type"].value_counts()
    print(f"  {len(table):,} mouse genes ({table['mgi_id'].nunique():,} MGI genes) -> {show_path(OUTPUT)}")
    print("  human orthologs: " + ", ".join(f"{kind} {kinds.get(kind, 0):,}" for kind in
                                           ["one-to-one", "one-to-many", "many-to-one", "many-to-many", "none"]))


if __name__ == "__main__":
    main()
