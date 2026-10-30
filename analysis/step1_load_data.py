"""
Step 1 - Load one dataset and save it as three standard tables.

Whatever the original file looks like, it ends up in the same shape here, so steps 2-4
never need to know the original layout. The dataset config says which column holds what.

Input:   the data file named in the config (e.g. data/GSE282182_gene_count.txt.gz)
Outputs: results/<dataset>/counts.tsv   raw read counts: one row per gene, one column per sample
         results/<dataset>/genes.tsv    gene annotation: symbol, description, chromosome, biotype
         results/<dataset>/samples.tsv  one row per sample: group, library size, fetal sex

Run:     uv run analysis/step1_load_data.py config/GSE282182.yaml
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from common import (RESULTS_ROOT, DataError, data_path, load_config, results_dir, show_path,
                    step_arguments, write_tsv)

# Optional annotation columns (config `columns:` keys) that the website knows how to show.
ANNOTATION_KEYS = ["description", "chromosome", "biotype"]

COMPRESSION_SUFFIXES = {".gz", ".bz2", ".zip", ".xz"}

# chrY genes that are expressed in male tissue and absent from female tissue. The placenta is
# fetal tissue, so reads on these genes reveal the sex of the baby.
CHR_Y_MARKERS = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY"]
MALE_CPM_THRESHOLD = 10  # chrY-marker reads per million: males are in the hundreds, females ~0


def read_table(path: Path, read_options: dict | None = None) -> pd.DataFrame:
    """Read a .csv/.tsv/.txt file (optionally compressed) or an .xlsx sheet.

    The delimiter follows the extension (.csv -> comma, otherwise tab). Anything in the
    config's `read_options` is passed straight to pandas, for layouts that need more.
    """
    options = dict(read_options or {})
    suffixes = [s.lower() for s in path.suffixes]
    if suffixes and suffixes[-1] in (".xlsx", ".xls"):
        return pd.read_excel(path, **options)
    extension = next((s for s in reversed(suffixes) if s not in COMPRESSION_SUFFIXES), "")
    options.setdefault("sep", "," if extension == ".csv" else "\t")
    return pd.read_csv(path, **options)


def load_dataset(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read the data file and return (counts, genes, samples), checking it matches the config."""
    path = data_path(cfg)
    if not path.exists():
        raise DataError(f"data file not found: {path}")
    table = read_table(path, cfg["read_options"])
    columns = cfg["columns"]
    sample_columns = list(cfg["samples"])

    # 1. Every column the config mentions must exist in the file.
    optional = [columns[key] for key in ["gene_symbol", *ANNOTATION_KEYS] if columns.get(key)]
    missing = [c for c in [columns["gene_id"], *sample_columns, *optional] if c not in table.columns]
    if missing:
        raise DataError(f"{path.name}: column(s) {missing} not found. Columns in the file: {list(table.columns)}")

    # 2. Gene IDs must be unique: they identify each gene through every later step.
    gene_ids = table[columns["gene_id"]].astype(str).str.strip()
    duplicated = gene_ids[gene_ids.duplicated()].unique()
    if len(duplicated):
        raise DataError(f"{path.name}: {len(duplicated)} gene IDs appear more than once "
                        f"(e.g. {list(duplicated[:3])}); each row needs a unique gene ID")

    # 3. Sample columns must hold raw read counts (whole numbers >= 0). DESeq2 models raw
    #    counts, so already-normalized values such as TPM or FPKM cannot be used.
    counts = table[sample_columns].apply(pd.to_numeric, errors="coerce")
    if counts.isna().to_numpy().any():
        bad = counts.columns[counts.isna().any()].tolist()
        raise DataError(f"{path.name}: non-numeric or empty values in sample column(s) {bad}")
    values = counts.to_numpy(dtype=float)
    if (values < 0).any() or not np.allclose(values, np.round(values)):
        raise DataError(f"{path.name}: sample columns must be raw read counts (whole numbers >= 0); "
                        "normalized values such as TPM/FPKM cannot be used with DESeq2")
    counts = counts.round().astype("int64")
    counts.index = pd.Index(gene_ids.to_numpy(), name="gene_id")

    # Gene annotation. A missing symbol falls back to the gene ID.
    genes = pd.DataFrame(index=counts.index)
    if columns.get("gene_symbol"):
        symbols = table[columns["gene_symbol"]].fillna("").astype(str).str.strip()
        genes["gene_symbol"] = symbols.where(symbols != "", gene_ids).to_numpy()
    else:
        genes["gene_symbol"] = gene_ids.to_numpy()
    for key in ANNOTATION_KEYS:
        if columns.get(key):
            genes[key] = table[columns[key]].to_numpy()

    samples = pd.DataFrame({"group": [cfg["samples"][s] for s in sample_columns]},
                           index=pd.Index(sample_columns, name="sample"))
    samples["library_size"] = counts.sum().to_numpy()
    return counts, genes, samples


def infer_fetal_sex(counts: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    """Call each sample male (M) or female (F) from its reads on chrY marker genes.

    Returns one row per sample with the evidence (reads per million on chrY markers, and on
    XIST, which is high in females) and the call.
    """
    markers = genes.index[genes["gene_symbol"].isin(CHR_Y_MARKERS)]
    if len(markers) == 0:
        raise DataError(f"fetal_sex: infer needs at least one chrY marker gene ({', '.join(CHR_Y_MARKERS)}), "
                        "but none are in this dataset. List each sample's sex in the config instead, "
                        "or remove the fetal_sex setting.")
    library_size = counts.sum()
    calls = pd.DataFrame(index=counts.columns)
    calls["chrY_per_million"] = (counts.loc[markers].sum() / library_size * 1e6).round(1)
    xist = genes.index[genes["gene_symbol"] == "XIST"]
    if len(xist):
        calls["XIST_per_million"] = (counts.loc[xist].sum() / library_size * 1e6).round(1)
    calls["fetal_sex"] = np.where(calls["chrY_per_million"] > MALE_CPM_THRESHOLD, "M", "F")
    return calls


def run(config_path: str | Path, results_root: Path = RESULTS_ROOT) -> None:
    cfg = load_config(config_path)
    counts, genes, samples = load_dataset(cfg)

    group_sizes = ", ".join(f"{g}: {n}" for g, n in samples["group"].value_counts(sort=False).items())
    print(f"{cfg['id']}: {cfg['title']}")
    print(f"  file:    {show_path(data_path(cfg))}")
    print(f"  genes:   {len(counts):,}")
    print(f"  samples: {len(samples)} ({group_sizes})")

    if cfg.get("fetal_sex") == "infer":
        calls = infer_fetal_sex(counts, genes)
        samples["fetal_sex"] = calls["fetal_sex"]
        print(f"\nFetal sex from chrY genes (male if more than {MALE_CPM_THRESHOLD} chrY reads per million):")
        print(samples[["group"]].join(calls).to_string())
    elif isinstance(cfg.get("fetal_sex"), dict):
        samples["fetal_sex"] = [cfg["fetal_sex"][s] for s in samples.index]
        print("\nFetal sex as listed in the config:")
        print(samples[["group", "fetal_sex"]].to_string())

    out = results_dir(cfg, results_root)
    write_tsv(counts, out / "counts.tsv")
    write_tsv(genes, out / "genes.tsv")
    write_tsv(samples, out / "samples.tsv")
    print(f"\nWrote {show_path(out)}/counts.tsv, genes.tsv, samples.tsv")


if __name__ == "__main__":
    args = step_arguments("Step 1: load a dataset into standard tables.")
    run(args.config, args.results_root)
