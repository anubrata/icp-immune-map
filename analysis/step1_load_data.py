"""
Step 1 - Load one dataset and save it as standard tables.

Whatever the original files look like, they end up in the same shape here, so the later steps
never need to know the original layout. The dataset config says what is where.

RNA-seq counts (data_type: counts), one table with a column per sample:
  Input:   the data file named in the config (e.g. data/GSE282182_gene_count.txt.gz)
  Outputs: results/<dataset>/counts.tsv          raw read counts: one row per gene, one column per sample
           results/<dataset>/genes.tsv           gene annotation: symbol, description, chromosome, biotype

Agilent one-color microarrays (data_type: microarray), one Feature Extraction file per sample:
  Input:   the files named in the config, plus the platform's probe annotation (from step 0)
  Outputs: results/<dataset>/intensities.tsv       processed signal: one row per probe, one column per array
           results/<dataset>/above_background.tsv  whether each probe is clearly above background on each array
           results/<dataset>/probes.tsv            probe annotation: gene symbol, gene name, chromosome

Both:      results/<dataset>/samples.tsv           one row per sample: group, (library size), fetal sex

Run:     uv run analysis/step1_load_data.py config/GSE282182.yaml
"""

from __future__ import annotations

import gzip
import io
from pathlib import Path

import numpy as np
import pandas as pd

from common import (RESULTS_ROOT, DataError, data_path, load_config, repo_path, results_dir, show_path,
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


def read_agilent(path: Path) -> pd.DataFrame:
    """One Agilent Feature Extraction file -> its FEATURES table (one row per spot on the array).

    The file has three tab-separated sections (FEPARAMS, STATS, FEATURES); only FEATURES holds
    the measurements: ProbeName, ControlType (0 = gene probe), gProcessedSignal (background-corrected
    signal) and gIsWellAboveBG (1 = clearly above background).
    """
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        lines = f.read().splitlines()
    start = next((i for i, line in enumerate(lines) if line.startswith("FEATURES\t")), None)
    if start is None:
        raise DataError(f"{path.name}: no FEATURES section; is this an Agilent Feature Extraction file?")
    spots = pd.read_csv(io.StringIO("\n".join(lines[start:])), sep="\t", dtype={"ProbeName": str}, low_memory=False)
    missing = [c for c in ("ProbeName", "ControlType", "gProcessedSignal", "gIsWellAboveBG") if c not in spots]
    if missing:
        raise DataError(f"{path.name}: missing column(s) {missing} in the FEATURES section")
    return spots


def load_microarray(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read every array; return (signal, above background, probe annotation, samples).

    Only gene probes are kept (control spots are dropped). Probes printed several times on the
    array are averaged; a probe counts as above background when most of its spots are.
    """
    signal, above = {}, {}
    for sample, file in cfg["files"].items():
        path = repo_path(file)
        if not path.exists():
            raise DataError(f"data file not found for sample {sample}: {path}")
        spots = read_agilent(path)
        spots = spots[spots["ControlType"] == 0]
        by_probe = spots.groupby("ProbeName")
        signal[sample] = by_probe["gProcessedSignal"].mean()
        above[sample] = by_probe["gIsWellAboveBG"].mean() >= 0.5
    signal, above = pd.DataFrame(signal), pd.DataFrame(above)
    if signal.isna().to_numpy().any():
        raise DataError("the arrays do not all have the same probes; are they from the same platform?")
    signal.index.name = above.index.name = "probe_id"

    annotation = pd.read_csv(repo_path(cfg["platform_annotation"]), sep="\t", dtype=str,
                             keep_default_na=False).set_index("probe_id")
    probes = annotation.reindex(signal.index).fillna("")
    samples = pd.DataFrame({"group": [cfg["samples"][s] for s in signal.columns]},
                           index=pd.Index(signal.columns, name="sample"))
    return signal, above, probes, samples


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


def run_microarray(cfg: dict, results_root: Path) -> None:
    signal, above, probes, samples = load_microarray(cfg)
    group_sizes = ", ".join(f"{g}: {n}" for g, n in samples["group"].value_counts(sort=False).items())
    with_gene = probes["gene_symbol"] != ""
    print(f"{cfg['id']}: {cfg['title']}")
    print(f"  arrays:  {len(samples)} ({group_sizes}), one file each")
    print(f"  probes:  {len(signal):,} gene probes, {with_gene.sum():,} of them for "
          f"{probes.loc[with_gene, 'gene_symbol'].nunique():,} named genes (platform annotation: "
          f"{show_path(repo_path(cfg['platform_annotation']))})")
    if isinstance(cfg.get("fetal_sex"), dict):
        samples["fetal_sex"] = [cfg["fetal_sex"][s] for s in samples.index]
        print("\nFetal sex as listed in the config:")
        print(samples[["group", "fetal_sex"]].to_string())

    out = results_dir(cfg, results_root)
    write_tsv(signal.round(3), out / "intensities.tsv")
    write_tsv(above, out / "above_background.tsv")
    write_tsv(probes, out / "probes.tsv")
    write_tsv(samples, out / "samples.tsv")
    print(f"\nWrote {show_path(out)}/intensities.tsv, above_background.tsv, probes.tsv, samples.tsv")


def run(config_path: str | Path, results_root: Path = RESULTS_ROOT) -> None:
    cfg = load_config(config_path)
    if cfg["data_type"] == "microarray":
        return run_microarray(cfg, results_root)
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
