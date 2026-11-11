"""Shared helpers for the analysis steps: reading a dataset config and knowing where files go.

Every step script takes the same argument, the path to a dataset config such as
config/GSE282182.yaml, and reads/writes its tables under results/<dataset id>/.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_ROOT = REPO_ROOT / "results"
RESOURCES_DIR = REPO_ROOT / "resources"
SITE_DATA_DIR = REPO_ROOT / "site" / "data"

# The kinds of data the pipeline can analyse (config `data_type`):
#   counts      RNA-seq raw read counts, one table          -> DESeq2
#   microarray  Agilent one-color arrays, one file per sample -> limma-style moderated t-test
DATA_TYPES = ("counts", "microarray")

# Species the pipeline can analyse (config `species`). The immune cell atlas lists human genes; for
# other species each gene needs its human ortholog, which comes from the config's `gene_annotation`
# table (for mouse: resources/mouse_genes.tsv.gz, made by analysis/step0_fetch_mouse_genes.py).
SPECIES = ("human", "mouse")

# Expression tables written by step 2: (normalized values on the original scale, log-scale values
# used for the heatmap, PCA and cell-type scores)
EXPRESSION_FILES = {"counts": ("normalized_counts.tsv", "vst.tsv"),
                    "microarray": ("normalized_intensities.tsv", "log2_intensities.tsv")}


class ConfigError(ValueError):
    """The dataset config is missing something or contradicts itself."""


class DataError(ValueError):
    """The data file does not match what the config describes."""


def load_config(path: str | Path) -> dict:
    """Read a dataset config (YAML), check it, and fill in defaults for optional settings."""
    path = Path(path)
    with open(path) as f:
        cfg = yaml.safe_load(f) or {}

    def require(condition: bool, message: str) -> None:
        if not condition:
            raise ConfigError(f"{path}: {message}")

    for key in ("id", "samples", "reference_group", "compare_group"):
        require(key in cfg, f"missing required setting '{key}'")
    require(isinstance(cfg["samples"], dict) and len(cfg["samples"]) > 0,
            "samples must map each sample to its group")
    data_type = cfg.setdefault("data_type", "counts")
    require(data_type in DATA_TYPES, f"data_type must be one of {DATA_TYPES}")
    if data_type == "counts":
        for key in ("file", "columns"):
            require(key in cfg, f"missing required setting '{key}' (needed for data_type: counts)")
        require("gene_id" in (cfg["columns"] or {}), "columns.gene_id is required")
    else:
        for key in ("files", "platform_annotation"):
            require(key in cfg, f"missing required setting '{key}' (needed for data_type: microarray)")
        require(isinstance(cfg["files"], dict) and set(cfg["files"]) == set(cfg["samples"]),
                "files must give one data file for every sample in samples")

    species = cfg.setdefault("species", "human")
    require(species in SPECIES, f"species must be one of {SPECIES}")
    cfg.setdefault("gene_annotation", None)
    if species != "human":
        require(data_type == "counts", f"species: {species} is only supported for RNA-seq counts")
        require(cfg["gene_annotation"], f"species: {species} needs gene_annotation (a table with each gene's "
                "human ortholog, e.g. resources/mouse_genes.tsv.gz) to match genes to the human immune cell atlas")
        require(cfg.get("fetal_sex") != "infer", "fetal_sex: infer only works for human placenta (chrY marker genes)")

    groups = set(cfg["samples"].values())
    for key in ("reference_group", "compare_group"):
        require(cfg[key] in groups, f"{key} '{cfg[key]}' is not one of the sample groups {sorted(groups)}")
    require(cfg["reference_group"] != cfg["compare_group"], "reference_group and compare_group must differ")

    fetal_sex = cfg.get("fetal_sex")
    if isinstance(fetal_sex, dict):
        require(set(fetal_sex) == set(cfg["samples"]), "fetal_sex must list every sample exactly once")
        require(set(fetal_sex.values()) <= {"M", "F"}, "fetal_sex values must be M or F")
    else:
        require(fetal_sex in (None, "infer"),
                "fetal_sex must be 'infer', a map of sample -> M/F, or left out")
    require(not (fetal_sex == "infer" and data_type == "microarray"),
            "fetal_sex: infer only works for RNA-seq counts; list each array's fetal sex instead")
    cfg.setdefault("adjust_for_fetal_sex", True)
    require(isinstance(cfg["adjust_for_fetal_sex"], bool), "adjust_for_fetal_sex must be true or false")

    cutoffs = cfg.setdefault("cutoffs", {}) or {}
    require(not ("pvalue" in cutoffs and "padj" in cutoffs),
            "cutoffs takes either pvalue (unadjusted) or padj (FDR-adjusted), not both")
    if "padj" not in cutoffs:
        cutoffs.setdefault("pvalue", 0.05)
    cutoffs.setdefault("abs_log2fc", 1.0)
    cfg["cutoffs"] = cutoffs

    cfg.setdefault("use_immport", False)
    require(isinstance(cfg["use_immport"], bool), "use_immport must be true or false")

    composition = cfg.setdefault("composition_check", None)
    if composition is not None:
        require(isinstance(composition, dict) and composition.get("name") and composition.get("markers"),
                "composition_check needs a name and a list of markers")

    cfg.setdefault("title", cfg["id"])
    cfg.setdefault("tissue", None)
    cfg.setdefault("geo_accession", None)
    cfg.setdefault("notes", None)  # one or two sentences on the design and its limits, shown on the site
    # where the data comes from, for datasets not on GEO: {label: "BioProject PRJNA846869", url: "https://..."}
    source = cfg.setdefault("source", None)
    if source is not None:
        require(isinstance(source, dict) and source.get("label") and source.get("url"), "source needs a label and a url")
    cfg.setdefault("min_count", 10)
    cfg.setdefault("read_options", {})
    return cfg


def repo_path(path: str | Path) -> Path:
    """Paths in the configs are relative to the repository root (absolute paths are kept)."""
    path = Path(path)
    return path if path.is_absolute() else REPO_ROOT / path


def data_path(cfg: dict) -> Path:
    """The dataset file of a counts dataset."""
    return repo_path(cfg["file"])


def expression_files(cfg: dict) -> tuple[str, str]:
    """Names of step 2's (normalized, log-scale) expression tables for this dataset's data type."""
    return EXPRESSION_FILES[cfg["data_type"]]


def results_dir(cfg: dict, results_root: Path = RESULTS_ROOT) -> Path:
    """results/<dataset id>/, created if needed."""
    folder = Path(results_root) / cfg["id"]
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def is_significant(table: pd.DataFrame, cutoffs: dict) -> pd.Series:
    """True for genes that pass the config cutoffs: p-value (or padj) below the threshold
    AND an absolute log2 fold change of at least abs_log2fc (1 = 2-fold)."""
    p_column = "padj" if "padj" in cutoffs else "pvalue"
    return (table[p_column] < cutoffs[p_column]) & (table["log2FoldChange"].abs() >= cutoffs["abs_log2fc"])


def row_zscores(table: pd.DataFrame) -> pd.DataFrame:
    """Scale each row (gene) to mean 0 and standard deviation 1 across the columns (samples)."""
    sd = table.std(axis=1).replace(0, float("nan"))
    return table.sub(table.mean(axis=1), axis=0).div(sd, axis=0).fillna(0.0)


def describe_cutoffs(cutoffs: dict) -> str:
    """Human-readable cutoffs, e.g. 'p < 0.05 and |log2FC| >= 1'."""
    p_label = f"FDR (padj) < {cutoffs['padj']}" if "padj" in cutoffs else f"p < {cutoffs['pvalue']}"
    return f"{p_label} and |log2FC| >= {cutoffs['abs_log2fc']:g}"


# Columns holding IDs or lists of IDs: always text (an atlas ID like "1.10" must not become the number 1.1)
TEXT_COLUMNS = {name: str for name in ["immune_lists", "atlas_markers", "marker_gene_ids", "markers", "markers_found"]}


def read_tsv(path: Path) -> pd.DataFrame:
    """Read a table written by one of the steps. The first column (gene IDs, sample names,
    population IDs) becomes the row index and, like the ID columns above, is read as text."""
    return pd.read_csv(path, sep="\t", index_col=0, converters={0: str}, dtype=TEXT_COLUMNS)


def write_tsv(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, sep="\t")


def show_path(path: Path) -> str:
    """Path relative to the repository root when possible, for tidy log messages."""
    try:
        return str(Path(path).resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def step_arguments(description: str, extra: dict | None = None) -> argparse.Namespace:
    """Command-line arguments shared by the step scripts."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("config", help="dataset config, e.g. config/GSE282182.yaml")
    parser.add_argument("--results-root", type=Path, default=RESULTS_ROOT,
                        help="where intermediate tables are written (default: results/)")
    for flag, (default, help_text) in (extra or {}).items():
        parser.add_argument(flag, type=Path, default=default, help=help_text)
    return parser.parse_args()
