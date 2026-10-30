import gzip

import pandas as pd
import pytest
import yaml

from common import REPO_ROOT, ConfigError, DataError, load_config
from conftest import MALE, SAMPLES, simulate_counts, write_config, write_gse_like_file
from step1_load_data import infer_fetal_sex, load_dataset, read_table


@pytest.fixture
def gse_like(tmp_path):
    counts, _, _ = simulate_counts(n_genes=30, n_changed=4)
    data_file = tmp_path / "counts.txt"
    write_gse_like_file(data_file, counts)
    return tmp_path, data_file


def test_loads_the_gse282182_layout(gse_like):
    folder, data_file = gse_like
    counts, genes, samples = load_dataset(load_config(write_config(folder / "c.yaml", data_file)))
    assert list(counts.columns) == SAMPLES
    assert len(counts) == 30 + 6
    assert (counts.dtypes == "int64").all()
    assert genes.loc["ENSG00000000002", "description"] == "duplicated symbol, quoted"  # quotes removed
    assert genes["gene_symbol"].tolist().count("GENE1") == 2  # duplicated symbols are fine, IDs are unique
    assert samples["group"].value_counts().to_dict() == {"Healthy": 5, "ICP": 5}
    assert samples.loc["H1", "library_size"] == counts["H1"].sum()


def test_infers_fetal_sex_from_chrY_genes(gse_like):
    folder, data_file = gse_like
    counts, genes, _ = load_dataset(load_config(write_config(folder / "c.yaml", data_file)))
    calls = infer_fetal_sex(counts, genes)
    assert set(calls.index[calls["fetal_sex"] == "M"]) == MALE


def test_fetal_sex_inference_needs_chrY_genes(gse_like):
    folder, data_file = gse_like
    counts, genes, _ = load_dataset(load_config(write_config(folder / "c.yaml", data_file)))
    no_y = ~genes["gene_symbol"].isin(["RPS4Y1", "DDX3Y"])
    with pytest.raises(DataError, match="chrY marker"):
        infer_fetal_sex(counts[no_y], genes[no_y])


def test_reads_gzipped_csv_with_other_column_names(tmp_path):
    data_file = tmp_path / "other_layout.csv.gz"
    with gzip.open(data_file, "wt") as f:
        f.write("Geneid,Symbol,ctrl_a,ctrl_b,icp_a,icp_b\n"
                "ENSG1,IL6,10,12,40,44\n"
                "ENSG2,,5,6,7,8\n")
    config = tmp_path / "c.yaml"
    config.write_text(yaml.safe_dump({
        "id": "OTHER", "file": str(data_file),
        "columns": {"gene_id": "Geneid", "gene_symbol": "Symbol"},
        "samples": {"ctrl_a": "Healthy", "ctrl_b": "Healthy", "icp_a": "ICP", "icp_b": "ICP"},
        "reference_group": "Healthy", "compare_group": "ICP",
    }))
    counts, genes, samples = load_dataset(load_config(config))
    assert counts.loc["ENSG1"].tolist() == [10, 12, 40, 44]
    assert genes.loc["ENSG2", "gene_symbol"] == "ENSG2"  # missing symbol falls back to the gene ID
    assert "fetal_sex" not in samples


def test_missing_sample_column_is_named_in_the_error(gse_like):
    folder, data_file = gse_like
    samples = {s: "Healthy" if s.startswith("H") else "ICP" for s in SAMPLES} | {"H9": "Healthy"}
    config = write_config(folder / "c.yaml", data_file, samples=samples)
    with pytest.raises(DataError, match="H9"):
        load_dataset(load_config(config))


def test_duplicate_gene_ids_are_rejected(tmp_path):
    data_file = tmp_path / "dup.tsv"
    data_file.write_text("id\ta\tb\nG1\t1\t2\nG1\t3\t4\n")
    config = tmp_path / "c.yaml"
    config.write_text(yaml.safe_dump({"id": "D", "file": str(data_file), "columns": {"gene_id": "id"},
                                      "samples": {"a": "Healthy", "b": "ICP"},
                                      "reference_group": "Healthy", "compare_group": "ICP"}))
    with pytest.raises(DataError, match="more than once"):
        load_dataset(load_config(config))


def test_normalized_values_are_rejected(tmp_path):
    data_file = tmp_path / "tpm.tsv"
    data_file.write_text("id\ta\tb\nG1\t1.5\t2.25\nG2\t3\t4\n")
    config = tmp_path / "c.yaml"
    config.write_text(yaml.safe_dump({"id": "T", "file": str(data_file), "columns": {"gene_id": "id"},
                                      "samples": {"a": "Healthy", "b": "ICP"},
                                      "reference_group": "Healthy", "compare_group": "ICP"}))
    with pytest.raises(DataError, match="raw read counts"):
        load_dataset(load_config(config))


def test_config_errors_are_explained(gse_like):
    folder, data_file = gse_like
    with pytest.raises(ConfigError, match="compare_group"):
        load_config(write_config(folder / "c.yaml", data_file, compare_group="Severe ICP"))
    with pytest.raises(ConfigError, match="not both"):
        load_config(write_config(folder / "c.yaml", data_file, cutoffs={"pvalue": 0.05, "padj": 0.05}))
    cfg = load_config(write_config(folder / "c.yaml", data_file, cutoffs={"abs_log2fc": 2}))
    assert cfg["cutoffs"] == {"abs_log2fc": 2, "pvalue": 0.05}


def test_real_config_matches_the_data_file_header():
    """config/GSE282182.yaml names columns that really exist (reads only the header line)."""
    cfg = load_config(REPO_ROOT / "config" / "GSE282182.yaml")
    data_file = REPO_ROOT / cfg["file"]
    if not data_file.exists():
        pytest.skip("data file not present")
    header = read_table(data_file, {"nrows": 0}).columns
    named = [cfg["columns"][k] for k in cfg["columns"]] + list(cfg["samples"])
    assert [c for c in named if c not in header] == []
    assert isinstance(header, pd.Index)
