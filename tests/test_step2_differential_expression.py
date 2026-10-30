import json

import numpy as np
import pandas as pd
import pytest

from common import DataError, is_significant, load_config, read_tsv
from step2_differential_expression import design_formula, filter_low_counts


def samples_table(groups, sexes=None):
    table = pd.DataFrame({"group": groups}, index=[f"S{i}" for i in range(len(groups))])
    if sexes is not None:
        table["fetal_sex"] = sexes
    return table


def test_design_adds_fetal_sex_only_when_it_can_be_estimated():
    groups = ["Healthy"] * 3 + ["ICP"] * 3
    assert design_formula(samples_table(groups)) == "~ condition"
    assert design_formula(samples_table(groups, ["F"] * 6)) == "~ condition"
    assert design_formula(samples_table(groups, ["F", "F", "M", "M", "F", "M"])) == "~ fetal_sex + condition"
    with pytest.raises(DataError, match="confounded"):
        design_formula(samples_table(groups, ["F"] * 3 + ["M"] * 3))


def test_filter_keeps_genes_with_enough_reads_in_enough_samples():
    samples = samples_table(["Healthy"] * 2 + ["ICP"] * 3)
    counts = pd.DataFrame({"S0": [10, 50, 0], "S1": [10, 0, 0], "S2": [0, 0, 0], "S3": [0, 0, 0], "S4": [0, 9, 100]},
                          index=["two_samples", "one_sample", "low"])
    assert filter_low_counts(counts, samples, min_count=10).index.tolist() == ["two_samples"]


def test_finds_the_planted_changes(pipeline_run):
    cfg = load_config(pipeline_run["config"])
    table = read_tsv(pipeline_run["results"] / "de_all_genes.tsv")
    hits = table[is_significant(table, cfg["cutoffs"])]
    higher, lower = pipeline_run["higher"], pipeline_run["lower"]

    found = hits.index.intersection(higher + lower)
    assert len(found) >= 0.9 * len(higher + lower)
    assert (hits.loc[hits.index.intersection(higher), "log2FoldChange"] > 0).all()
    assert (hits.loc[hits.index.intersection(lower), "log2FoldChange"] < 0).all()

    unchanged = table.index.difference(higher + lower)
    false_hits = hits.index.intersection(unchanged)
    assert len(false_hits) <= 0.05 * len(unchanged)


def test_writes_the_expected_tables(pipeline_run):
    folder = pipeline_run["results"]
    table = read_tsv(folder / "de_all_genes.tsv")
    assert {"gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"} <= set(table.columns)
    assert "ENSG00000000001" not in table.index  # all-zero gene was not tested
    assert table["pvalue"].dropna().is_monotonic_increasing

    normalized = read_tsv(folder / "normalized_counts.tsv")
    vst = read_tsv(folder / "vst.tsv")
    assert normalized.shape == vst.shape == (len(table), 10)
    assert read_tsv(folder / "pca.tsv").shape == (10, 2)

    info = json.loads((folder / "step2_info.json").read_text())
    assert info["design"] == "~ fetal_sex + condition"
    assert info["genes_tested"] == len(table)


def test_group_estimates_match_the_fold_change(pipeline_run):
    """The arrow chart's levels: their ratio is exactly DESeq2's fold change for every gene."""
    folder = pipeline_run["results"]
    table = read_tsv(folder / "de_all_genes.tsv")
    estimates = read_tsv(folder / "group_estimates.tsv").loc[table.index]
    ratio = np.log2(estimates["ICP"] / estimates["Healthy"])
    assert np.allclose(ratio, table["log2FoldChange"], atol=1e-6)
