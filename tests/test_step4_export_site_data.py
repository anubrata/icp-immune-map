import json
import math

import pandas as pd

import step4_export_site_data
from common import is_significant, load_config, read_tsv
from conftest import HEALTHY, ICP, RENAMED_ID, gene_id, write_config
from step4_export_site_data import clean_description, row_zscores


def strict_json(path):
    """Parse JSON, failing on NaN/Infinity (which browsers cannot read)."""
    def reject(constant):
        raise ValueError(f"invalid JSON constant {constant}")
    return json.loads(path.read_text(), parse_constant=reject)


def test_exports_only_significant_genes(pipeline_run):
    cfg = load_config(pipeline_run["config"])
    data = strict_json(pipeline_run["site_data"] / "SYNTH.json")
    annotated = read_tsv(pipeline_run["results"] / "de_all_genes_annotated.tsv")
    expected = annotated.index[is_significant(annotated, cfg["cutoffs"])]

    assert sorted(g["id"] for g in data["genes"]) == sorted(expected)
    assert all(g["pvalue"] < 0.05 and abs(g["log2fc"]) >= 1 for g in data["genes"])
    assert [s["id"] for s in data["samples"]] == HEALTHY + ICP
    assert all(len(g["normalized"]) == len(g["z"]) == 10 for g in data["genes"])
    assert data["genes_tested"] == len(annotated)
    assert data["design"] == "~ fetal_sex + condition"

    # the arrow chart's two levels are exactly the DESeq2 fold change apart
    for gene in data["genes"]:
        levels = gene["group_estimates"]
        assert set(levels) == {"Healthy", "ICP"}
        assert abs(math.log2(levels["ICP"] / levels["Healthy"]) - gene["log2fc"]) < 0.01


def test_immune_cells_and_their_changed_marker_genes(pipeline_run):
    data = strict_json(pipeline_run["site_data"] / "SYNTH.json")
    assert data["immport"] is False and data["immune_lists"] == []      # ImmPort is off by default
    genes = {g["id"]: g for g in data["genes"]}
    renamed = genes[RENAMED_ID]
    assert renamed["current_symbol"] == "ZGRF1" and renamed["atlas_markers"] == ["1.3"]
    assert renamed["immune"] and renamed["immune_lists"] == []

    grid = {c["id"]: c for c in data["cell_marker_changes"]}
    assert "13.1" not in grid                                             # functional states are not cell types
    assert grid["1.3"]["name"] == "Renamed-marker cells" and grid["1.3"]["genes"] == [RENAMED_ID]
    assert grid["1.1"]["genes"] and set(grid["1.1"]["genes"]) <= {gene_id(i) for i in (0, 2, 3, 4)}
    assert all(set(c["genes"]) <= set(genes) for c in grid.values())     # only significant genes
    assert list(grid) == sorted(grid, key=lambda i: tuple(map(int, i.split("."))))  # atlas order


def test_cell_types_reach_the_site(pipeline_run):
    data = strict_json(pipeline_run["site_data"] / "SYNTH.json")
    assert data["atlas_sections"] == ["1. Major lineages", "13. Common immune-cell functional states across lineages"]
    assert [c["id"] for c in data["cell_types"]] == ["1.1", "1.2", "1.3", "13.1"]  # "1.10" has too few markers
    planted = data["cell_types"][0]
    assert [m["symbol"] for m in planted["markers"]] == ["GENE0", "GENE2", "GENE3", "GENE4"]
    assert planted["effect"] > 1 and planted["ci"][0] > 0 and len(planted["scores"]) == 10
    assert data["atlas_populations"]["1.10"]["name"] == "One-marker cells"
    assert data["composition_check"]["name"] == "Decidua" and len(data["composition_check"]["scores"]) == 10

    genes = {g["id"]: g for g in data["genes"]}
    assert genes[gene_id(2)]["atlas_markers"] == ["1.1"] and genes[gene_id(2)]["immune"]
    assert genes[gene_id(5)]["atlas_markers"] == ["13.1"] and not genes[gene_id(5)]["immune"]

    cells = pd.read_csv(pipeline_run["site_data"] / "SYNTH_cell_type_scores.csv", dtype={"population_id": str})
    assert cells["population_id"].tolist()[0] == "1.1" and "r_with_Decidua" in cells


def test_csv_and_dataset_index(pipeline_run):
    data = strict_json(pipeline_run["site_data"] / "SYNTH.json")
    table = pd.read_csv(pipeline_run["site_data"] / "SYNTH_significant_genes.csv")
    assert len(table) == len(data["genes"])
    assert "H1 (Healthy, normalized)" in table.columns
    assert "immune_cell_markers" in table and "immport_lists" not in table
    index = strict_json(pipeline_run["site_data"] / "datasets.json")
    assert index == [{"id": "SYNTH", "title": "Synthetic test data", "species": "human", "tissue": "Placenta",
                      "file": "SYNTH.json", "significant_genes": len(data["genes"])}]
    assert data["species"] == "human" and data["notes"] is None
    assert all(g["human_ortholog"] is None for g in data["genes"])


def test_no_significant_genes_still_exports(pipeline_run, tmp_path):
    strict = write_config(pipeline_run["root"] / "SYNTH_strict.yaml", pipeline_run["data_file"],
                          cutoffs={"padj": 1e-300, "abs_log2fc": 50})
    step4_export_site_data.run(strict, pipeline_run["results_root"], tmp_path)
    data = strict_json(tmp_path / "SYNTH.json")
    assert data["genes"] == [] and data["immune_lists"] == [] and data["cell_marker_changes"] == []
    assert data["cutoffs"] == {"padj": 1e-300, "abs_log2fc": 50}


def test_helpers():
    assert clean_description("fibronectin 1 [Source:HGNC Symbol;Acc:HGNC:3778]") == "fibronectin 1"
    assert clean_description(float("nan")) is None
    z = row_zscores(pd.DataFrame({"a": [1.0, 5.0], "b": [3.0, 5.0]}))
    assert z.loc[1].tolist() == [0.0, 0.0]  # a gene with no variation gets zeros, not NaN
    assert round(z.loc[0, "a"], 3) == -0.707
