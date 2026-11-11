"""Mouse datasets: a header-less count file, gene annotation from a reference table, atlas markers
matched only through one-to-one human orthologs, the site export and the cross-check; and the
step 0 builder of the mouse gene table."""

import json

import numpy as np
import pandas as pd
import pytest

import step0_fetch_mouse_genes
import step1_load_data
import step2_differential_expression
import step3_immune_annotation
import step4_export_site_data
import step5_check_site_data
from common import ConfigError, load_config, read_tsv
from conftest import write_config, write_resources

VEHICLE, ESTROGEN = ["V1", "V2", "V3", "V4"], ["E1", "E2", "E3", "E4"]
N_GENES = 320
DECOY = 2          # mouse gene spelled like the human marker GENE2, but without a one-to-one ortholog
UNANNOTATED = 319  # a gene missing from the reference table


def mouse_id(i: int) -> str:
    return f"ENSMUSG9{i:010d}"


@pytest.fixture(scope="module")
def mouse_run(tmp_path_factory):
    """Synthetic GSE314547-like data: genes 0-9 4x higher with estrogen, 10-19 4x lower."""
    root = tmp_path_factory.mktemp("mouse")
    rng = np.random.default_rng(5)
    mean = rng.lognormal(6.0, 0.8, N_GENES)
    effect = np.ones(N_GENES)
    effect[:10], effect[10:20] = 4.0, 0.25
    columns = VEHICLE + ESTROGEN + ["CM1"]  # CM1: a group the config leaves out
    mu = np.column_stack([mean * (effect if c in ESTROGEN else 1.0) for c in columns])
    counts = rng.negative_binomial(n=20, p=20 / (20 + mu))

    # Header with one name fewer than the data rows: the gene ID column has no name (as in GSE314547).
    data_file = root / "mouse_counts.txt"
    lines = ["\t".join(columns)] + ["\t".join([mouse_id(i), *map(str, row)]) for i, row in enumerate(counts)]
    data_file.write_text("\n".join(lines) + "\n")

    annotation = pd.DataFrame({
        "gene_id": [mouse_id(i) for i in range(N_GENES)],
        "gene_symbol": [f"Gene{i}" for i in range(N_GENES)],
        "description": [f"made-up mouse gene {i}" for i in range(N_GENES)],
        "chromosome": [str(1 + i % 19) for i in range(N_GENES)],
        "biotype": "protein coding",
        "mgi_id": [f"MGI:{i}" for i in range(N_GENES)],
        "human_ortholog": [f"GENE{i}" for i in range(N_GENES)],
        "ortholog_type": "one-to-one",
        "human_genes": [f"GENE{i}" for i in range(N_GENES)],
    })
    annotation.loc[DECOY, ["gene_symbol", "human_ortholog", "ortholog_type", "human_genes"]] = [
        "GENE2", "", "many-to-one", "GENE2"]
    annotation = annotation.drop(index=UNANNOTATED)
    annotation_file = root / "mouse_genes.tsv.gz"
    annotation.to_csv(annotation_file, sep="\t", index=False)

    config = write_config(
        root / "MOUSE.yaml", data_file, id="MOUSE", species="mouse", tissue="Liver",
        notes="A synthetic mouse model.", columns={"gene_id": "gene_id"}, gene_annotation=str(annotation_file),
        samples={**{s: "Vehicle" for s in VEHICLE}, **{s: "Estrogen" for s in ESTROGEN}},
        reference_group="Vehicle", compare_group="Estrogen", fetal_sex=None, composition_check=None)
    resources, results_root, site_data = write_resources(root / "resources"), root / "results", root / "site_data"
    step1_load_data.run(config, results_root)
    step2_differential_expression.run(config, results_root)
    step3_immune_annotation.run(config, results_root, resources)
    step4_export_site_data.run(config, results_root, site_data)
    return {"config": config, "results": results_root / "MOUSE", "results_root": results_root,
            "site_data": site_data, "resources": resources}


def test_headerless_file_and_gene_annotation(mouse_run):
    counts = read_tsv(mouse_run["results"] / "counts.tsv")
    assert list(counts.columns) == VEHICLE + ESTROGEN and len(counts) == N_GENES  # CM1 left out
    genes = read_tsv(mouse_run["results"] / "genes.tsv").fillna("")
    assert genes.loc[mouse_id(0), "gene_symbol"] == "Gene0" and genes.loc[mouse_id(0), "human_ortholog"] == "GENE0"
    assert genes.loc[mouse_id(0), "chromosome"] in (1, "1")
    assert genes.loc[mouse_id(UNANNOTATED), "gene_symbol"] == mouse_id(UNANNOTATED)  # not in the reference: keeps its ID
    assert genes.loc[mouse_id(DECOY), "human_ortholog"] == ""


def test_atlas_markers_only_through_orthologs(mouse_run):
    annotated = read_tsv(mouse_run["results"] / "de_all_genes_annotated.tsv")
    annotated["atlas_markers"] = annotated["atlas_markers"].fillna("")
    assert annotated.loc[mouse_id(0), "atlas_markers"] == "1.1" and annotated.loc[mouse_id(0), "immune"]
    # spelled exactly like the marker GENE2, but no one-to-one ortholog: not a marker
    assert annotated.loc[mouse_id(DECOY), "atlas_markers"] == "" and not annotated.loc[mouse_id(DECOY), "immune"]
    cells = read_tsv(mouse_run["results"] / "cell_type_scores.tsv")
    assert cells.loc["1.1", "markers_found"] == "GENE0,GENE3,GENE4"  # GENE2 cannot be measured
    assert cells.loc["1.1", "effect"] > 0 and cells.loc["1.1", "pvalue"] < 0.05
    info = json.loads((mouse_run["results"] / "step2_info.json").read_text())
    assert info["design"] == "~ condition" and info["chr_y_genes_excluded"] == 0


def test_site_data_and_cross_check(mouse_run, tmp_path):
    data = json.loads((mouse_run["site_data"] / "MOUSE.json").read_text())
    assert data["species"] == "mouse" and data["notes"] == "A synthetic mouse model."
    assert all(s["fetal_sex"] is None for s in data["samples"]) and data["sex_counts"] is None
    by_id = {g["id"]: g for g in data["genes"]}
    assert by_id[mouse_id(0)]["symbol"] == "Gene0" and by_id[mouse_id(0)]["human_ortholog"] == "GENE0"
    index = json.loads((mouse_run["site_data"] / "datasets.json").read_text())
    assert index[0]["species"] == "mouse"
    assert step5_check_site_data.run(mouse_run["config"], mouse_run["results_root"], mouse_run["site_data"],
                                     mouse_run["resources"]) == 0

    # A wrong ortholog on the site must fail the cross-check.
    tampered = tmp_path / "site_data"
    tampered.mkdir()
    by_id[mouse_id(0)]["human_ortholog"] = "GENE1"
    (tampered / "MOUSE.json").write_text(json.dumps(data))
    assert step5_check_site_data.run(mouse_run["config"], mouse_run["results_root"], tampered,
                                     mouse_run["resources"]) > 0


def test_config_rules(tmp_path):
    data_file = tmp_path / "x.txt"
    with pytest.raises(ConfigError, match="gene_annotation"):
        load_config(write_config(tmp_path / "a.yaml", data_file, species="mouse", fetal_sex=None))
    with pytest.raises(ConfigError, match="fetal_sex: infer"):
        load_config(write_config(tmp_path / "b.yaml", data_file, species="mouse", gene_annotation="g.tsv"))
    with pytest.raises(ConfigError, match="species must be"):
        load_config(write_config(tmp_path / "c.yaml", data_file, species="zebrafish"))
    assert load_config(write_config(tmp_path / "d.yaml", data_file))["species"] == "human"


def test_mouse_gene_table_builder():
    """One-to-one only: counted over the whole homology report, in both directions."""
    genes = pd.DataFrame({
        "mgi_id": ["MGI:1", "MGI:1", "MGI:2", "MGI:3", "MGI:4", "MGI:5", "MGI:6"],
        "gene_symbol": ["A", "A", "B", "C", "D", "E", "F"],
        "description": "x", "cm_position": "", "chromosome": "1",
        "gene_id": ["ENSMUSG1", "ENSMUSG1b", "ENSMUSG2", "ENSMUSG3", "ENSMUSG4", "ENSMUSG5", "ENSMUSG6"],
        "ensembl_transcripts": "", "ensembl_proteins": "",
        "feature_type": ["protein coding gene"] * 6 + ["lncRNA gene"],
        "start": "", "end": "", "strand": "", "biotypes": "",
    })
    rows = [("c1", "10090", "A", "MGI:1"), ("c1", "9606", "AA", ""),        # one-to-one
            ("c2", "10090", "B", "MGI:2"), ("c2", "9606", "BB1", ""), ("c2", "9606", "BB2", ""),  # one-to-many
            ("c3", "10090", "C", "MGI:3"), ("c3", "9606", "CC", ""),        # CC also pairs with D...
            ("c4", "10090", "D", "MGI:4"), ("c4", "9606", "CC", ""),        # ...in another class: many-to-one
            ("c6", "10090", "F", "MGI:6"), ("c6", "10090", "F", "MGI:6"), ("c6", "9606", "FF", "")]  # F listed twice
    homology = pd.DataFrame(rows, columns=["DB Class Key", "NCBI Taxon ID", "Symbol", "Mouse MGI ID"])
    table = step0_fetch_mouse_genes.build(genes, homology).set_index("gene_id")
    assert list(table.columns) == step0_fetch_mouse_genes.OUTPUT_COLUMNS[1:]
    assert table.loc["ENSMUSG1", "human_ortholog"] == table.loc["ENSMUSG1b", "human_ortholog"] == "AA"
    assert table.loc["ENSMUSG2", ["human_ortholog", "ortholog_type", "human_genes"]].tolist() == ["", "one-to-many", "BB1;BB2"]
    assert table.loc["ENSMUSG3", "ortholog_type"] == table.loc["ENSMUSG4", "ortholog_type"] == "many-to-one"
    assert table.loc["ENSMUSG3", "human_ortholog"] == ""
    assert table.loc["ENSMUSG5", "ortholog_type"] == "none"
    assert table.loc["ENSMUSG6", ["human_ortholog", "ortholog_type", "biotype"]].tolist() == ["FF", "one-to-one", "lncRNA"]
