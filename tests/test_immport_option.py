"""use_immport: true adds ImmPort's immune gene lists on top of the atlas (steps 3-5 re-run with it on)."""

import json
import shutil

import step3_immune_annotation
import step4_export_site_data
import step5_check_site_data
from common import read_tsv
from conftest import RENAMED_ID, gene_id, write_config


def test_immport_on(pipeline_run, tmp_path):
    results_root, site_data = tmp_path / "results", tmp_path / "site_data"
    shutil.copytree(pipeline_run["results"], results_root / "SYNTH")  # reuse steps 1-2
    config = write_config(tmp_path / "SYNTH.yaml", pipeline_run["data_file"], use_immport=True)

    step3_immune_annotation.run(config, results_root, pipeline_run["resources"])
    step4_export_site_data.run(config, results_root, site_data)
    assert step5_check_site_data.run(config, results_root, site_data, pipeline_run["resources"]) == 0

    annotated = read_tsv(results_root / "SYNTH" / "de_all_genes_annotated.tsv")
    assert annotated.loc[RENAMED_ID, "immune_lists"] == "GO:0002218"
    assert annotated.loc[gene_id(1), "immune"]           # only in an ImmPort list: immune-related now
    lists = read_tsv(results_root / "SYNTH" / "immune_lists.tsv")
    assert lists.loc["GO:0002218", "genes_in_list"] == 4
    assert lists.loc["GO:0002218", "genes_tested"] == 4  # GENE0, both GENE1 rows, renamed gene

    data = json.loads((site_data / "SYNTH.json").read_text())
    assert data["immport"] is True
    genes = {g["id"]: g for g in data["genes"]}
    assert genes[RENAMED_ID]["immune_lists"] == ["GO:0002218"]
    bars = {entry["id"]: entry for entry in data["immune_lists"]}
    assert bars["GO:0002218"]["higher"] >= 2             # GENE0 and the renamed gene are planted "higher" genes
