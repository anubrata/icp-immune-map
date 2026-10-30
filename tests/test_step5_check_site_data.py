import json
import shutil

import step5_check_site_data
from conftest import SEXFLIP_ID


def test_all_checks_pass_on_correct_site_data(pipeline_run):
    failures = step5_check_site_data.run(pipeline_run["config"], pipeline_run["results_root"],
                                         pipeline_run["site_data"], pipeline_run["resources"])
    assert failures == 0


def test_a_wrong_number_is_caught(pipeline_run, tmp_path):
    site_data = tmp_path / "site_data"
    shutil.copytree(pipeline_run["site_data"], site_data)
    path = site_data / "SYNTH.json"
    data = json.loads(path.read_text())
    data["immune_lists"][0]["higher"] += 1           # one bar off by one
    data["genes"][0]["log2fc"] = -data["genes"][0]["log2fc"]  # one gene pointing the wrong way
    path.write_text(json.dumps(data))
    failures = step5_check_site_data.run(pipeline_run["config"], pipeline_run["results_root"],
                                         site_data, pipeline_run["resources"])
    assert failures >= 2


def test_sex_dependent_gene_is_flagged(pipeline_run):
    data = json.loads((pipeline_run["site_data"] / "SYNTH.json").read_text())
    assert data["sex_counts"] == {"F": {"Healthy": 4, "ICP": 2}, "M": {"Healthy": 1, "ICP": 3}}
    genes = {g["id"]: g for g in data["genes"]}
    flip = genes[SEXFLIP_ID]
    assert flip["log2fc"] < 0 < flip["plain_log2fc"]     # model: lower; plain averages: higher
    assert flip["raw_disagrees"]
    assert flip["sex_log2fc"]["F"] < 0 < flip["sex_log2fc"]["M"]
    assert all(g["raw_disagrees"] == (g["plain_log2fc"] * g["log2fc"] < 0) for g in genes.values())
