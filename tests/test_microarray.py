"""The microarray path on synthetic Agilent files: reading, normalization, moderated t-test, steps 1-5."""

import json

import numpy as np
import pandas as pd
import pytest
from scipy import special

import step1_load_data
import step2_differential_expression
import step3_immune_annotation
import step4_export_site_data
import step5_check_site_data
from common import read_tsv
from conftest import write_config, write_resources
from moderated_t import fit_prior, moderated_t_test, trigamma_inverse
from step2_differential_expression import quantile_normalize

ARRAYS = {"H1": ("Healthy", "F"), "H2": ("Healthy", "F"), "I1": ("ICP", "M"), "I2": ("ICP", "M"),
          "I3": ("ICP", "F"), "I4": ("ICP", "F")}  # like GSE46157: the male placentas are all ICP


def write_agilent_file(path, spots: pd.DataFrame) -> None:
    """A minimal Agilent Feature Extraction file: FEPARAMS, STATS and FEATURES sections."""
    columns = ["FeatureNum", "ControlType", "ProbeName", "SystematicName", "gProcessedSignal", "gIsWellAboveBG"]
    lines = ["TYPE\ttext\tinteger", "FEPARAMS\tProtocol_Name\tScan_NumChannels", "DATA\tGE1_test\t1", "*",
             "TYPE\tfloat", "STATS\tgDarkOffsetAverage", "DATA\t1.0", "*",
             "TYPE\t" + "\t".join(["integer", "integer", "text", "text", "float", "boolean"]),
             "FEATURES\t" + "\t".join(columns)]
    lines += ["DATA\t" + "\t".join(str(row[c]) for c in columns) for _, row in spots.iterrows()]
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture(scope="module")
def microarray_run(tmp_path_factory):
    """Synthetic GSE46157-like dataset: 200 genes (GENE0-9 4x higher in ICP, GENE10-19 4x lower),
    a chrY gene, a gene with a bright and a dim probe, replicate spots and control spots."""
    root = tmp_path_factory.mktemp("microarray")
    rng = np.random.default_rng(11)
    genes = [f"GENE{i}" for i in range(200)]
    level = rng.uniform(6, 12, len(genes))
    annotation = pd.DataFrame({"probe_id": [f"P{i}" for i in range(200)], "gene_symbol": genes,
                               "gene_name": [f"made-up gene {i}" for i in range(200)], "ensembl_id": "",
                               "chromosome": "1"})
    extra = pd.DataFrame({"probe_id": ["PY", "P5dim"], "gene_symbol": ["DDX3Y", "GENE5"], "gene_name": ["DEAD-box Y", "dim probe"],
                          "ensembl_id": "", "chromosome": ["Y", "1"]})
    annotation = pd.concat([annotation, extra])
    annotation_file = root / "platform_probes.tsv"
    annotation.to_csv(annotation_file, sep="\t", index=False)

    files = {}
    for k, (array, (group, sex)) in enumerate(ARRAYS.items()):
        log2 = level + rng.normal(0, 0.25, len(genes))
        if group == "ICP":
            log2[:10] += 2
            log2[10:20] -= 2
        brightness = 1.5 if array == "I2" else 1.0  # one array scanned brighter: quantile normalization fixes it
        signal = list(2 ** log2 * brightness) + [300.0 if sex == "M" else 1.5, 2 ** (level[5] - 3)]
        spots = pd.DataFrame({"ProbeName": annotation["probe_id"].tolist(), "SystematicName": "NM_0",
                              "gProcessedSignal": signal, "ControlType": 0,
                              "gIsWellAboveBG": [1] * 200 + [1 if sex == "M" else 0, 1]})
        replicate = spots.iloc[[0]].assign(gProcessedSignal=spots["gProcessedSignal"].iloc[0])  # P0 printed twice
        controls = pd.DataFrame({"ProbeName": ["(+)E1A_r60_1", "NegativeControl"], "SystematicName": "ctrl",
                                 "gProcessedSignal": [50000.0, 2.0], "ControlType": [1, -1], "gIsWellAboveBG": [1, 0]})
        spots = pd.concat([spots, replicate, controls], ignore_index=True)
        spots["FeatureNum"] = range(1, len(spots) + 1)
        files[array] = root / f"{array}.txt"
        write_agilent_file(files[array], spots)

    config = write_config(
        root / "ARRAY.yaml", files["H1"], id="ARRAY", data_type="microarray",
        files={a: str(f) for a, f in files.items()}, samples={a: g for a, (g, _) in ARRAYS.items()},
        platform_annotation=str(annotation_file), fetal_sex={a: s for a, (_, s) in ARRAYS.items()},
        adjust_for_fetal_sex=False, composition_check={"name": "Decidua", "markers": ["GENE150", "GENE151"]})
    resources, results_root, site_data = write_resources(root / "resources"), root / "results", root / "site_data"
    step1_load_data.run(config, results_root)
    step2_differential_expression.run(config, results_root)
    step3_immune_annotation.run(config, results_root, resources)
    step4_export_site_data.run(config, results_root, site_data)
    return {"config": config, "results": results_root / "ARRAY", "results_root": results_root,
            "site_data": site_data, "resources": resources}


def test_reads_agilent_files(microarray_run):
    signal = read_tsv(microarray_run["results"] / "intensities.tsv")
    assert signal.shape == (202, 6)                       # control spots dropped, replicate spots averaged
    assert "NegativeControl" not in signal.index
    probes = read_tsv(microarray_run["results"] / "probes.tsv")
    assert probes.loc["PY", "chromosome"] == "Y" and probes.loc["P0", "gene_symbol"] == "GENE0"
    samples = read_tsv(microarray_run["results"] / "samples.tsv")
    assert samples["fetal_sex"].to_dict() == {a: s for a, (_, s) in ARRAYS.items()}


def test_microarray_statistics(microarray_run):
    table = read_tsv(microarray_run["results"] / "de_all_genes.tsv")
    info = json.loads((microarray_run["results"] / "step2_info.json").read_text())
    assert info["design"] == "~ condition" and info["measure"] == "intensity"
    assert "DDX3Y" not in table.index and info["chr_y_genes_excluded"] == 1   # sex not adjusted -> chrY left out
    assert table.loc["GENE5", "probe_id"] == "P5"                               # brighter of GENE5's two probes
    hits = table[(table["pvalue"] < 0.05) & (table["log2FoldChange"].abs() >= 1)]
    planted_up, planted_down = {f"GENE{i}" for i in range(10)}, {f"GENE{i}" for i in range(10, 20)}
    assert len(planted_up & set(hits.index)) >= 9 and len(planted_down & set(hits.index)) >= 9
    assert (table.loc[sorted(planted_up), "log2FoldChange"] > 1).all()
    assert len(set(hits.index) - planted_up - planted_down) <= 10
    estimates = read_tsv(microarray_run["results"] / "group_estimates.tsv").loc[table.index]
    assert np.allclose(np.log2(estimates["ICP"] / estimates["Healthy"]), table["log2FoldChange"])


def test_site_data_passes_the_cross_check(microarray_run):
    data = json.loads((microarray_run["site_data"] / "ARRAY.json").read_text())
    assert data["data_type"] == "microarray" and data["measure"] == "intensity" and data["chr_y_genes_excluded"] == 1
    assert all(s["library_size"] is None for s in data["samples"])
    assert step5_check_site_data.run(microarray_run["config"], microarray_run["results_root"],
                                     microarray_run["site_data"], microarray_run["resources"]) == 0


def test_quantile_normalization_equalizes_arrays():
    raw = pd.DataFrame({"a": [1.0, 5.0, 3.0], "b": [2.0, 10.0, 6.0]})
    normalized = quantile_normalize(raw)
    assert sorted(normalized["a"]) == sorted(normalized["b"]) == [1.5, 4.5, 7.5]
    assert normalized["a"].tolist() == [1.5, 7.5, 4.5]    # each array keeps its own ordering


def test_prior_is_recovered_from_simulated_variances():
    rng = np.random.default_rng(3)
    d0, s0_2, dof = 4.0, 0.05, 4
    true_variances = s0_2 * d0 / rng.chisquare(d0, 50000)      # scaled inverse chi-square prior
    observed = true_variances * rng.chisquare(dof, 50000) / dof
    fitted_d0, fitted_s0_2 = fit_prior(observed, dof)
    assert abs(fitted_d0 - d0) / d0 < 0.15 and abs(fitted_s0_2 - s0_2) / s0_2 < 0.1
    assert abs(trigamma_inverse(float(special.polygamma(1, 2.5))) - 2.5) < 1e-6


def test_moderated_t_equals_ordinary_t_without_a_prior():
    """If every gene has the same variance, the prior is infinitely strong and every gene's
    variance becomes the common one; with genes of different variances it lies in between."""
    expression = pd.DataFrame([[1.0, 2.0, 3.0, 4.0], [2.0, 2.5, 3.0, 3.5]], columns=list("abcd"))
    design = pd.DataFrame({"intercept": 1.0, "group: B": [0.0, 0.0, 1.0, 1.0]}, index=list("abcd"))
    table, prior = moderated_t_test(expression, design, "group: B")
    assert table["log2FoldChange"].round(6).tolist() == [2.0, 1.0]
    assert prior["residual_df"] == 2 and table["pvalue"].between(0, 1).all()
