"""Synthetic test data. The tests never touch the real dataset (apart from checking its config)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

import step1_load_data
import step2_differential_expression
import step3_immune_annotation
import step4_export_site_data

HEALTHY = ["H1", "H2", "H3", "H4", "H5"]
ICP = ["I1", "I2", "I3", "I4", "I5"]
SAMPLES = HEALTHY + ICP
MALE = {"H4", "I1", "I4", "I5"}  # like GSE282182: 1 of 5 healthy and 3 of 5 ICP samples are male

# A gene that was renamed: the data file uses the old symbol, ImmPort and HGNC the new one.
RENAMED_ID, OLD_SYMBOL, NEW_SYMBOL = "ENSG00000138658", "C4orf21", "ZGRF1"
SEXFLIP_ID = "ENSG00000000003"  # changes in opposite directions in female and male placentas


def gene_id(i: int) -> str:
    """ID of simulated gene number i (its symbol is GENE<i>)."""
    return f"ENSG9{i:010d}"


def simulate_counts(n_genes: int = 400, n_changed: int = 40, fold: float = 4.0, seed: int = 7):
    """Negative-binomial read counts, genes x samples. The first n_changed genes differ in ICP:
    the first half `fold` times higher, the second half `fold` times lower.
    Returns (counts, ids of higher genes, ids of lower genes)."""
    rng = np.random.default_rng(seed)
    mean = rng.lognormal(mean=6.0, sigma=0.8, size=n_genes)  # a typical gene gets a few hundred reads
    effect = np.ones(n_genes)
    effect[: n_changed // 2] = fold
    effect[n_changed // 2 : n_changed] = 1 / fold
    in_icp = np.array([s in ICP for s in SAMPLES])
    mu = mean[:, None] * np.where(in_icp, effect[:, None], 1.0)
    dispersion = 0.05
    counts = rng.negative_binomial(n=1 / dispersion, p=1 / (1 + mu * dispersion))
    ids = [gene_id(i) for i in range(n_genes)]
    return pd.DataFrame(counts, index=ids, columns=SAMPLES), ids[: n_changed // 2], ids[n_changed // 2 : n_changed]


def write_gse_like_file(path: Path, counts: pd.DataFrame) -> None:
    """Write counts in the GSE282182 layout: gene_id, sample counts, then annotation columns.
    Also adds sex-chromosome genes, an all-zero gene, a duplicated symbol and a renamed gene."""
    header = ["gene_id", *SAMPLES, "gene_name", "gene_chr", "gene_start", "gene_end", "gene_strand",
              "gene_length", "gene_biotype", "gene_description", "tf_family"]
    lines = ["\t".join(header)]

    def add(gene_id, values, name, chrom, description, biotype="protein_coding"):
        lines.append("\t".join([gene_id, *map(str, values), name, chrom, "1", "100", "+", "100",
                                biotype, description, "-"]))

    for i, (gene_id, row) in enumerate(counts.iterrows()):
        add(gene_id, row.tolist(), f"GENE{i}", str(1 + i % 22), f"made-up gene {i} [Source:HGNC Symbol;Acc:HGNC:{i}]")
    male = [s in MALE for s in SAMPLES]
    add("ENSG00000129824", [300 if m else 0 for m in male], "RPS4Y1", "Y", "ribosomal protein S4 Y-linked 1")
    add("ENSG00000067048", [250 if m else 0 for m in male], "DDX3Y", "Y", "DEAD-box helicase 3 Y-linked")
    add("ENSG00000229807", [3 if m else 600 for m in male], "XIST", "X", "X inactive specific transcript", "lincRNA")
    add("ENSG00000000001", [0] * 10, "SILENT", "2", "NA")
    add("ENSG00000000002", counts.iloc[len(counts) // 2].tolist(), "GENE1", "3", '"duplicated symbol, quoted"')
    add(RENAMED_ID, counts.iloc[1].tolist(), OLD_SYMBOL, "4", "zinc finger GRF-type containing 1")
    # Like PTX3 in GSE282182: lower in ICP in female placentas, higher in male placentas. Since most
    # ICP placentas are male, the plain averages say "higher" while the sex-adjusted model says "lower".
    sexflip = {"H1": 600, "H2": 650, "H3": 700, "H4": 150, "H5": 580, "I1": 1200, "I2": 30, "I3": 35, "I4": 1100, "I5": 1250}
    add(SEXFLIP_ID, [sexflip[s] for s in SAMPLES], "SEXFLIP", "5", "sex-dependent change")
    path.write_text("\n".join(lines) + "\n")


def write_config(path: Path, data_file: Path, **overrides) -> Path:
    cfg = {
        "id": "SYNTH",
        "title": "Synthetic test data",
        "tissue": "Placenta",
        "geo_accession": None,
        "file": str(data_file),
        "columns": {"gene_id": "gene_id", "gene_symbol": "gene_name", "description": "gene_description",
                    "chromosome": "gene_chr", "biotype": "gene_biotype"},
        "samples": {s: "Healthy" if s in HEALTHY else "ICP" for s in SAMPLES},
        "reference_group": "Healthy",
        "compare_group": "ICP",
        "fetal_sex": "infer",
        "composition_check": {"name": "Decidua", "markers": ["GENE310", "GENE311"]},
        "cutoffs": {"pvalue": 0.05, "abs_log2fc": 1.0},
    }
    cfg.update(overrides)
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return path


def write_resources(folder: Path) -> Path:
    """Tiny stand-ins for resources/immport_gene_lists.gmt and resources/hgnc_symbols.tsv."""
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "immport_gene_lists.gmt").write_text(
        "Activation of innate immune response\thttps://amigo.geneontology.org/amigo/term/GO:0002218"
        f"\tGENE0\tGENE1\t{NEW_SYMBOL}\tNOT_IN_DATA\n"
        "Interleukin-6 signaling\thttps://reactome.org/PathwayBrowser/#/R-HSA-1059683\tGENE0\tGENE25\n"
    )
    pd.DataFrame({"symbol": [NEW_SYMBOL, "IL6"], "ensembl_gene_id": [RENAMED_ID, "ENSG00000136244"],
                  "prev_symbol": [OLD_SYMBOL, "IFNB2"]}).to_csv(folder / "hgnc_symbols.tsv", sep="\t", index=False)
    # Stand-in for the atlas marker table: GENE0-19 are planted "higher in ICP" genes, GENE300+ are unchanged.
    major, states = "1. Major lineages", "13. Common immune-cell functional states across lineages"
    pd.DataFrame([
        ["1.1", major, "Planted-up cells", "cell population", "GENE0,GENE2,GENE3,GENE4"],
        ["1.2", major, "Unchanged cells", "cell population", "GENE300,GENE301,GENE302,GENE303"],
        ["1.3", major, "Renamed-marker cells", "cell population", f"{NEW_SYMBOL},GENE304"],
        ["1.10", major, "One-marker cells", "cell population", "GENE305,NOT_IN_DATA"],
        ["13.1", states, "Cytotoxicity", "functional state", "GENE5,GENE6"],
    ], columns=["population_id", "section", "population", "kind", "markers"]).assign(
        low_markers="", source_text="").to_csv(folder / "immune_cell_markers.tsv", sep="\t", index=False)
    return folder


@pytest.fixture(scope="session")
def pipeline_run(tmp_path_factory):
    """Run steps 1-4 once on synthetic data; tests inspect the outputs."""
    root = tmp_path_factory.mktemp("pipeline")
    counts, higher, lower = simulate_counts()
    data_file = root / "synthetic_gene_count.txt"
    write_gse_like_file(data_file, counts)
    config = write_config(root / "SYNTH.yaml", data_file)
    resources = write_resources(root / "resources")
    results_root, site_data = root / "results", root / "site_data"

    step1_load_data.run(config, results_root)
    step2_differential_expression.run(config, results_root)
    step3_immune_annotation.run(config, results_root, resources)
    step4_export_site_data.run(config, results_root, site_data)
    return {"root": root, "config": config, "data_file": data_file, "results": results_root / "SYNTH",
            "results_root": results_root, "site_data": site_data, "resources": resources,
            "higher": higher, "lower": lower}
