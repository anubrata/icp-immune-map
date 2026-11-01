import numpy as np
import pandas as pd

from common import read_tsv
from conftest import NEW_SYMBOL, OLD_SYMBOL, RENAMED_ID, gene_id, write_resources
from step3_immune_annotation import annotate_immune, read_atlas, read_gmt, read_hgnc, symbol_lookup


def test_marker_lookup_prefers_the_exact_symbol():
    """MRC1 must use the row named MRC1, not the more expressed former duplicate MRC1L1 (now also MRC1)."""
    genes = pd.DataFrame({"gene_symbol": ["MRC1", "MRC1L1", "C4orf21"], "current_symbol": ["MRC1", "MRC1", "ZGRF1"],
                          "baseMean": [123.0, 446.0, 50.0]}, index=["exact", "former_duplicate", "renamed"])
    lookup = symbol_lookup(genes)
    assert lookup["MRC1"] == "exact"
    assert lookup["ZGRF1"] == "renamed"          # renamed genes are still found through the current symbol


def test_reads_gmt_lists_with_their_ids(tmp_path):
    lists = read_gmt(write_resources(tmp_path) / "immport_gene_lists.gmt")
    assert set(lists) == {"GO:0002218", "R-HSA-1059683"}
    assert lists["GO:0002218"]["source"] == "GO"
    assert lists["R-HSA-1059683"]["source"] == "Reactome"
    assert lists["R-HSA-1059683"]["name"] == "Interleukin-6 signaling"
    assert lists["GO:0002218"]["genes"] == {"GENE0", "GENE1", NEW_SYMBOL, "NOT_IN_DATA"}


def test_matches_by_symbol_and_by_current_symbol(tmp_path):
    resources = write_resources(tmp_path)
    atlas = read_atlas(resources / "immune_cell_markers.tsv")
    hgnc = read_hgnc(resources / "hgnc_symbols.tsv")
    table = pd.DataFrame({"gene_symbol": ["GENE0", OLD_SYMBOL, OLD_SYMBOL, "IFNB2", "GENE7"]},
                         index=["G0", f"{RENAMED_ID}.5", "no_ensembl_id", "other", "G7"])

    atlas_only = annotate_immune(table, atlas, hgnc)
    assert "immune_lists" not in atlas_only                       # ImmPort off by default
    assert atlas_only.loc["G0", "atlas_markers"] == "1.1"
    # renamed gene, matched through its Ensembl ID (version suffix ignored)...
    assert atlas_only.loc[f"{RENAMED_ID}.5", "current_symbol"] == NEW_SYMBOL
    assert atlas_only.loc[f"{RENAMED_ID}.5", "atlas_markers"] == "1.3"
    # ...or through HGNC's previous-symbol table when there is no Ensembl ID
    assert atlas_only.loc["no_ensembl_id", "current_symbol"] == NEW_SYMBOL
    assert atlas_only.loc["other", "current_symbol"] == "IL6"
    assert atlas_only["immune"].tolist() == [True, True, True, False, False]

    with_immport = annotate_immune(table, atlas, hgnc, read_gmt(resources / "immport_gene_lists.gmt"))
    assert with_immport.loc["G0", "immune_lists"] == "GO:0002218;R-HSA-1059683"
    assert with_immport.loc[f"{RENAMED_ID}.5", "immune_lists"] == "GO:0002218"
    assert with_immport.loc["G7", "immune_lists"] == ""


def test_atlas_markers_tag_genes(pipeline_run):
    annotated = read_tsv(pipeline_run["results"] / "de_all_genes_annotated.tsv")
    assert annotated.loc[gene_id(2), "atlas_markers"] == "1.1"
    assert annotated.loc[gene_id(2), "immune"]                       # marker of an atlas cell population
    assert annotated.loc[gene_id(5), "atlas_markers"] == "13.1"
    assert not annotated.loc[gene_id(5), "immune"]                   # a functional-state marker alone doesn't count
    assert annotated.loc[gene_id(305), "atlas_markers"] == "1.10"    # read back as text, not the number 1.1
    assert annotated.loc[RENAMED_ID, "gene_symbol"] == OLD_SYMBOL
    assert annotated.loc[RENAMED_ID, "current_symbol"] == NEW_SYMBOL
    assert annotated.loc[RENAMED_ID, "atlas_markers"] == "1.3" and annotated.loc[RENAMED_ID, "immune"]
    assert not annotated.loc[gene_id(1), "immune"]                   # only in an ImmPort list, and ImmPort is off
    assert not (pipeline_run["results"] / "immune_lists.tsv").exists()


def test_cell_type_scores_find_the_planted_population(pipeline_run):
    scores = read_tsv(pipeline_run["results"] / "cell_type_scores.tsv")
    assert scores.index.tolist() == ["1.1", "1.2", "1.3", "1.10", "13.1"]

    planted = scores.loc["1.1"]
    assert planted["n_found"] == 4 and planted["markers_higher"] == 4
    assert planted["effect"] > 1 and planted["ci_low"] > 0 and planted["pvalue"] < 0.01
    assert abs(scores.loc["1.2", "effect"]) < planted["effect"]
    assert scores.loc["1.3", "n_found"] == 2                          # the renamed marker is found
    assert np.isnan(scores.loc["1.10", "pvalue"])                     # only 1 marker tested: not scored
    assert -1 <= planted["composition_r"] <= 1 and not np.isnan(planted["pvalue_composition_adjusted"])

    patients = read_tsv(pipeline_run["results"] / "cell_type_patient_scores.tsv")
    assert patients.index.tolist() == ["1.1", "1.2", "1.3", "13.1"] and patients.shape[1] == 10
    assert read_tsv(pipeline_run["results"] / "composition_score.tsv").shape == (10, 1)
