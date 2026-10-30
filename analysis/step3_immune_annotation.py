"""
Step 3 - Immune annotation: which genes are immune-related, and which immune cell types change?

Two references:
  * ImmPort immune gene lists: ~150 lists from GO "immune system process" and Reactome "Immune
    System", one per immune process (e.g. "Activation of innate immune response").
  * The Human Immune Cell Atlas (resources/human_immune_cell_atlas.docx, made into
    resources/immune_cell_markers.tsv by step 0): marker genes for ~200 immune cell populations
    and states (e.g. "NK cells: NKG7, GNLY, KLRD1, NCR1, PRF1").

Part A - tag genes. Each tested gene gets the ImmPort lists it belongs to and the atlas
populations it is a marker of. A gene counts as immune-related if it is in an ImmPort list or
marks an atlas cell population. (The atlas's cross-lineage "functional state" signatures, such
as proliferation or hypoxia, don't count: the atlas notes they are not specific to immune cells.)
Symbols change over time (e.g. C4orf21 is now ZGRF1), so each gene is matched both by its
symbol in the data file and by its current official symbol, looked up from its Ensembl ID.

Part B - score immune cell types. For every atlas population with at least 2 marker genes among
the tested genes, each patient gets a score: the average of those markers' expression, after
scaling each gene to mean 0 and standard deviation (SD) 1 across patients (variance-stabilized
counts from step 2). Scores are compared between groups with the same adjustment as step 2
(score ~ fetal sex + group), giving the difference in SD, a 95% confidence interval, a p-value
and an FDR across all scored populations. With `composition_check` in the config, each score is
also compared with a tissue-composition signature (in placenta: the amount of decidua, maternal
tissue rich in NK cells, which varies between samples and can mimic a change in immune cells).

Inputs:  results/<dataset>/de_all_genes.tsv, vst.tsv, samples.tsv, step2_info.json  (steps 1-2)
         resources/immport_gene_lists.gmt, hgnc_symbols.tsv, immune_cell_markers.tsv (step 0)
Outputs: results/<dataset>/de_all_genes_annotated.tsv     step 2's table + current_symbol,
                                                         immune_lists, atlas_markers, immune
         results/<dataset>/immune_lists.tsv              one row per ImmPort list
         results/<dataset>/cell_type_scores.tsv          one row per atlas population: markers found,
                                                         ICP difference, 95% CI, p, FDR, composition check
         results/<dataset>/cell_type_patient_scores.tsv  each scored population's score per patient
         results/<dataset>/composition_score.tsv         the composition signature per patient (if set)

Run:     uv run analysis/step3_immune_annotation.py config/GSE282182.yaml
"""

from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from common import (RESOURCES_DIR, RESULTS_ROOT, describe_cutoffs, is_significant, load_config, read_tsv,
                    results_dir, row_zscores, show_path, step_arguments, write_tsv)

LIST_ID = re.compile(r"(GO:\d{7}|R-HSA-\d+)")
MIN_MARKERS = 2  # a population is scored when at least this many of its markers were tested


# ------------------------------------------------------------------ reference files

def read_gmt(path: Path) -> dict[str, dict]:
    """Read a GMT file (one list per line: name, link, then gene symbols).

    Returns {list ID: {"name", "source", "url", "genes"}}, where the ID is the GO or Reactome
    identifier found in the link (e.g. GO:0002218 or R-HSA-879415).
    """
    lists: dict[str, dict] = {}
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                continue
            name, url, symbols = fields[0], fields[1], fields[2:]
            match = LIST_ID.search(url)
            list_id = match.group(1) if match else name
            source = "GO" if list_id.startswith("GO:") else "Reactome" if list_id.startswith("R-HSA-") else ""
            entry = lists.setdefault(list_id, {"name": name, "source": source, "url": url, "genes": set()})
            entry["genes"].update(s.strip() for s in symbols if s.strip())
    return lists


def read_atlas(path: Path) -> pd.DataFrame:
    """The atlas marker table, one row per population, with its markers as a list in `genes`."""
    atlas = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).set_index("population_id")
    atlas["genes"] = atlas["markers"].map(lambda text: [g for g in text.split(",") if g])
    return atlas


def population_order(population_id: str) -> tuple:
    """Sort atlas IDs like '1.2' < '1.10' < '2.1'."""
    return tuple(int(part) for part in population_id.split("."))


@dataclass
class HgncLookup:
    """Lookups from the HGNC table to today's official gene symbol."""
    by_ensembl: dict[str, str]   # Ensembl gene ID -> current symbol
    by_previous: dict[str, str]  # previous symbol -> current symbol (only when unambiguous)
    approved: set[str]           # all current symbols


def read_hgnc(path: Path) -> HgncLookup:
    hgnc = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    by_ensembl = {e: s for e, s in zip(hgnc["ensembl_gene_id"], hgnc["symbol"]) if e}
    renamed_to = defaultdict(set)
    for symbol, previous in zip(hgnc["symbol"], hgnc["prev_symbol"]):
        for old in filter(None, (p.strip() for p in previous.split("|"))):
            renamed_to[old].add(symbol)
    # A previous symbol later split between several genes is ambiguous, so it is not used.
    by_previous = {old: next(iter(new)) for old, new in renamed_to.items() if len(new) == 1}
    return HgncLookup(by_ensembl, by_previous, set(hgnc["symbol"]))


def current_symbol(gene_id: str, symbol: str, hgnc: HgncLookup) -> str:
    """Today's official symbol: via the Ensembl ID if HGNC knows it, otherwise via the
    previous-symbol table (unless the symbol is still an approved one), otherwise unchanged."""
    ensembl_id = gene_id.split(".")[0]  # drop a version suffix such as ENSG00000136488.5
    if ensembl_id in hgnc.by_ensembl:
        return hgnc.by_ensembl[ensembl_id]
    if symbol in hgnc.approved:
        return symbol
    return hgnc.by_previous.get(symbol, symbol)


# ---------------------------------------------------------------- Part A: gene tags

def memberships(symbol_pairs: list[tuple[str, str]], sets: dict[str, set[str]], order=str) -> list[str]:
    """For each gene (symbol in the file, current symbol): ';'-joined IDs of the sets containing it."""
    sets_by_symbol = defaultdict(set)
    for set_id, genes in sets.items():
        for gene in genes:
            sets_by_symbol[gene].add(set_id)
    return [";".join(sorted(sets_by_symbol.get(old, set()) | sets_by_symbol.get(new, set()), key=order))
            for old, new in symbol_pairs]


def annotate_immune(table: pd.DataFrame, lists: dict[str, dict], atlas: pd.DataFrame, hgnc: HgncLookup) -> pd.DataFrame:
    """Add current_symbol, immune_lists (ImmPort IDs), atlas_markers (atlas IDs) and immune (True/False)."""
    symbols = table["gene_symbol"].fillna("").astype(str)
    current = [current_symbol(gene_id, symbol, hgnc) for gene_id, symbol in zip(table.index, symbols)]
    pairs = list(zip(symbols, current))

    annotated = table.copy()
    annotated.insert(annotated.columns.get_loc("gene_symbol") + 1, "current_symbol", current)
    annotated["immune_lists"] = memberships(pairs, {list_id: info["genes"] for list_id, info in lists.items()})
    annotated["atlas_markers"] = memberships(pairs, {pid: set(genes) for pid, genes in atlas["genes"].items()},
                                             order=population_order)
    cell_populations = set(atlas.index[atlas["kind"] == "cell population"])
    annotated["immune"] = [bool(in_lists) or any(pid in cell_populations for pid in markers.split(";"))
                           for in_lists, markers in zip(annotated["immune_lists"], annotated["atlas_markers"])]
    return annotated


def lists_table(lists: dict[str, dict], annotated: pd.DataFrame) -> pd.DataFrame:
    """One row per ImmPort list, with its size and how many of its genes were tested here."""
    tested = Counter(list_id for m in annotated["immune_lists"] if m for list_id in m.split(";"))
    rows = [{"list_id": list_id, "name": info["name"], "source": info["source"], "url": info["url"],
             "genes_in_list": len(info["genes"]), "genes_tested": tested.get(list_id, 0)}
            for list_id, info in lists.items()]
    return pd.DataFrame(rows).set_index("list_id").sort_values("name")


# ------------------------------------------------------- Part B: cell-type scores

def symbol_lookup(annotated: pd.DataFrame) -> dict[str, str]:
    """Marker symbol -> the tested gene used for it in the cell-type scores.

    A gene whose symbol in the data file matches exactly wins; a gene that only matches through
    its current official symbol (renamed, or a former duplicate such as MRC1L1 -> MRC1) is the
    fallback. If several genes qualify at the same level, the most highly expressed one is used.
    """
    ranked = annotated.sort_values("baseMean", ascending=False)
    lookup: dict[str, str] = {}
    for gene_id, original in zip(ranked.index, ranked["gene_symbol"]):
        lookup.setdefault(original, gene_id)
    for gene_id, current in zip(ranked.index, ranked["current_symbol"]):
        lookup.setdefault(current, gene_id)
    return lookup


def group_design(samples: pd.DataFrame, reference: str, adjust_for_sex: bool) -> pd.DataFrame:
    """Columns of the model 'score ~ fetal sex + group', the same adjustment as step 2."""
    design = pd.DataFrame({"intercept": 1.0}, index=samples.index)
    if adjust_for_sex:
        design["male fetus"] = (samples["fetal_sex"] == "M").astype(float)
    for group in samples["group"].unique():
        if group != reference:
            design[f"group: {group}"] = (samples["group"] == group).astype(float)
    return design


def ols_term(y: np.ndarray, design: pd.DataFrame, term: str) -> dict:
    """Least-squares fit of y on the design; the estimate of one term with its 95% CI and p-value."""
    X = design.to_numpy(dtype=float)
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    residual = y - X @ beta
    dof = X.shape[0] - X.shape[1]
    j = list(design.columns).index(term)
    se = np.sqrt(residual @ residual / dof * np.linalg.inv(X.T @ X)[j, j])
    half_width = stats.t.ppf(0.975, dof) * se
    return {"effect": beta[j], "ci_low": beta[j] - half_width, "ci_high": beta[j] + half_width,
            "pvalue": 2 * stats.t.sf(abs(beta[j] / se), dof)}


def score_cell_types(annotated: pd.DataFrame, vst: pd.DataFrame, samples: pd.DataFrame, atlas: pd.DataFrame,
                     reference: str, compare: str, adjust_for_sex: bool, composition: dict | None = None):
    """Score every atlas population per patient and test the difference between groups.

    Returns (summary table, patient scores [populations x samples], composition score or None).
    """
    z = row_zscores(vst[samples.index])
    lookup = symbol_lookup(annotated)
    design = group_design(samples, reference, adjust_for_sex)
    term = f"group: {compare}"

    composition_score = None
    if composition:
        ids = [lookup[s] for s in composition["markers"] if s in lookup]
        composition_score = z.loc[ids].mean() if ids else None

    rows, scores = [], {}
    for population_id, population in atlas.iterrows():
        found = {}  # tested gene ID -> marker symbol as written in the atlas
        for symbol in population["genes"]:
            gene_id = lookup.get(symbol)
            if gene_id and gene_id not in found:
                found[gene_id] = symbol
        row = {"population_id": population_id, "section": population["section"],
               "population": population["population"], "kind": population["kind"],
               "markers": ",".join(population["genes"]), "markers_found": ",".join(found.values()),
               "marker_gene_ids": ";".join(found), "n_markers": len(population["genes"]), "n_found": len(found)}
        if len(found) >= MIN_MARKERS:
            score = z.loc[list(found)].mean()
            log2fc = annotated.loc[list(found), "log2FoldChange"]
            row |= {"markers_lower": int((log2fc < 0).sum()), "markers_higher": int((log2fc > 0).sum()),
                    "median_log2fc": float(log2fc.median()), **ols_term(score.to_numpy(), design, term)}
            if composition_score is not None:
                adjusted = ols_term(score.to_numpy(), design.assign(composition=composition_score.to_numpy()), term)
                row |= {"composition_r": float(np.corrcoef(score, composition_score)[0, 1]),
                        "effect_composition_adjusted": adjusted["effect"],
                        "pvalue_composition_adjusted": adjusted["pvalue"]}
            scores[population_id] = score
        rows.append(row)

    columns = ["section", "population", "kind", "markers", "markers_found", "marker_gene_ids", "n_markers", "n_found",
               "markers_lower", "markers_higher", "median_log2fc", "effect", "ci_low", "ci_high", "pvalue", "fdr",
               "composition_r", "effect_composition_adjusted", "pvalue_composition_adjusted"]
    summary = pd.DataFrame(rows).set_index("population_id").reindex(columns=columns)
    scored = summary["pvalue"].notna()
    if scored.any():
        summary.loc[scored, "fdr"] = stats.false_discovery_control(summary.loc[scored, "pvalue"].astype(float))
    patient_scores = pd.DataFrame(scores).T.reindex(columns=samples.index)
    patient_scores.index.name = "population_id"
    return summary, patient_scores, composition_score


# ------------------------------------------------------------------------------ run

def run(config_path: str | Path, results_root: Path = RESULTS_ROOT, resources_dir: Path = RESOURCES_DIR) -> None:
    cfg = load_config(config_path)
    folder = results_dir(cfg, results_root)
    resources_dir = Path(resources_dir)
    table = read_tsv(folder / "de_all_genes.tsv")
    vst = read_tsv(folder / "vst.tsv")
    samples = read_tsv(folder / "samples.tsv")
    info = json.loads((folder / "step2_info.json").read_text())
    lists = read_gmt(resources_dir / "immport_gene_lists.gmt")
    hgnc = read_hgnc(resources_dir / "hgnc_symbols.tsv")
    atlas = read_atlas(resources_dir / "immune_cell_markers.tsv")

    # Part A - tag genes
    annotated = annotate_immune(table, lists, atlas, hgnc)
    write_tsv(annotated, folder / "de_all_genes_annotated.tsv")
    write_tsv(lists_table(lists, annotated), folder / "immune_lists.tsv")

    significant = is_significant(annotated, cfg["cutoffs"])
    in_lists = annotated["immune_lists"] != ""
    renamed = int((annotated["current_symbol"] != annotated["gene_symbol"]).sum())
    print(f"References: {len(lists)} ImmPort immune gene lists; {len(atlas)} atlas populations/states")
    print(f"Tested genes whose official symbol changed since the data file was made: {renamed:,}")
    print(f"Immune-related: {annotated['immune'].sum():,} of {len(annotated):,} tested genes "
          f"({in_lists.sum():,} in ImmPort lists, {(annotated['immune'] & ~in_lists).sum():,} more from atlas markers)")
    print(f"  and {annotated.loc[significant, 'immune'].sum()} of the {significant.sum()} genes passing "
          f"{describe_cutoffs(cfg['cutoffs'])}")

    # Part B - score immune cell types
    summary, patient_scores, composition = score_cell_types(
        annotated, vst, samples, atlas, cfg["reference_group"], cfg["compare_group"],
        adjust_for_sex="fetal_sex" in info["design"], composition=cfg["composition_check"])
    write_tsv(summary, folder / "cell_type_scores.tsv")
    write_tsv(patient_scores.round(4), folder / "cell_type_patient_scores.tsv")
    composition_path = folder / "composition_score.tsv"
    if composition is not None:
        write_tsv(composition.rename("score").to_frame().round(4), composition_path)
    elif composition_path.exists():
        composition_path.unlink()

    scored = summary.dropna(subset=["pvalue"]).sort_values("pvalue")
    print(f"\nCell-type scores: {len(scored)} of {len(summary)} atlas populations have >= {MIN_MARKERS} markers "
          f"among the tested genes. Strongest differences, {cfg['compare_group']} vs {cfg['reference_group']} "
          f"(score ~ {'fetal sex + ' if 'fetal_sex' in info['design'] else ''}group, in SD):")
    top = scored.head(10).assign(
        markers=lambda t: t["n_found"].astype(int).astype(str) + "/" + t["n_markers"].astype(int).astype(str),
        difference=lambda t: t["effect"].map("{:+.2f}".format),
        ci95=lambda t: "[" + t["ci_low"].map("{:+.2f}".format) + ", " + t["ci_high"].map("{:+.2f}".format) + "]",
        p=lambda t: t["pvalue"].map("{:.3g}".format), FDR=lambda t: t["fdr"].map("{:.2g}".format))
    shown = ["population", "markers", "difference", "ci95", "p", "FDR"]
    if composition is not None:
        name = cfg["composition_check"]["name"]
        top[f"r with {name}"] = top["composition_r"].map("{:+.2f}".format)
        top[f"p adjusted for {name}"] = top["pvalue_composition_adjusted"].map("{:.3g}".format)
        shown += [f"r with {name}", f"p adjusted for {name}"]
    print(top[shown].to_string())
    print(f"\nWrote {show_path(folder)}/de_all_genes_annotated.tsv, immune_lists.tsv, cell_type_scores.tsv, "
          "cell_type_patient_scores.tsv" + (", composition_score.tsv" if composition is not None else ""))


if __name__ == "__main__":
    args = step_arguments("Step 3: immune gene tags (ImmPort + atlas) and immune cell-type scores.",
                          {"--resources-dir": (RESOURCES_DIR, "folder with the step 0 files (default: resources/)")})
    run(args.config, args.results_root, args.resources_dir)
