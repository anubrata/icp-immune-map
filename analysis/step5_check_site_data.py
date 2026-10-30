"""
Step 5 - Cross-check the website data against an independent recomputation.

Recomputes every number the dashboard shows, straight from the tables of steps 1-2 and the
reference files, with its own code (it does not reuse the functions of steps 3-4), and compares
the results with site/data/<dataset>.json. Any difference beyond rounding is a FAIL and the script
exits with an error, so GitHub Actions stops before publishing a broken site.

Checked: which genes are significant; their fold change, p-value, FDR and mean count; per-patient
values; heatmap values; arrow-chart levels; ImmPort and atlas tags and the immune-related flag; the
immune-process bars; the summary-tile inputs; every immune cell-type score with its p-value and FDR;
and the flags for genes whose plain averages point the other way than the model.

Run:  uv run analysis/step5_check_site_data.py config/GSE282182.yaml
"""

from __future__ import annotations

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from common import RESOURCES_DIR, RESULTS_ROOT, SITE_DATA_DIR, load_config, read_tsv, results_dir, step_arguments

EPS = 1e-9  # room for floating-point noise on top of the rounding used in the export


class Checks:
    """Prints PASS/FAIL lines and counts the failures."""

    def __init__(self) -> None:
        self.passed = self.failed = 0

    def __call__(self, name: str, ok: bool, detail: str = "") -> None:
        print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))
        self.passed += bool(ok)
        self.failed += not ok


def relative(a: float, b: float) -> float:
    return abs(a - b) / max(abs(b), 1e-300)


def worst(values) -> float:
    return max(values, default=0.0)


def zscores(table: pd.DataFrame) -> pd.DataFrame:
    sd = table.std(axis=1).replace(0, np.nan)
    return table.sub(table.mean(axis=1), axis=0).div(sd, axis=0).fillna(0.0)


# ------------------------------------------------------------ gene symbols and reference sets

def official_symbols(de: pd.DataFrame, hgnc: pd.DataFrame) -> dict[str, str]:
    """Current official symbol of every tested gene: from its Ensembl ID if HGNC knows it, else its
    own symbol if still approved, else the symbol it was renamed to (if unambiguous)."""
    by_ensembl = {e: s for e, s in zip(hgnc["ensembl_gene_id"], hgnc["symbol"]) if e}
    approved = set(hgnc["symbol"])
    renamed = defaultdict(set)
    for symbol, previous in zip(hgnc["symbol"], hgnc["prev_symbol"]):
        for old in filter(None, (p.strip() for p in previous.split("|"))):
            renamed[old].add(symbol)
    result = {}
    for gene_id, symbol in zip(de.index, de["gene_symbol"].fillna("").astype(str)):
        if gene_id.split(".")[0] in by_ensembl:
            result[gene_id] = by_ensembl[gene_id.split(".")[0]]
        elif symbol in approved or len(renamed.get(symbol, ())) != 1:
            result[gene_id] = symbol
        else:
            result[gene_id] = next(iter(renamed[symbol]))
    return result


def read_immport(path: Path) -> dict[str, tuple[str, set[str]]]:
    """{list ID: (name, gene symbols)} from the ImmPort GMT file."""
    lists = {}
    for line in open(path):
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) >= 3:
            found = re.search(r"(GO:\d{7}|R-HSA-\d+)", fields[1])
            list_id = found.group(1) if found else fields[0]
            name, genes = lists.get(list_id, (fields[0], set()))
            lists[list_id] = (name, genes | {g.strip() for g in fields[2:] if g.strip()})
    return lists


def members(gene_symbols: dict[str, set[str]], sets: dict[str, set[str]]) -> dict[str, set[str]]:
    """{gene ID: IDs of the sets that contain any of the gene's symbols}."""
    sets_of_symbol = defaultdict(set)
    for set_id, genes in sets.items():
        for gene in genes:
            sets_of_symbol[gene].add(set_id)
    return {gene_id: set().union(*(sets_of_symbol.get(s, set()) for s in symbols))
            for gene_id, symbols in gene_symbols.items()}


# ------------------------------------------------------------------------------- checks

def run(config_path: str | Path, results_root: Path = RESULTS_ROOT, site_data_dir: Path = SITE_DATA_DIR,
        resources_dir: Path = RESOURCES_DIR) -> int:
    """Run every check; returns the number of failed checks."""
    cfg = load_config(config_path)
    folder, resources_dir = results_dir(cfg, results_root), Path(resources_dir)
    site = json.loads((Path(site_data_dir) / f"{cfg['id']}.json").read_text())
    de = read_tsv(folder / "de_all_genes.tsv")
    normalized, vst = read_tsv(folder / "normalized_counts.tsv"), read_tsv(folder / "vst.tsv")
    samples, estimates = read_tsv(folder / "samples.tsv"), read_tsv(folder / "group_estimates.tsv")
    info = json.loads((folder / "step2_info.json").read_text())
    reference, compare, cutoffs = cfg["reference_group"], cfg["compare_group"], cfg["cutoffs"]
    genes = site["genes"]
    by_id = {g["id"]: g for g in genes}
    order = [s["id"] for s in site["samples"]]
    check = Checks()
    print(f"Checking site/data/{cfg['id']}.json against an independent recomputation\n")

    # 1. Which genes are significant ---------------------------------------------------------
    p_column = "padj" if "padj" in cutoffs else "pvalue"
    significant = de[(de[p_column] < cutoffs[p_column]) & (de["log2FoldChange"].abs() >= cutoffs["abs_log2fc"])]
    check("Significant genes = DESeq2 results filtered by the config cutoffs",
          set(significant.index) == set(by_id), f"{len(significant)} expected, {len(by_id)} on the site")
    check("Genes are listed by p-value", [g["pvalue"] for g in genes] == sorted(g["pvalue"] for g in genes))
    rank = {reference: 0, compare: 1}
    groups_in_order = [rank.get(s["group"], 2) for s in site["samples"]]
    check("Patients: all of them, reference group first", sorted(order) == sorted(samples.index)
          and groups_in_order == sorted(groups_in_order))
    shared = [g for g in genes if g["id"] in significant.index]

    # 2. Gene statistics ---------------------------------------------------------------------
    d_lfc = worst(abs(g["log2fc"] - de.at[g["id"], "log2FoldChange"]) for g in shared)
    d_p = worst(relative(g["pvalue"], de.at[g["id"], "pvalue"]) for g in shared)
    d_padj = worst(relative(g["padj"], de.at[g["id"], "padj"]) for g in shared if g["padj"] is not None)
    padj_gaps_match = all((g["padj"] is None) == pd.isna(de.at[g["id"], "padj"]) for g in shared)
    d_mean = worst(abs(g["base_mean"] - de.at[g["id"], "baseMean"]) for g in shared)
    check("Fold change, p-value, FDR and mean count (up to rounding)",
          d_lfc <= 5e-4 + EPS and d_p <= 5e-3 + EPS and d_padj <= 5e-3 + EPS and d_mean <= 0.05 + EPS and padj_gaps_match,
          f"largest differences: log2FC {d_lfc:.1g}, p {d_p:.1%}, FDR {d_padj:.1%}, mean count {d_mean:.2g}")

    # 3. Per-patient values (gene panel) and heatmap ------------------------------------------
    d_norm = worst(np.abs(np.array(g["normalized"]) - normalized.loc[g["id"], order].to_numpy()).max() for g in shared)
    check("Per-patient normalized counts", d_norm <= 0.05 + EPS, f"largest difference {d_norm:.2g}")
    z = zscores(vst.loc[[g["id"] for g in shared], order])
    d_z = worst(np.abs(np.array(g["z"]) - z.loc[g["id"]].to_numpy()).max() for g in shared)
    check("Heatmap z-scores", d_z <= 5e-3 + EPS, f"largest difference {d_z:.2g}")

    # 4. Arrow chart --------------------------------------------------------------------------
    d_arrow = worst(abs(np.log2(g["group_estimates"][compare] / g["group_estimates"][reference]) - g["log2fc"]) for g in genes)
    d_level = worst(relative(g["group_estimates"][grp], estimates.at[g["id"], grp]) for g in shared for grp in (reference, compare))
    check("Arrow chart: each arrow spans exactly the gene's fold change", d_arrow < 0.01, f"largest difference {d_arrow:.2g} log2")
    check("Arrow chart: levels match the model's group estimates", d_level <= 5e-4 + EPS, f"largest difference {d_level:.1%}")

    # 5. Plain averages and the "points the other way" flag ----------------------------------------
    by_group = {grp: samples.index[samples["group"] == grp] for grp in (reference, compare)}
    wrong_flag, d_plain = [], 0.0
    for g in shared:
        values = normalized.loc[g["id"]]
        plain = np.log2((values[by_group[compare]].mean() + 0.5) / (values[by_group[reference]].mean() + 0.5))
        d_plain = max(d_plain, abs(plain - g["plain_log2fc"]))
        if (plain * g["log2fc"] < 0) != g["raw_disagrees"]:
            wrong_flag.append(g["symbol"])
    check("Plain-average change and the ⚠ flag (plain averages vs model)", not wrong_flag and d_plain <= 5e-4 + EPS,
          f"{sum(g['raw_disagrees'] for g in genes)} flagged: "
          f"{[g['symbol'] for g in genes if g['raw_disagrees']]}; wrong: {wrong_flag}")

    # 6. Immune tags --------------------------------------------------------------------------
    hgnc = pd.read_csv(resources_dir / "hgnc_symbols.tsv", sep="\t", dtype=str, keep_default_na=False)
    official = official_symbols(de, hgnc)
    in_file = de["gene_symbol"].fillna("").astype(str)
    symbols = {i: {in_file[i], official[i]} for i in de.index}
    immport = read_immport(resources_dir / "immport_gene_lists.gmt")
    atlas = pd.read_csv(resources_dir / "immune_cell_markers.tsv", sep="\t", dtype=str, keep_default_na=False)
    atlas_sets = {pid: {m for m in markers.split(",") if m} for pid, markers in zip(atlas["population_id"], atlas["markers"])}
    in_lists = members(symbols, {i: genes_ for i, (_, genes_) in immport.items()})
    in_atlas = members(symbols, atlas_sets)
    cell_populations = set(atlas.loc[atlas["kind"] == "cell population", "population_id"])
    immune = {i: bool(in_lists[i]) or bool(in_atlas[i] & cell_populations) for i in de.index}
    check("ImmPort lists of each significant gene",
          all(set(g["immune_lists"]) == in_lists[g["id"]] for g in shared),
          f"mismatches: {[g['symbol'] for g in shared if set(g['immune_lists']) != in_lists[g['id']]][:5]}")
    check("Atlas cell-type markers of each significant gene",
          all(set(g["atlas_markers"]) == in_atlas[g["id"]] for g in shared),
          f"mismatches: {[g['symbol'] for g in shared if set(g['atlas_markers']) != in_atlas[g['id']]][:5]}")
    check("Immune-related flag (ImmPort list or atlas cell-population marker)",
          all(g["immune"] == immune[g["id"]] for g in shared))

    # 7. Immune-process bars -----------------------------------------------------------------
    expected_bars = defaultdict(lambda: [0, 0])
    for gene_id in significant.index:
        for list_id in in_lists[gene_id]:
            expected_bars[list_id][0 if significant.at[gene_id, "log2FoldChange"] > 0 else 1] += 1
    shown_bars = {l["id"]: [l["higher"], l["lower"]] for l in site["immune_lists"]}
    check("Immune-process bars: significant genes per list, higher and lower",
          shown_bars == dict(expected_bars), f"{len(expected_bars)} lists")

    # 8. Summary tiles ---------------------------------------------------------------------------
    check("Tiles: genes in file and genes tested", site["genes_in_file"] == info["genes_in_file"]
          and site["genes_tested"] == len(de), f"{site['genes_tested']:,} of {site['genes_in_file']:,}")
    check("Tiles: immune-related genes among those tested", site["immune_genes_tested"] == sum(immune.values()),
          f"{sum(immune.values()):,}")
    higher = sum(g["log2fc"] > 0 for g in genes)
    print(f"      (tiles show {len(genes)} significant genes: {higher} higher, {len(genes) - higher} lower; "
          f"{sum(g['immune'] for g in genes)} immune-related)")

    # 9. Immune cell-type scores ---------------------------------------------------------------
    lookup = {}  # marker symbol -> gene: exact symbol in the data file first, else current symbol
    ranked = de.sort_values("baseMean", ascending=False).index
    for gene_id in ranked:
        lookup.setdefault(in_file[gene_id], gene_id)
    for gene_id in ranked:
        lookup.setdefault(official[gene_id], gene_id)
    z_all = zscores(vst[order])
    design = [np.ones(len(order))]
    if "fetal_sex" in info["design"]:
        design.append(samples.loc[order, "fetal_sex"].eq("M").to_numpy(float))
    design += [samples.loc[order, "group"].eq(grp).to_numpy(float)
               for grp in dict.fromkeys(samples.loc[order, "group"]) if grp not in (reference, compare)]
    design.append(samples.loc[order, "group"].eq(compare).to_numpy(float))  # last column = the ICP effect
    X = np.column_stack(design)
    inverse = np.linalg.pinv(X.T @ X)
    recomputed = {}
    for pid, markers in atlas_sets.items():
        ids = sorted({lookup[m] for m in markers if m in lookup})  # each tested gene counts once
        if len(ids) < 2:
            continue
        score = z_all.loc[ids].mean().to_numpy()
        beta = inverse @ X.T @ score
        residual = score - X @ beta
        dof = len(score) - X.shape[1]
        se = np.sqrt(residual @ residual / dof * inverse[-1, -1])
        recomputed[pid] = (beta[-1], 2 * stats.t.sf(abs(beta[-1] / se), dof))
    cells = {c["id"]: c for c in site["cell_types"]}
    check("Cell types scored (at least 2 markers tested)", set(recomputed) == set(cells), f"{len(recomputed)} populations")
    common_ids = set(recomputed) & set(cells)
    d_effect = worst(abs(recomputed[k][0] - cells[k]["effect"]) for k in common_ids)
    d_cell_p = worst(relative(cells[k]["pvalue"], recomputed[k][1]) for k in common_ids)
    check("Cell-type score differences and p-values", d_effect <= 5e-4 + 1e-6 and d_cell_p <= 5e-3 + 1e-6,
          f"largest differences: {d_effect:.1g} SD, p {d_cell_p:.1%}")
    ranked_ids = sorted(common_ids, key=lambda k: recomputed[k][1])
    p_sorted = np.array([recomputed[k][1] for k in ranked_ids])
    fdr = np.minimum.accumulate((p_sorted * len(p_sorted) / np.arange(1, len(p_sorted) + 1))[::-1])[::-1].clip(max=1)
    d_fdr = worst(relative(cells[k]["fdr"], fdr[i]) for i, k in enumerate(ranked_ids))
    check("Cell-type FDR (Benjamini-Hochberg across all scored populations)", d_fdr <= 5e-3 + 1e-6, f"largest difference {d_fdr:.1%}")

    print(f"\n{check.passed} checks passed, {check.failed} failed")
    return check.failed


if __name__ == "__main__":
    args = step_arguments("Step 5: cross-check the website data against an independent recomputation.",
                          {"--site-data-dir": (SITE_DATA_DIR, "where the website data is (default: site/data/)"),
                           "--resources-dir": (RESOURCES_DIR, "folder with the reference files (default: resources/)")})
    sys.exit(1 if run(args.config, args.results_root, args.site_data_dir, args.resources_dir) else 0)
