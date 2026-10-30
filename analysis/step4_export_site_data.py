"""
Step 4 - Export the significant genes for the website.

A gene is significant when it passes both cutoffs in the config (here: unadjusted p < 0.05
and |log2 fold change| >= 1, i.e. at least 2-fold higher or lower in ICP). Only those genes
are written for the website.

Immune cell-type scores (step 3) are exported for every scored atlas population, since they
are computed from all tested genes, not only the significant ones.

Inputs:  results/<dataset>/de_all_genes_annotated.tsv, immune_lists.tsv, cell_type_scores.tsv,
         cell_type_patient_scores.tsv, composition_score.tsv (if any)       (from step 3)
         results/<dataset>/group_estimates.tsv, normalized_counts.tsv, vst.tsv, pca.tsv,
         step2_info.json                                                    (from step 2)
         results/<dataset>/samples.tsv                                      (from step 1)
Outputs: site/data/<dataset>.json                   everything the dashboard shows for this dataset
         site/data/<dataset>_significant_genes.csv  the same genes as a spreadsheet (downloadable from the site)
         site/data/<dataset>_cell_type_scores.csv   the immune cell-type scores as a spreadsheet
         site/data/datasets.json                    list of all exported datasets, for the site's dataset picker

Run:     uv run analysis/step4_export_site_data.py config/GSE282182.yaml
"""

from __future__ import annotations

import json
import math
import re
from collections import Counter
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from common import (RESULTS_ROOT, SITE_DATA_DIR, describe_cutoffs, is_significant, load_config, read_tsv,
                    results_dir, row_zscores, show_path, step_arguments)


def number(value, digits: int = 3, significant_digits: bool = False) -> float | None:
    """Round for the website; missing values become null in the JSON."""
    if value is None or pd.isna(value) or math.isinf(value):
        return None
    return float(f"{value:.{digits}g}") if significant_digits else round(float(value), digits)


def text(value) -> str | None:
    """Annotation value as text (chromosome 17 read back as the number 17.0 becomes '17')."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    return str(value).strip() or None


def clean_description(value) -> str | None:
    """Drop the '[Source:HGNC Symbol;Acc:HGNC:2440]' tail that Ensembl adds to gene descriptions."""
    description = text(value)
    return re.sub(r"\s*\[Source:[^\]]*\]$", "", description) if description else None


def ordered_samples(samples: pd.DataFrame, reference: str, compare: str) -> pd.DataFrame:
    """Reference group first, then the compare group, then any others (file order within a group)."""
    rank = samples["group"].map({reference: 0, compare: 1}).fillna(2)
    return samples.iloc[np.argsort(rank.to_numpy(), kind="stable")]


def group_means(normalized: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    """Plain average of each gene's normalized counts in each group (genes x groups)."""
    return normalized[samples.index].T.groupby(samples["group"], sort=False).mean().T


def plain_changes(normalized: pd.DataFrame, samples: pd.DataFrame, reference: str, compare: str) -> pd.DataFrame:
    """log2 ratio of plain averages (compare vs reference) for every gene: over all patients ("all")
    and, when fetal sex is known, within female ("F") and within male ("M") placentas.
    0.5 is added to both averages so a group average of 0 does not divide by zero."""
    def log2_ratio(keep: pd.Series) -> pd.Series:
        ref = samples.index[(samples["group"] == reference) & keep]
        cmp = samples.index[(samples["group"] == compare) & keep]
        if len(ref) == 0 or len(cmp) == 0:
            return pd.Series(np.nan, index=normalized.index)
        return np.log2((normalized[cmp].mean(axis=1) + 0.5) / (normalized[ref].mean(axis=1) + 0.5))

    changes = pd.DataFrame({"all": log2_ratio(pd.Series(True, index=samples.index))})
    if "fetal_sex" in samples:
        for sex in ("F", "M"):
            changes[sex] = log2_ratio(samples["fetal_sex"] == sex)
    return changes


def gene_records(significant: pd.DataFrame, normalized: pd.DataFrame, zscores: pd.DataFrame,
                 estimates: pd.DataFrame, plain_means: pd.DataFrame, plain: pd.DataFrame,
                 sample_order: list[str]) -> list[dict]:
    records = []
    for gene_id, gene in significant.iterrows():
        records.append({
            "id": gene_id,
            "symbol": gene["gene_symbol"],
            "current_symbol": gene["current_symbol"] if gene["current_symbol"] != gene["gene_symbol"] else None,
            "description": clean_description(gene.get("description")),
            "chromosome": text(gene.get("chromosome")),
            "biotype": text(gene.get("biotype")),
            "base_mean": number(gene["baseMean"], 1),
            "log2fc": number(gene["log2FoldChange"], 3),
            "lfc_se": number(gene["lfcSE"], 3),
            "pvalue": number(gene["pvalue"], 3, significant_digits=True),
            "padj": number(gene["padj"], 3, significant_digits=True),
            "immune": bool(gene["immune"]),
            "immune_lists": gene["immune_lists"].split(";") if gene["immune_lists"] else [],
            "atlas_markers": gene["atlas_markers"].split(";") if gene["atlas_markers"] else [],
            # model-estimated level per group: the arrow chart draws healthy -> ICP from these
            "group_estimates": {group: number(v, 4, significant_digits=True) for group, v in estimates.loc[gene_id].items()},
            # plain averages, overall and within each fetal sex; flagged when they point the other way than the model
            "plain_means": {group: number(v, 1) for group, v in plain_means.loc[gene_id].items()},
            "plain_log2fc": number(plain.at[gene_id, "all"], 3),
            "sex_log2fc": {sex: number(plain.at[gene_id, sex], 3) for sex in ("F", "M") if sex in plain},
            "raw_disagrees": bool(plain.at[gene_id, "all"] * gene["log2FoldChange"] < 0),
            "normalized": [number(v, 1) for v in normalized.loc[gene_id, sample_order]],
            "z": [number(v, 2) for v in zscores.loc[gene_id, sample_order]],
        })
    return records


def immune_list_summary(significant: pd.DataFrame, lists: pd.DataFrame) -> list[dict]:
    """Each ImmPort list holding at least one significant gene: how many are higher / lower."""
    higher, lower = Counter(), Counter()
    for memberships, log2fc in zip(significant["immune_lists"], significant["log2FoldChange"]):
        for list_id in filter(None, memberships.split(";")):
            (higher if log2fc > 0 else lower)[list_id] += 1
    rows = [{"id": list_id, "name": lists.at[list_id, "name"], "source": lists.at[list_id, "source"],
             "url": lists.at[list_id, "url"], "genes_in_list": int(lists.at[list_id, "genes_in_list"]),
             "genes_tested": int(lists.at[list_id, "genes_tested"]),
             "higher": higher[list_id], "lower": lower[list_id]}
            for list_id in set(higher) | set(lower)]
    return sorted(rows, key=lambda r: (-(r["higher"] + r["lower"]), r["name"]))


def cell_type_records(scores: pd.DataFrame, patient_scores: pd.DataFrame, annotated: pd.DataFrame,
                      sample_order: list[str]) -> list[dict]:
    """Every scored atlas population: its markers (with their own fold change), the difference
    between groups and each patient's score."""
    records = []
    for population_id, cell in scores.dropna(subset=["pvalue"]).iterrows():
        gene_ids = cell["marker_gene_ids"].split(";")
        found = cell["markers_found"].split(",")
        markers = [{"symbol": symbol, "gene_id": gene_id,
                    "log2fc": number(annotated.at[gene_id, "log2FoldChange"], 3),
                    "pvalue": number(annotated.at[gene_id, "pvalue"], 3, significant_digits=True)}
                   for symbol, gene_id in zip(found, gene_ids)]
        records.append({
            "id": population_id,
            "section": cell["section"],
            "name": cell["population"],
            "kind": cell["kind"],
            "markers": markers,
            "markers_missing": [m for m in cell["markers"].split(",") if m not in found],
            "effect": number(cell["effect"], 3),
            "ci": [number(cell["ci_low"], 3), number(cell["ci_high"], 3)],
            "pvalue": number(cell["pvalue"], 3, significant_digits=True),
            "fdr": number(cell["fdr"], 3, significant_digits=True),
            "composition_r": number(cell["composition_r"], 3),
            "effect_composition_adjusted": number(cell["effect_composition_adjusted"], 3),
            "pvalue_composition_adjusted": number(cell["pvalue_composition_adjusted"], 3, significant_digits=True),
            "scores": [number(v, 3) for v in patient_scores.loc[population_id, sample_order]],
        })
    return records


def cell_type_csv(scores: pd.DataFrame, composition: dict | None) -> pd.DataFrame:
    """Spreadsheet version of the cell-type scores (all atlas populations; unscored ones have no p-value)."""
    table = scores.reset_index()[["population_id", "section", "population", "kind", "markers", "markers_found",
                                  "n_markers", "n_found", "markers_lower", "markers_higher", "effect", "ci_low",
                                  "ci_high", "pvalue", "fdr"]]
    table = table.rename(columns={"effect": "score_difference_sd", "ci_low": "ci95_low", "ci_high": "ci95_high"})
    if composition:
        name = composition["name"]
        table[f"r_with_{name}"] = scores["composition_r"].to_numpy()
        table[f"score_difference_adjusted_for_{name}"] = scores["effect_composition_adjusted"].to_numpy()
        table[f"pvalue_adjusted_for_{name}"] = scores["pvalue_composition_adjusted"].to_numpy()
    return table.sort_values("pvalue", na_position="last")


def significant_genes_csv(significant: pd.DataFrame, normalized: pd.DataFrame, estimates: pd.DataFrame,
                          plain: pd.DataFrame, samples: pd.DataFrame, lists: pd.DataFrame,
                          atlas_names: dict[str, str], compare: str) -> pd.DataFrame:
    """Spreadsheet version of the significant genes, including each sample's normalized count."""
    log2fc = significant["log2FoldChange"]
    table = pd.DataFrame({"gene_id": significant.index, "symbol": significant["gene_symbol"],
                          "current_symbol": significant["current_symbol"]})
    if "description" in significant:
        table["description"] = [clean_description(v) for v in significant["description"]]
    if "chromosome" in significant:
        table["chromosome"] = [text(v) for v in significant["chromosome"]]
    table = table.assign(**{
        "direction": np.where(log2fc > 0, f"higher in {compare}", f"lower in {compare}"),
        "log2_fold_change": log2fc.round(3),
        "fold_change": (2.0 ** log2fc).round(3),
        "mean_normalized_count": significant["baseMean"].round(1),
        "pvalue": significant["pvalue"],
        "padj": significant["padj"],
        "immune_related": significant["immune"],
        "immune_lists": ["; ".join(lists.at[i, "name"] for i in m.split(";")) if m else ""
                         for m in significant["immune_lists"]],
        "immune_cell_markers": ["; ".join(atlas_names[i] for i in m.split(";")) if m else ""
                                for m in significant["atlas_markers"]],
    })
    for group, level in estimates.loc[significant.index].items():
        table[f"{group} level (model estimate)"] = level.round(2)
    for group, mean in group_means(normalized.loc[significant.index], samples).items():
        table[f"{group} plain mean (normalized)"] = mean.round(1)
    shown = plain.loc[significant.index]
    table["log2 ratio of plain means"] = shown["all"].round(3)
    for sex, label in (("F", "female"), ("M", "male")):
        if sex in shown:
            table[f"log2 ratio of plain means, {label} placentas"] = shown[sex].round(3)
    table["plain means disagree with model"] = (shown["all"] * log2fc < 0)
    for sample, group in samples["group"].items():
        table[f"{sample} ({group}, normalized)"] = normalized.loc[significant.index, sample].round(1)
    return table.reset_index(drop=True)


def write_dataset_index(site_data_dir: Path) -> None:
    """Rebuild datasets.json from every dataset JSON in site/data/."""
    entries = []
    for path in sorted(site_data_dir.glob("*.json")):
        if path.name == "datasets.json":
            continue
        data = json.loads(path.read_text())
        entries.append({"id": data["id"], "title": data["title"], "tissue": data["tissue"],
                        "file": path.name, "significant_genes": len(data["genes"])})
    (site_data_dir / "datasets.json").write_text(json.dumps(entries, indent=2) + "\n")


def run(config_path: str | Path, results_root: Path = RESULTS_ROOT, site_data_dir: Path = SITE_DATA_DIR) -> None:
    cfg = load_config(config_path)
    folder = results_dir(cfg, results_root)
    reference, compare = cfg["reference_group"], cfg["compare_group"]

    annotated = read_tsv(folder / "de_all_genes_annotated.tsv")
    annotated[["immune_lists", "atlas_markers"]] = annotated[["immune_lists", "atlas_markers"]].fillna("")
    lists = read_tsv(folder / "immune_lists.tsv")
    cell_scores = read_tsv(folder / "cell_type_scores.tsv")
    cell_scores[["markers", "markers_found", "marker_gene_ids"]] = cell_scores[["markers", "markers_found", "marker_gene_ids"]].fillna("")
    patient_scores = read_tsv(folder / "cell_type_patient_scores.tsv")
    composition_path = folder / "composition_score.tsv"
    composition = read_tsv(composition_path)["score"] if composition_path.exists() else None
    normalized = read_tsv(folder / "normalized_counts.tsv")
    estimates = read_tsv(folder / "group_estimates.tsv")
    vst = read_tsv(folder / "vst.tsv")
    pca = read_tsv(folder / "pca.tsv")
    samples = ordered_samples(read_tsv(folder / "samples.tsv"), reference, compare)
    info = json.loads((folder / "step2_info.json").read_text())

    significant = annotated[is_significant(annotated, cfg["cutoffs"])].sort_values("pvalue")
    sample_order = list(samples.index)
    zscores = row_zscores(vst.loc[significant.index])
    plain_means = group_means(normalized.loc[significant.index], samples)
    plain = plain_changes(normalized.loc[significant.index], samples, reference, compare)
    sex_counts = (samples.groupby(["fetal_sex", "group"]).size().unstack(fill_value=0)
                  if "fetal_sex" in samples else None)

    data = {
        "id": cfg["id"],
        "title": cfg["title"],
        "tissue": cfg["tissue"],
        "geo_accession": cfg["geo_accession"],
        "reference_group": reference,
        "compare_group": compare,
        "method": info["method"],
        "software": f"pydeseq2 {info['pydeseq2_version']}",
        "design": info["design"],
        "cutoffs": cfg["cutoffs"],
        "genes_in_file": info["genes_in_file"],
        "genes_tested": info["genes_tested"],
        "immune_genes_tested": int(annotated["immune"].sum()),
        "pca_variance": info["pca_variance"],
        "samples": [{"id": sample,
                     "group": samples.at[sample, "group"],
                     "fetal_sex": samples.at[sample, "fetal_sex"] if "fetal_sex" in samples else None,
                     "library_size": int(samples.at[sample, "library_size"]),
                     "pc1": number(pca.at[sample, "PC1"]),
                     "pc2": number(pca.at[sample, "PC2"])} for sample in sample_order],
        "genes": gene_records(significant, normalized, zscores, estimates, plain_means, plain, sample_order),
        # patients per fetal sex and group, e.g. {"F": {"Healthy": 4, "ICP": 2}, "M": {...}}
        "sex_counts": None if sex_counts is None else {sex: {g: int(n) for g, n in row.items()}
                                                       for sex, row in sex_counts.iterrows()},
        "immune_lists": immune_list_summary(significant, lists),
        "atlas_sections": list(dict.fromkeys(cell_scores["section"])),
        "atlas_populations": {pid: {"name": cell["population"], "section": cell["section"], "kind": cell["kind"]}
                              for pid, cell in cell_scores.iterrows()},
        "cell_types": cell_type_records(cell_scores, patient_scores, annotated, sample_order),
        "composition_check": None if composition is None else {
            "name": cfg["composition_check"]["name"], "markers": cfg["composition_check"]["markers"],
            "scores": [number(composition[s], 3) for s in sample_order]},
        "generated": date.today().isoformat(),
    }

    site_data_dir = Path(site_data_dir)
    site_data_dir.mkdir(parents=True, exist_ok=True)
    json_path = site_data_dir / f"{cfg['id']}.json"
    json_path.write_text(json.dumps(data, allow_nan=False, separators=(",", ":")) + "\n")
    csv_path = site_data_dir / f"{cfg['id']}_significant_genes.csv"
    atlas_names = cell_scores["population"].to_dict()
    significant_genes_csv(significant, normalized, estimates, plain, samples, lists, atlas_names, compare).to_csv(csv_path, index=False)
    cell_csv_path = site_data_dir / f"{cfg['id']}_cell_type_scores.csv"
    cell_type_csv(cell_scores, cfg["composition_check"] if composition is not None else None).to_csv(cell_csv_path, index=False)
    write_dataset_index(site_data_dir)

    higher = int((significant["log2FoldChange"] > 0).sum())
    print(f"Significant genes ({describe_cutoffs(cfg['cutoffs'])}): {len(significant)} "
          f"({higher} higher in {compare}, {len(significant) - higher} lower), "
          f"{int(significant['immune'].sum())} immune-related")
    flagged = [g["symbol"] for g in data["genes"] if g["raw_disagrees"]]
    print(f"Genes whose plain averages point the other way than the model (flagged on the site): "
          f"{len(flagged)} {flagged}")
    print(f"Immune cell types: {len(data['cell_types'])} scored atlas populations exported")
    print(f"Wrote {show_path(json_path)}, {csv_path.name}, {cell_csv_path.name}, datasets.json")


if __name__ == "__main__":
    args = step_arguments("Step 4: export significant genes for the website.",
                          {"--site-data-dir": (SITE_DATA_DIR, "where the website data goes (default: site/data/)")})
    run(args.config, args.results_root, args.site_data_dir)
