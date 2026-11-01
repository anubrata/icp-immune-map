"""
Step 2 - Differential expression: which genes are expressed differently in ICP vs healthy?

RNA-seq counts (data_type: counts): DESeq2 (Love, Huber & Anders 2014), via pydeseq2.
  1. Drop genes with too few reads to test reliably.
  2. Fit a negative binomial model per gene, design "~ fetal_sex + condition" (see fetal sex below).
  3. Wald test of condition (compare group vs reference group).

Agilent microarrays (data_type: microarray): limma-style moderated t-test (analysis/moderated_t.py).
  1. log2 of each probe's signal, then quantile normalization: every array gets the same
     distribution of values, which removes technical differences in overall brightness.
  2. Keep probes that are clearly above background on at least as many arrays as the smallest
     group and that measure a named gene; for genes with several probes, keep the brightest one.
  3. Fit a linear model per gene and run the moderated t-test on the condition term.

Either way, each gene gets:
  log2FoldChange  > 0 higher in the compare group (ICP), < 0 lower. 1 = 2-fold, 2 = 4-fold.
  pvalue          unadjusted p-value
  padj            Benjamini-Hochberg adjusted p-value (false discovery rate)

Fetal sex: added to the design when the config has adjust_for_fetal_sex: true (the default) and
both sexes are present, so sex differences are not mistaken for ICP effects. When fetal sex is
recorded but not adjusted for, genes on chromosome Y are left out of the test, because they
would only reflect which placentas are male.

Inputs:  results/<dataset>/ tables from step 1
Outputs: results/<dataset>/de_all_genes.tsv       statistics for every tested gene, sorted by p-value
         results/<dataset>/group_estimates.tsv    each group's expression level estimated by the model
                                                  (other terms held at the cohort average)
         results/<dataset>/normalized_counts.tsv  (counts) counts adjusted for sequencing depth
           or normalized_intensities.tsv          (microarray) quantile-normalized signal
         results/<dataset>/vst.tsv                (counts) variance-stabilized expression, log2-like
           or log2_intensities.tsv                (microarray) the same on the log2 scale
         results/<dataset>/pca.tsv                each sample's position on the first two principal components
         results/<dataset>/step2_info.json        method, design, numbers of genes

Run:     uv run analysis/step2_differential_expression.py config/GSE282182.yaml
"""

from __future__ import annotations

import json
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from common import (RESULTS_ROOT, DataError, describe_cutoffs, expression_files, is_significant, load_config,
                    read_tsv, results_dir, show_path, step_arguments, write_tsv)
from moderated_t import fit_linear_models, moderated_t_test

PCA_TOP_GENES = 500  # like DESeq2's plotPCA: use the 500 most variable genes


# ------------------------------------------------------------------ shared by both paths

def design_formula(samples: pd.DataFrame, adjust_for_sex: bool = True) -> str:
    """'~ fetal_sex + condition' when adjusting and both sexes are present, otherwise '~ condition'."""
    if not adjust_for_sex or "fetal_sex" not in samples or samples["fetal_sex"].nunique() < 2:
        return "~ condition"
    # Sex and group can only be told apart if at least one group contains both sexes.
    if all(sexes.nunique() == 1 for _, sexes in samples.groupby("group")["fetal_sex"]):
        raise DataError("fetal sex and group are completely confounded (each group has a single sex), so their "
                        "effects cannot be separated. Set adjust_for_fetal_sex: false in the config to run without it.")
    return "~ fetal_sex + condition"


def chr_y_genes(genes: pd.DataFrame, samples: pd.DataFrame, design: str) -> pd.Index:
    """Genes on chromosome Y, to leave out when the samples' fetal sexes differ but the design
    does not adjust for sex (their differences would only reflect which placentas are male)."""
    sexes_differ = "fetal_sex" in samples and samples["fetal_sex"].nunique() > 1
    if not sexes_differ or "fetal_sex" in design or "chromosome" not in genes:
        return pd.Index([])
    return genes.index[genes["chromosome"].astype(str).str.upper().isin(["Y", "CHRY"])]


def sample_pca(log_scale: pd.DataFrame, n_top: int = PCA_TOP_GENES) -> tuple[pd.DataFrame, list[float]]:
    """Samples' positions on the first two principal components of the most variable genes."""
    top_genes = log_scale.var(axis=1).nlargest(n_top).index
    pca = PCA(n_components=2)
    coords = pca.fit_transform(log_scale.loc[top_genes].T)
    table = pd.DataFrame(coords, index=pd.Index(log_scale.columns, name="sample"), columns=["PC1", "PC2"])
    return table, pca.explained_variance_ratio_.tolist()


def model_terms(samples: pd.DataFrame, reference: str, design: str) -> pd.DataFrame:
    """Design matrix (samples x terms): intercept, fetal sex if in the design, one column per non-reference group."""
    terms = pd.DataFrame({"intercept": 1.0}, index=samples.index)
    if "fetal_sex" in design:
        terms["male fetus"] = (samples["fetal_sex"] == "M").astype(float)
    for group in samples["group"].unique():
        if group != reference:
            terms[f"group: {group}"] = (samples["group"] == group).astype(float)
    return terms


def group_levels(coefficients: np.ndarray, terms: pd.DataFrame, samples: pd.DataFrame, to_linear) -> pd.DataFrame:
    """Each group's level from fitted coefficients (genes x terms), with every term other than the
    group terms held at its average over all samples, so two groups' levels differ by exactly
    the fold change. `to_linear` turns the model's scale back into expression units."""
    group_columns = [c for c in terms.columns if c.startswith("group") or c.startswith("condition")]
    other_columns = [c for c in terms.columns if c not in group_columns]
    levels = {}
    for group, members in samples.groupby("group", sort=False).groups.items():
        x = terms.loc[members].mean()
        x[other_columns] = terms[other_columns].mean()
        levels[group] = to_linear(coefficients @ x.to_numpy())
    return pd.DataFrame(levels)


# ------------------------------------------------------------------ RNA-seq counts: DESeq2

def filter_low_counts(counts: pd.DataFrame, samples: pd.DataFrame, min_count: int) -> pd.DataFrame:
    """Keep genes with at least `min_count` reads in at least as many samples as the smallest group."""
    min_samples = int(samples["group"].value_counts().min())
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def run_deseq2(counts: pd.DataFrame, samples: pd.DataFrame, reference: str, compare: str, design: str):
    """Fit DESeq2 and test compare vs reference. Returns (results table, fitted model)."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    metadata = pd.DataFrame({"condition": samples["group"]}, index=samples.index)
    if "fetal_sex" in design:
        metadata["fetal_sex"] = samples["fetal_sex"]
    model = DeseqDataSet(counts=counts.T, metadata=metadata, design=design, quiet=True)
    model.deseq2()  # size factors, dispersions, fold changes, outlier (Cook's distance) handling
    stats = DeseqStats(model, contrast=["condition", compare, reference], quiet=True)
    stats.summary()  # Wald test p-values + Benjamini-Hochberg adjustment
    return stats.results_df.copy(), model


def analyze_counts(cfg: dict, folder: Path, samples: pd.DataFrame) -> dict:
    counts, genes = read_tsv(folder / "counts.tsv"), read_tsv(folder / "genes.tsv")
    reference, compare = cfg["reference_group"], cfg["compare_group"]
    design = design_formula(samples, cfg["adjust_for_fetal_sex"])
    tested = filter_low_counts(counts, samples, cfg["min_count"])
    min_samples = int(samples["group"].value_counts().min())
    print(f"Genes tested: {len(tested):,} of {len(counts):,} (the others have fewer than "
          f"{cfg['min_count']} reads in at least {min_samples} samples)")
    on_y = chr_y_genes(genes.loc[tested.index], samples, design)
    tested = tested.drop(on_y)

    results, model = run_deseq2(tested, samples, reference, compare, design)
    # The model's own design matrix and natural-log coefficients give each group's level.
    terms = model.obsm["design_matrix"]
    estimates = group_levels(model.varm["LFC"][terms.columns].to_numpy(), terms, samples, np.exp)
    estimates.index = pd.Index(model.var_names, name="gene_id")
    normalized = pd.DataFrame(model.layers["normed_counts"].T, index=tested.index, columns=tested.columns)
    model.vst(use_design=False)  # last: it refits dispersions on the same model object
    log_scale = pd.DataFrame(model.layers["vst_counts"].T, index=tested.index, columns=tested.columns)
    return {"table": genes.loc[results.index].join(results), "estimates": estimates, "normalized": normalized,
            "log_scale": log_scale, "design": design, "chr_y_excluded": len(on_y),
            "info": {"method": "DESeq2 Wald test", "software": f"pydeseq2 {version('pydeseq2')}",
                     "measure": "read count", "genes_in_file": len(counts), "min_count": cfg["min_count"],
                     "min_samples": min_samples}}


# ------------------------------------------------------------------ microarrays: moderated t-test

def quantile_normalize(log2: pd.DataFrame) -> pd.DataFrame:
    """Give every array the same distribution: the k-th smallest value on each array becomes the
    average of the k-th smallest values across arrays (limma's normalizeBetweenArrays, "quantile")."""
    mean_quantiles = np.sort(log2.to_numpy(), axis=0).mean(axis=1)
    ranks = log2.rank(axis=0, method="first").astype(int) - 1
    return pd.DataFrame({array: mean_quantiles[ranks[array].to_numpy()] for array in log2.columns}, index=log2.index)


def probes_to_genes(normalized: pd.DataFrame, above: pd.DataFrame, probes: pd.DataFrame, min_samples: int):
    """Keep probes clearly above background on >= min_samples arrays that measure a named gene, then
    one probe per gene: the brightest on average. Returns (genes x arrays expression, gene annotation)."""
    keep = (above.sum(axis=1) >= min_samples) & (probes["gene_symbol"] != "")
    candidates = normalized.loc[keep[keep].index]
    brightest = candidates.mean(axis=1).groupby(probes.loc[candidates.index, "gene_symbol"]).idxmax()
    expression = candidates.loc[brightest.to_numpy()]
    expression.index = pd.Index(brightest.index, name="gene_id")  # on arrays, genes are identified by symbol
    genes = pd.DataFrame({"gene_symbol": brightest.index, "probe_id": brightest.to_numpy(),
                          "description": probes.loc[brightest.to_numpy(), "gene_name"].to_numpy(),
                          "chromosome": probes.loc[brightest.to_numpy(), "chromosome"].to_numpy()},
                         index=expression.index)
    return expression, genes


def analyze_microarray(cfg: dict, folder: Path, samples: pd.DataFrame) -> dict:
    signal = read_tsv(folder / "intensities.tsv")[samples.index]
    above = read_tsv(folder / "above_background.tsv")[samples.index].astype(bool)
    probes = read_tsv(folder / "probes.tsv").fillna("")
    reference, compare = cfg["reference_group"], cfg["compare_group"]
    design = design_formula(samples, cfg["adjust_for_fetal_sex"])
    min_samples = int(samples["group"].value_counts().min())

    log2 = quantile_normalize(np.log2(signal.clip(lower=1.0)))  # signals below 1 are floored at 1 (log2 = 0)
    expression, genes = probes_to_genes(log2, above, probes, min_samples)
    genes_on_array = probes.loc[probes["gene_symbol"] != "", "gene_symbol"].nunique()
    print(f"Genes tested: {len(expression):,} of {genes_on_array:,} named genes on the array (the others are not "
          f"clearly above background on at least {min_samples} arrays)")
    on_y = chr_y_genes(genes, samples, design)
    expression, genes = expression.drop(on_y), genes.drop(on_y)

    terms = model_terms(samples, reference, design)
    results, prior = moderated_t_test(expression, terms, f"group: {compare}")
    coefficients = fit_linear_models(expression.to_numpy(), terms.loc[expression.columns].to_numpy())[0]
    estimates = group_levels(coefficients, terms, samples, np.exp2).set_index(expression.index)
    normalized = np.exp2(expression)
    results.insert(0, "baseMean", normalized.mean(axis=1))  # mean normalized intensity, like DESeq2's baseMean
    print(f"Moderated t-test: prior degrees of freedom {prior['prior_df']:.1f} added to the "
          f"{prior['residual_df']} residual degrees of freedom of each gene")
    return {"table": genes.join(results), "estimates": estimates, "normalized": normalized, "log_scale": expression,
            "design": design, "chr_y_excluded": len(on_y),
            "info": {"method": "limma-style moderated t-test on quantile-normalized log2 intensities",
                     "software": "analysis/moderated_t.py", "measure": "intensity",
                     "genes_in_file": int(genes_on_array), "min_samples": min_samples, **prior}}


# ------------------------------------------------------------------------------- run

def run(config_path: str | Path, results_root: Path = RESULTS_ROOT) -> None:
    cfg = load_config(config_path)
    folder = results_dir(cfg, results_root)
    samples = read_tsv(folder / "samples.tsv")
    reference, compare = cfg["reference_group"], cfg["compare_group"]

    result = (analyze_microarray if cfg["data_type"] == "microarray" else analyze_counts)(cfg, folder, samples)
    design = result["design"]
    print(f"Design: {design}" + ("" if "fetal_sex" in design or "fetal_sex" not in samples
                                 else "  (fetal sex recorded but not adjusted for)"))
    if "fetal_sex" in samples:
        print(pd.crosstab(samples["group"], samples["fetal_sex"]).to_string())
    if result["chr_y_excluded"]:
        print(f"Left out {result['chr_y_excluded']} chromosome Y genes (fetal sex differs but is not adjusted for)")

    table = result["table"].sort_values("pvalue")
    pca, pca_variance = sample_pca(result["log_scale"])
    normalized_file, log_file = expression_files(cfg)
    write_tsv(table, folder / "de_all_genes.tsv")
    write_tsv(result["estimates"], folder / "group_estimates.tsv")
    write_tsv(result["normalized"].round(2), folder / normalized_file)
    write_tsv(result["log_scale"].round(4), folder / log_file)
    write_tsv(pca.round(4), folder / "pca.tsv")
    info = {**result["info"], "design": design, "reference_group": reference, "compare_group": compare,
            "genes_tested": len(table), "chr_y_genes_excluded": result["chr_y_excluded"],
            "pca_variance": [round(v, 4) for v in pca_variance]}
    (folder / "step2_info.json").write_text(json.dumps(info, indent=2) + "\n")

    significant = is_significant(table, cfg["cutoffs"])
    higher = int((table.loc[significant, "log2FoldChange"] > 0).sum())
    print(f"\n{compare} vs {reference}:")
    print(f"  p < 0.05:            {(table['pvalue'] < 0.05).sum():>6,} genes")
    print(f"  FDR (padj) < 0.05:   {(table['padj'] < 0.05).sum():>6,} genes")
    print(f"  {describe_cutoffs(cfg['cutoffs'])}: {significant.sum():,} genes "
          f"({higher} higher in {compare}, {significant.sum() - higher} lower)  <- shown on the website")
    print(f"\nWrote {show_path(folder)}/de_all_genes.tsv, group_estimates.tsv, {normalized_file}, {log_file}, "
          "pca.tsv, step2_info.json")


if __name__ == "__main__":
    args = step_arguments("Step 2: differential expression (DESeq2 for counts, moderated t-test for microarrays).")
    run(args.config, args.results_root)
