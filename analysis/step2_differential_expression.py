"""
Step 2 - Differential expression: which genes are expressed differently in ICP vs healthy?

Method: DESeq2 (Love, Huber & Anders 2014, Genome Biology), via pydeseq2, its Python version.
  1. Drop genes with too few reads to test reliably.
  2. Fit a negative binomial model for each gene with the design "~ fetal_sex + condition".
     The fetal_sex term absorbs differences caused by the baby's sex, so they are not mistaken
     for ICP effects. It is left out when fetal sex is unknown or the same for every sample.
  3. Wald test of condition (compare group vs reference group). For each gene this gives:
       log2FoldChange  > 0 higher in the compare group (ICP), < 0 lower. 1 = 2-fold, 2 = 4-fold.
       pvalue          unadjusted p-value (left empty when one outlier sample drives the gene)
       padj            Benjamini-Hochberg adjusted p-value (false discovery rate)

Inputs:  results/<dataset>/counts.tsv, genes.tsv, samples.tsv   (from step 1)
Outputs: results/<dataset>/de_all_genes.tsv       statistics for every tested gene, sorted by p-value
         results/<dataset>/group_estimates.tsv    each group's expression level estimated by the model
                                                  (fetal sex held at the cohort average)
         results/<dataset>/normalized_counts.tsv  counts adjusted for sequencing depth (DESeq2 size factors)
         results/<dataset>/vst.tsv                variance-stabilized expression (log2-like), for heatmap/PCA
         results/<dataset>/pca.tsv                each sample's position on the first two principal components
         results/<dataset>/step2_info.json        design, software version, number of genes

Run:     uv run analysis/step2_differential_expression.py config/GSE282182.yaml
"""

from __future__ import annotations

import json
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA

from common import (RESULTS_ROOT, DataError, describe_cutoffs, is_significant, load_config, read_tsv,
                    results_dir, show_path, step_arguments, write_tsv)

PCA_TOP_GENES = 500  # like DESeq2's plotPCA: use the 500 most variable genes


def filter_low_counts(counts: pd.DataFrame, samples: pd.DataFrame, min_count: int) -> pd.DataFrame:
    """Keep genes with at least `min_count` reads in at least as many samples as the smallest group."""
    min_samples = int(samples["group"].value_counts().min())
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def design_formula(samples: pd.DataFrame) -> str:
    """'~ fetal_sex + condition' when both sexes are present, otherwise '~ condition'."""
    if "fetal_sex" not in samples or samples["fetal_sex"].nunique() < 2:
        return "~ condition"
    # Sex and group can only be told apart if at least one group contains both sexes.
    if all(sexes.nunique() == 1 for _, sexes in samples.groupby("group")["fetal_sex"]):
        raise DataError("fetal sex and group are completely confounded (each group has a single sex), so "
                        "their effects cannot be separated. Remove fetal_sex from the config to run without it.")
    return "~ fetal_sex + condition"


def run_deseq2(counts: pd.DataFrame, samples: pd.DataFrame, reference: str, compare: str):
    """Fit DESeq2 and test compare vs reference. Returns (results table, fitted model, design)."""
    design = design_formula(samples)
    metadata = pd.DataFrame({"condition": samples["group"]}, index=samples.index)
    if "fetal_sex" in design:
        metadata["fetal_sex"] = samples["fetal_sex"]

    model = DeseqDataSet(counts=counts.T, metadata=metadata, design=design, quiet=True)
    model.deseq2()  # size factors, dispersions, fold changes, outlier (Cook's distance) handling
    stats = DeseqStats(model, contrast=["condition", compare, reference], quiet=True)
    stats.summary()  # Wald test p-values + Benjamini-Hochberg adjustment
    return stats.results_df.copy(), model, design


def group_estimates(model: DeseqDataSet, samples: pd.DataFrame) -> pd.DataFrame:
    """Each group's expression level as estimated by the fitted model (normalized-count scale).

    Other model terms (fetal sex) are set to their average over all samples, so the ratio
    between two groups' estimates is exactly the model's fold change, without the sex imbalance.
    These are what the website's arrow chart draws.
    """
    design = model.obsm["design_matrix"]         # samples x model terms
    coefficients = model.varm["LFC"][design.columns].to_numpy()  # genes x terms, natural log
    group_terms = [term for term in design.columns if term.startswith("condition")]
    other_terms = [term for term in design.columns if term not in group_terms]
    estimates = {}
    for group, members in samples.groupby("group", sort=False).groups.items():
        terms = design.loc[members].mean()
        terms[other_terms] = design[other_terms].mean()
        estimates[group] = np.exp(coefficients @ terms.to_numpy())
    return pd.DataFrame(estimates, index=pd.Index(model.var_names, name="gene_id"))


def sample_pca(vst: pd.DataFrame, n_top: int = PCA_TOP_GENES) -> tuple[pd.DataFrame, list[float]]:
    """Samples' positions on the first two principal components of the most variable genes."""
    top_genes = vst.var(axis=1).nlargest(n_top).index
    pca = PCA(n_components=2)
    coords = pca.fit_transform(vst.loc[top_genes].T)
    table = pd.DataFrame(coords, index=pd.Index(vst.columns, name="sample"), columns=["PC1", "PC2"])
    return table, pca.explained_variance_ratio_.tolist()


def run(config_path: str | Path, results_root: Path = RESULTS_ROOT) -> None:
    cfg = load_config(config_path)
    folder = results_dir(cfg, results_root)
    counts = read_tsv(folder / "counts.tsv")
    genes = read_tsv(folder / "genes.tsv")
    samples = read_tsv(folder / "samples.tsv")
    reference, compare = cfg["reference_group"], cfg["compare_group"]

    tested = filter_low_counts(counts, samples, cfg["min_count"])
    min_samples = int(samples["group"].value_counts().min())
    print(f"Genes tested: {len(tested):,} of {len(counts):,} (the others have fewer than "
          f"{cfg['min_count']} reads in at least {min_samples} samples)")

    results, model, design = run_deseq2(tested, samples, reference, compare)
    print(f"Design: {design}")
    if "fetal_sex" in design:
        print(pd.crosstab(samples["group"], samples["fetal_sex"]).to_string())

    # Group estimates, normalized counts and the variance-stabilizing transform (VST) for plots.
    # The VST is computed last because it refits dispersions on the same model object.
    estimates = group_estimates(model, samples)
    normalized = pd.DataFrame(model.layers["normed_counts"].T, index=tested.index, columns=tested.columns)
    model.vst(use_design=False)
    vst = pd.DataFrame(model.layers["vst_counts"].T, index=tested.index, columns=tested.columns)
    pca, pca_variance = sample_pca(vst)

    table = genes.loc[results.index].join(results).sort_values("pvalue")
    write_tsv(table, folder / "de_all_genes.tsv")
    write_tsv(estimates, folder / "group_estimates.tsv")
    write_tsv(normalized.round(2), folder / "normalized_counts.tsv")
    write_tsv(vst.round(4), folder / "vst.tsv")
    write_tsv(pca.round(4), folder / "pca.tsv")
    info = {
        "method": "DESeq2 (pydeseq2), Wald test",
        "pydeseq2_version": version("pydeseq2"),
        "design": design,
        "reference_group": reference,
        "compare_group": compare,
        "genes_in_file": len(counts),
        "genes_tested": len(tested),
        "min_count": cfg["min_count"],
        "min_samples": min_samples,
        "pca_variance": [round(v, 4) for v in pca_variance],
    }
    (folder / "step2_info.json").write_text(json.dumps(info, indent=2) + "\n")

    significant = is_significant(table, cfg["cutoffs"])
    higher = int((table.loc[significant, "log2FoldChange"] > 0).sum())
    print(f"\n{compare} vs {reference}:")
    print(f"  p < 0.05:            {(table['pvalue'] < 0.05).sum():>6,} genes")
    print(f"  FDR (padj) < 0.05:   {(table['padj'] < 0.05).sum():>6,} genes")
    print(f"  {describe_cutoffs(cfg['cutoffs'])}: {significant.sum():,} genes "
          f"({higher} higher in {compare}, {significant.sum() - higher} lower)  <- shown on the website")
    print(f"\nWrote {show_path(folder)}/de_all_genes.tsv, group_estimates.tsv, normalized_counts.tsv, vst.tsv, "
          "pca.tsv, step2_info.json")


if __name__ == "__main__":
    args = step_arguments("Step 2: differential expression with DESeq2 (pydeseq2).")
    run(args.config, args.results_root)
