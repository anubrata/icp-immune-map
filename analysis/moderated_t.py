"""
Moderated t-test for microarrays: the empirical Bayes method of limma (Smyth 2004), in Python.

With only a few arrays per group, each gene's own variance estimate is very noisy, so ordinary
t-tests swing wildly. The moderated t-test borrows strength across all genes: it fits one
distribution to all gene variances (a scaled inverse chi-square "prior"), then shrinks each
gene's variance toward it. Genes keep their own fold change; only the noise estimate is
stabilized, which also adds the prior's degrees of freedom to every test.

Smyth GK (2004). Linear models and empirical Bayes methods for assessing differential
expression in microarray experiments. Statistical Applications in Genetics and Molecular
Biology 3(1), Article 3. The functions mirror limma's lmFit, fitFDist, squeezeVar and eBayes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import special, stats


def fit_linear_models(y: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, np.ndarray, int, np.ndarray]:
    """Least squares for every gene at once (y: genes x samples, X: samples x terms).

    Returns (coefficients [genes x terms], residual variances, residual degrees of freedom,
    unscaled coefficient covariance (X'X)^-1).
    """
    unscaled = np.linalg.inv(X.T @ X)
    beta = y @ X @ unscaled
    residual = y - beta @ X.T
    dof = X.shape[0] - X.shape[1]
    return beta, (residual ** 2).sum(axis=1) / dof, dof, unscaled


def trigamma_inverse(x: float) -> float:
    """Solve trigamma(y) = x for y, by Newton's method (as limma's trigammaInverse)."""
    if x > 1e7:
        return 1 / np.sqrt(x)
    if x < 1e-6:
        return 1 / x
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        step = tri * (1 - tri / x) / special.polygamma(2, y)
        y += step
        if -step / y < 1e-8:
            break
    return y


def fit_prior(variances: np.ndarray, dof: int) -> tuple[float, float]:
    """Fit a scaled inverse chi-square distribution to the gene variances (limma's fitFDist).

    Returns (prior degrees of freedom d0, prior variance s0^2). d0 is infinite when the genes'
    variances are no more spread out than chance alone would make them.
    """
    s2 = np.maximum(variances, 0)
    floor = np.median(s2) or 1.0
    s2 = np.maximum(s2, 1e-5 * floor)   # avoid log(0) for genes with no variation
    e = np.log(s2) - special.digamma(dof / 2) + np.log(dof / 2)
    e_mean = e.mean()
    e_var = e.var(ddof=1) - special.polygamma(1, dof / 2)
    if e_var > 0:
        d0 = 2 * trigamma_inverse(e_var)
        return d0, float(np.exp(e_mean + special.digamma(d0 / 2) - np.log(d0 / 2)))
    return np.inf, float(np.exp(e_mean))


def moderated_t_test(expression: pd.DataFrame, design: pd.DataFrame, term: str) -> tuple[pd.DataFrame, dict]:
    """Test one design term for every gene.

    expression: genes x samples, log2 scale. design: samples x terms (same sample order).
    Returns a table with log2FoldChange (the term's coefficient), lfcSE, stat (moderated t),
    pvalue and padj (Benjamini-Hochberg), plus a dict with the prior that was used.
    """
    X = design.loc[expression.columns].to_numpy(dtype=float)
    beta, s2, dof, unscaled = fit_linear_models(expression.to_numpy(dtype=float), X)
    d0, s0_2 = fit_prior(s2, dof)
    shrunk = np.full_like(s2, s0_2) if np.isinf(d0) else (d0 * s0_2 + dof * s2) / (d0 + dof)
    total_dof = min(dof + d0, dof * len(s2))  # as limma: never more than all genes' residual df pooled

    j = list(design.columns).index(term)
    se = np.sqrt(shrunk * unscaled[j, j])
    t = beta[:, j] / se
    p = 2 * stats.t.sf(np.abs(t), total_dof)
    table = pd.DataFrame({"log2FoldChange": beta[:, j], "lfcSE": se, "stat": t, "pvalue": p,
                          "padj": stats.false_discovery_control(p)}, index=expression.index)
    return table, {"prior_df": float(d0), "prior_variance": s0_2, "residual_df": dof, "test_df": float(total_dof)}
