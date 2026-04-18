"""Dirichlet-model fitting and differential subisoform usage testing."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import digamma, gammaln
from scipy.stats import chi2


def _as_2d_array(values: np.ndarray | Sequence[Sequence[float]]) -> np.ndarray:
    """Coerce an input matrix to a two-dimensional float array."""
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 2:
        raise ValueError("values must be two-dimensional.")
    return arr


def _smooth_compositions(values: np.ndarray, pseudocount: float) -> np.ndarray:
    """Convert non-negative counts into smoothed per-sample compositions."""
    if pseudocount <= 0:
        raise ValueError("pseudocount must be positive.")
    values = np.asarray(values, dtype=float)
    if np.any(values < 0):
        raise ValueError("values must be non-negative.")
    totals = values.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise ValueError("Each sample must have positive total abundance.")
    return (values + pseudocount) / (totals + pseudocount * values.shape[1])


def fit_dirichlet(compositions: np.ndarray, maxiter: int = 500) -> tuple[np.ndarray, float, bool]:
    """Fit Dirichlet concentration parameters to observed compositions."""
    compositions = _as_2d_array(compositions)
    log_x = np.log(compositions)
    log_x_sum = log_x.sum(axis=0)
    mean = compositions.mean(axis=0)
    variance = compositions.var(axis=0)
    precision = np.median((mean * (1.0 - mean) / np.maximum(variance, 1e-8)) - 1.0)
    precision = float(max(precision, 1.0))
    alpha0 = np.clip(mean * precision, 1e-3, None)
    theta0 = np.log(alpha0)

    n_samples = compositions.shape[0]

    def objective(theta: np.ndarray) -> tuple[float, np.ndarray]:
        alpha = np.exp(theta)
        log_likelihood = n_samples * (gammaln(alpha.sum()) - gammaln(alpha).sum()) + ((alpha - 1.0) * log_x_sum).sum()
        grad_alpha = n_samples * (digamma(alpha.sum()) - digamma(alpha)) + log_x_sum
        grad_theta = -(grad_alpha * alpha)
        return -float(log_likelihood), grad_theta

    result = minimize(
        lambda theta: objective(theta)[0],
        theta0,
        jac=lambda theta: objective(theta)[1],
        method="L-BFGS-B",
        options={"maxiter": maxiter},
    )
    alpha = np.exp(result.x)
    log_likelihood = -float(result.fun)
    return alpha, log_likelihood, bool(result.success)


def benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    """Adjust p-values with the Benjamini-Hochberg FDR procedure."""
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / (np.arange(len(ranked)) + 1.0)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result


@dataclass
class DirichletLrtResult:
    """Summary statistics for a Dirichlet likelihood-ratio test."""

    alpha_null: np.ndarray
    alpha_by_condition: dict[object, np.ndarray]
    log_likelihood_null: float
    log_likelihood_alt: float
    statistic: float
    degrees_of_freedom: int
    p_value: float
    success: bool


def dirichlet_lrt(
    values: np.ndarray,
    condition_labels: Sequence[object],
    pseudocount: float = 1e-6,
) -> DirichletLrtResult:
    """Compare shared and condition-specific Dirichlet composition models."""
    values = _as_2d_array(values)
    condition_labels = np.asarray(condition_labels, dtype=object)
    if values.shape[0] != len(condition_labels):
        raise ValueError("condition_labels has incompatible length.")

    compositions = _smooth_compositions(values, pseudocount=pseudocount)
    unique_conditions = pd.unique(condition_labels)
    if len(unique_conditions) < 2:
        raise ValueError("At least two conditions are required for a differential test.")

    alpha_null, ll_null, success_null = fit_dirichlet(compositions)
    alpha_by_condition: dict[object, np.ndarray] = {}
    ll_alt = 0.0
    success = success_null
    for condition in unique_conditions:
        subset = compositions[condition_labels == condition]
        alpha_condition, ll_condition, success_condition = fit_dirichlet(subset)
        alpha_by_condition[condition] = alpha_condition
        ll_alt += ll_condition
        success = success and success_condition

    statistic = max(0.0, 2.0 * (ll_alt - ll_null))
    dof = compositions.shape[1] * (len(unique_conditions) - 1)
    p_value = float(chi2.sf(statistic, dof))
    return DirichletLrtResult(
        alpha_null=alpha_null,
        alpha_by_condition=alpha_by_condition,
        log_likelihood_null=ll_null,
        log_likelihood_alt=ll_alt,
        statistic=statistic,
        degrees_of_freedom=dof,
        p_value=p_value,
        success=success,
    )


def test_differential_subisoform_usage(
    subisoform_values: np.ndarray,
    subisoform_table: pd.DataFrame,
    condition_labels: Sequence[object],
    sample_ids: Sequence[object] | None = None,
    pseudocount: float = 1e-6,
    min_subisoforms: int = 2,
    min_samples_per_condition: int = 2,
) -> pd.DataFrame:
    """Run event-wise Dirichlet LRTs across groups of subisoforms."""
    subisoform_values = _as_2d_array(subisoform_values)
    condition_labels = np.asarray(condition_labels, dtype=object)
    if subisoform_values.shape[0] != len(condition_labels):
        raise ValueError("condition_labels has incompatible length.")
    if sample_ids is None:
        sample_ids = np.arange(subisoform_values.shape[0], dtype=object)
    sample_ids = np.asarray(sample_ids, dtype=object)

    results: list[dict[str, object]] = []
    for event_id, event_df in subisoform_table.groupby("event_id", sort=False):
        if len(event_df) < min_subisoforms:
            continue
        event_indices = event_df.index.to_numpy()
        event_matrix = subisoform_values[:, event_indices]
        positive = event_matrix.sum(axis=1) > 0
        event_conditions = condition_labels[positive]
        if len(pd.unique(event_conditions)) < 2:
            continue
        counts_per_condition = pd.Series(event_conditions).value_counts()
        keep_conditions = counts_per_condition[counts_per_condition >= min_samples_per_condition].index.to_numpy()
        keep = positive & np.isin(condition_labels, keep_conditions)
        if len(pd.unique(condition_labels[keep])) < 2:
            continue

        event_matrix = event_matrix[keep]
        event_conditions = condition_labels[keep]
        lrt = dirichlet_lrt(event_matrix, event_conditions, pseudocount=pseudocount)
        mean_compositions = _smooth_compositions(event_matrix, pseudocount=pseudocount)
        mean_df = pd.DataFrame(mean_compositions, columns=event_df["subisoform_id"].tolist())
        mean_df["condition"] = event_conditions
        condition_means = mean_df.groupby("condition", sort=False).mean()

        results.append(
            {
                "event_id": event_id,
                "gene_id": event_df["gene_id"].iloc[0],
                "n_subisoforms": len(event_df),
                "n_samples": int(keep.sum()),
                "n_conditions": len(pd.unique(event_conditions)),
                "subisoform_ids": tuple(event_df["subisoform_id"].tolist()),
                "p_value": lrt.p_value,
                "lrt_statistic": lrt.statistic,
                "degrees_of_freedom": lrt.degrees_of_freedom,
                "log_likelihood_null": lrt.log_likelihood_null,
                "log_likelihood_alt": lrt.log_likelihood_alt,
                "success": lrt.success,
                "condition_mean_composition": condition_means.to_dict(orient="index"),
            }
        )

    if not results:
        return pd.DataFrame(
            columns=[
                "event_id",
                "gene_id",
                "n_subisoforms",
                "n_samples",
                "n_conditions",
                "subisoform_ids",
                "p_value",
                "lrt_statistic",
                "degrees_of_freedom",
                "log_likelihood_null",
                "log_likelihood_alt",
                "success",
                "condition_mean_composition",
                "fdr",
            ]
        )

    result_df = pd.DataFrame.from_records(results).sort_values("p_value").reset_index(drop=True)
    result_df["fdr"] = benjamini_hochberg(result_df["p_value"].to_numpy())
    return result_df
