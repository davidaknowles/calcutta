from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Optional, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp


def sparse_sum(matrix: sp.spmatrix, axis: int) -> np.ndarray:
    return np.squeeze(np.asarray(matrix.sum(axis)))


def _as_label_array(labels: Sequence[object], expected_length: int, name: str) -> np.ndarray:
    arr = np.asarray(labels, dtype=object)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional.")
    if len(arr) != expected_length:
        raise ValueError(f"{name} has length {len(arr)} but expected {expected_length}.")
    return arr


@dataclass
class EmRunResult:
    fractions: np.ndarray
    tpm: np.ndarray
    expected_counts: np.ndarray
    converged: bool
    iterations: int
    delta: float


@dataclass
class GroupedEmResult:
    sample_ids: np.ndarray
    transcript_ids: np.ndarray
    fractions: np.ndarray
    tpm: np.ndarray
    expected_counts: np.ndarray
    sample_metadata: pd.DataFrame
    transcript_mask: np.ndarray
    ec_mask: np.ndarray
    converged: np.ndarray
    iterations: np.ndarray
    deltas: np.ndarray

    def tpm_frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.tpm, index=self.sample_ids, columns=self.transcript_ids)

    def expected_count_frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.expected_counts, index=self.sample_ids, columns=self.transcript_ids)


def aggregate_rows(matrix: sp.spmatrix, row_labels: Sequence[object]) -> tuple[np.ndarray, sp.csr_matrix, np.ndarray]:
    matrix = sp.csr_matrix(matrix)
    row_labels = _as_label_array(row_labels, matrix.shape[0], "row_labels")
    codes, uniques = pd.factorize(row_labels, sort=False)
    aggregator = sp.csr_matrix(
        (np.ones(matrix.shape[0], dtype=float), (codes, np.arange(matrix.shape[0]))),
        shape=(len(uniques), matrix.shape[0]),
    )
    counts = np.bincount(codes, minlength=len(uniques))
    return np.asarray(uniques, dtype=object), (aggregator @ matrix).tocsr(), counts


def make_metacell_assignments(
    group_labels: Sequence[object],
    batch_labels: Optional[Sequence[object]] = None,
    target_size: int = 200,
    random_state: int = 0,
) -> pd.DataFrame:
    if target_size <= 0:
        raise ValueError("target_size must be positive.")

    group_labels = _as_label_array(group_labels, len(group_labels), "group_labels")
    if batch_labels is None:
        batch_labels = np.full(len(group_labels), "all", dtype=object)
        has_batch = False
    else:
        batch_labels = _as_label_array(batch_labels, len(group_labels), "batch_labels")
        has_batch = True

    rng = np.random.default_rng(random_state)
    metacell_labels = np.empty(len(group_labels), dtype=object)
    records: list[dict[str, object]] = []
    metacell_number = 0

    for group_value in pd.unique(group_labels):
        group_mask = group_labels == group_value
        for batch_value in pd.unique(batch_labels[group_mask]):
            indices = np.where(group_mask & (batch_labels == batch_value))[0]
            if len(indices) == 0:
                continue
            shuffled = indices.copy()
            rng.shuffle(shuffled)
            n_chunks = max(1, math.ceil(len(shuffled) / target_size))
            for chunk in np.array_split(shuffled, n_chunks):
                if len(chunk) == 0:
                    continue
                if has_batch:
                    metacell_id = f"{group_value}__{batch_value}__mc{metacell_number:03d}"
                else:
                    metacell_id = f"{group_value}__mc{metacell_number:03d}"
                metacell_labels[chunk] = metacell_id
                records.append(
                    {
                        "sample_id": metacell_id,
                        "condition": group_value,
                        "batch": batch_value if has_batch else None,
                        "n_cells": len(chunk),
                    }
                )
                metacell_number += 1

    assignment = pd.DataFrame(
        {
            "cell_index": np.arange(len(group_labels)),
            "condition": group_labels,
            "batch": batch_labels if has_batch else np.full(len(group_labels), None, dtype=object),
            "sample_id": metacell_labels,
        }
    )
    metadata = pd.DataFrame.from_records(records).set_index("sample_id")
    return assignment.merge(metadata.reset_index(), on=["sample_id", "condition", "batch"], how="left")


def build_grouped_ec_matrix(
    cell_ec_matrix: sp.spmatrix,
    group_labels: Sequence[object],
    batch_labels: Optional[Sequence[object]] = None,
    metacell_size: Optional[int] = None,
    random_state: int = 0,
) -> tuple[np.ndarray, sp.csr_matrix, pd.DataFrame]:
    cell_ec_matrix = sp.csr_matrix(cell_ec_matrix)
    group_labels = _as_label_array(group_labels, cell_ec_matrix.shape[0], "group_labels")

    if metacell_size is not None:
        assignments = make_metacell_assignments(
            group_labels=group_labels,
            batch_labels=batch_labels,
            target_size=metacell_size,
            random_state=random_state,
        )
        sample_ids, grouped_matrix, _ = aggregate_rows(cell_ec_matrix, assignments["sample_id"].to_numpy())
        sample_metadata = (
            assignments.groupby("sample_id", sort=False)
            .agg(condition=("condition", "first"), batch=("batch", "first"), n_cells=("cell_index", "size"))
            .reindex(sample_ids)
        )
    else:
        sample_ids, grouped_matrix, n_cells = aggregate_rows(cell_ec_matrix, group_labels)
        sample_metadata = pd.DataFrame({"condition": sample_ids, "batch": None, "n_cells": n_cells}, index=sample_ids)

    sample_metadata["total_count"] = sparse_sum(grouped_matrix, 1)
    return sample_ids, grouped_matrix.tocsr(), sample_metadata


def run_em(
    counts: Sequence[float],
    ec_transcript_mat: sp.spmatrix,
    transcript_weights: Sequence[float],
    iterations: int = 100,
    tol: float = 1e-8,
    init: Optional[Sequence[float]] = None,
) -> EmRunResult:
    counts = np.asarray(counts, dtype=float).ravel()
    weights = np.asarray(transcript_weights, dtype=float).ravel()
    compat = sp.coo_matrix(ec_transcript_mat)
    original_n_transcripts = compat.shape[1]
    active_columns = np.ones(original_n_transcripts, dtype=bool)

    if compat.shape[0] != counts.size:
        raise ValueError("counts and ec_transcript_mat have incompatible shapes.")
    if compat.shape[1] != weights.size:
        raise ValueError("transcript_weights and ec_transcript_mat have incompatible shapes.")

    positive_rows = counts > 0
    if not positive_rows.any():
        raise ValueError("At least one EC count must be positive.")
    if not np.all(positive_rows):
        compat = compat.tocsr()[positive_rows, :].tocoo()
        counts = counts[positive_rows]

    positive_cols = np.bincount(compat.col, minlength=compat.shape[1]) > 0
    if not np.all(positive_cols):
        active_columns &= positive_cols
        compat = compat.tocsc()[:, positive_cols].tocoo()
        weights = weights[positive_cols]
        if init is not None:
            init = np.asarray(init, dtype=float).ravel()[positive_cols]

    valid_weight_mask = np.isfinite(weights) & (weights > 0)
    if not np.all(valid_weight_mask):
        filtered_active = active_columns.copy()
        filtered_active[active_columns] = valid_weight_mask
        active_columns = filtered_active
        compat = compat.tocsc()[:, valid_weight_mask].tocoo()
        weights = weights[valid_weight_mask]
        if init is not None:
            init = np.asarray(init, dtype=float).ravel()[valid_weight_mask]

    if compat.shape[1] == 0:
        raise ValueError("No transcripts remain after filtering EM inputs.")

    if init is None:
        alpha = np.full(compat.shape[1], 1.0 / compat.shape[1], dtype=float)
    else:
        alpha = np.asarray(init, dtype=float).ravel()
        if alpha.size != compat.shape[1]:
            raise ValueError("init has incompatible length.")
        alpha = np.clip(alpha, np.finfo(float).tiny, None)
        alpha /= alpha.sum()

    final_expected = np.zeros_like(alpha)
    final_delta = np.inf
    converged = False

    for iteration in range(iterations):
        weighted = alpha * weights
        row_sums = np.bincount(compat.row, weights=weighted[compat.col], minlength=compat.shape[0])
        nonzero_entries = row_sums[compat.row] > 0
        expected = np.bincount(
            compat.col[nonzero_entries],
            weights=counts[compat.row[nonzero_entries]]
            * weighted[compat.col[nonzero_entries]]
            / row_sums[compat.row[nonzero_entries]],
            minlength=compat.shape[1],
        )
        total = expected.sum()
        if total <= 0:
            raise ValueError("EM update collapsed to zero total expected count.")
        alpha_new = expected / total
        final_delta = float(np.max(np.abs(alpha_new - alpha)))
        alpha = alpha_new
        final_expected = expected
        if final_delta <= tol:
            converged = True
            break

    fractions = np.zeros(original_n_transcripts, dtype=float)
    tpm = np.zeros(original_n_transcripts, dtype=float)
    expected_counts = np.zeros(original_n_transcripts, dtype=float)
    fractions[active_columns] = alpha
    tpm[active_columns] = alpha * 1e6
    expected_counts[active_columns] = final_expected

    return EmRunResult(
        fractions=fractions,
        tpm=tpm,
        expected_counts=expected_counts,
        converged=converged,
        iterations=iteration + 1,
        delta=final_delta,
    )


def run_grouped_em(
    cell_ec_matrix: sp.spmatrix,
    ec_transcript_mat: sp.spmatrix,
    transcript_weights: Sequence[float],
    group_labels: Sequence[object],
    transcript_ids: Optional[Sequence[object]] = None,
    batch_labels: Optional[Sequence[object]] = None,
    metacell_size: Optional[int] = None,
    random_state: int = 0,
    min_group_count: int = 1,
    min_ec_total: int = 1,
    iterations: int = 100,
    tol: float = 1e-8,
) -> GroupedEmResult:
    sample_ids, grouped_ec, sample_metadata = build_grouped_ec_matrix(
        cell_ec_matrix=cell_ec_matrix,
        group_labels=group_labels,
        batch_labels=batch_labels,
        metacell_size=metacell_size,
        random_state=random_state,
    )

    keep_samples = sample_metadata["total_count"].to_numpy() >= min_group_count
    sample_ids = sample_ids[keep_samples]
    grouped_ec = grouped_ec[keep_samples, :].tocsr()
    sample_metadata = sample_metadata.iloc[np.where(keep_samples)[0]].copy()
    sample_metadata.index = sample_ids

    ec_mask = sparse_sum(grouped_ec, 0) >= min_ec_total
    compat = sp.csr_matrix(ec_transcript_mat)[ec_mask, :]

    transcript_weights = np.asarray(transcript_weights, dtype=float).ravel()
    if transcript_ids is None:
        transcript_ids = np.arange(compat.shape[1])
    transcript_ids = np.asarray(transcript_ids, dtype=object)
    if transcript_ids.size != compat.shape[1]:
        raise ValueError("transcript_ids has incompatible length.")

    transcript_mask = sparse_sum(compat, 0) > 0
    transcript_mask &= np.isfinite(transcript_weights) & (transcript_weights > 0)
    compat = compat[:, transcript_mask].tocoo()
    transcript_weights = transcript_weights[transcript_mask]
    transcript_ids = transcript_ids[transcript_mask]

    fractions = np.zeros((grouped_ec.shape[0], compat.shape[1]), dtype=float)
    tpm = np.zeros_like(fractions)
    expected_counts = np.zeros_like(fractions)
    converged = np.zeros(grouped_ec.shape[0], dtype=bool)
    em_iterations = np.zeros(grouped_ec.shape[0], dtype=int)
    deltas = np.zeros(grouped_ec.shape[0], dtype=float)

    grouped_ec = grouped_ec[:, ec_mask].tocsr()
    for row_index in range(grouped_ec.shape[0]):
        counts = grouped_ec.getrow(row_index).toarray().ravel()
        result = run_em(
            counts=counts,
            ec_transcript_mat=compat,
            transcript_weights=transcript_weights,
            iterations=iterations,
            tol=tol,
        )
        fractions[row_index] = result.fractions
        tpm[row_index] = result.tpm
        expected_counts[row_index] = result.expected_counts
        converged[row_index] = result.converged
        em_iterations[row_index] = result.iterations
        deltas[row_index] = result.delta

    return GroupedEmResult(
        sample_ids=sample_ids,
        transcript_ids=transcript_ids,
        fractions=fractions,
        tpm=tpm,
        expected_counts=expected_counts,
        sample_metadata=sample_metadata,
        transcript_mask=transcript_mask,
        ec_mask=ec_mask,
        converged=converged,
        iterations=em_iterations,
        deltas=deltas,
    )
