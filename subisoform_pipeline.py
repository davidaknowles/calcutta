from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp

from differential_usage import test_differential_subisoform_usage
from grouped_em import GroupedEmResult, run_grouped_em
from subisoforms import SubisoformModel, build_subisoform_model


@dataclass
class SubisoformPipelineResult:
    grouped_em: GroupedEmResult
    subisoform_model: SubisoformModel
    subisoform_tpm: np.ndarray
    differential_usage: pd.DataFrame

    def subisoform_tpm_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            self.subisoform_tpm,
            index=self.grouped_em.sample_ids,
            columns=self.subisoform_model.subisoform_ids,
        )


def run_subisoform_pipeline(
    cell_ec_matrix: sp.spmatrix,
    ec_transcript_mat: sp.spmatrix,
    transcript_ids: Sequence[object],
    transcript_weights: Sequence[float],
    group_labels: Sequence[object],
    gtf_filename: str | Path,
    batch_labels: Optional[Sequence[object]] = None,
    metacell_size: Optional[int] = None,
    random_state: int = 0,
    min_group_count: int = 1,
    min_ec_total: int = 1,
    em_iterations: int = 100,
    em_tol: float = 1e-8,
    dirichlet_pseudocount: float = 1e-6,
    min_subisoforms: int = 2,
    min_samples_per_condition: int = 2,
) -> SubisoformPipelineResult:
    grouped = run_grouped_em(
        cell_ec_matrix=cell_ec_matrix,
        ec_transcript_mat=ec_transcript_mat,
        transcript_weights=transcript_weights,
        group_labels=group_labels,
        transcript_ids=transcript_ids,
        batch_labels=batch_labels,
        metacell_size=metacell_size,
        random_state=random_state,
        min_group_count=min_group_count,
        min_ec_total=min_ec_total,
        iterations=em_iterations,
        tol=em_tol,
    )

    subisoform_model = build_subisoform_model(gtf_filename=gtf_filename, transcript_ids=grouped.transcript_ids)
    subisoform_tpm = subisoform_model.collapse(grouped.tpm, transcript_ids=grouped.transcript_ids)
    differential_usage = test_differential_subisoform_usage(
        subisoform_values=subisoform_tpm,
        subisoform_table=subisoform_model.subisoform_table,
        condition_labels=grouped.sample_metadata["condition"].to_numpy(),
        sample_ids=grouped.sample_ids,
        pseudocount=dirichlet_pseudocount,
        min_subisoforms=min_subisoforms,
        min_samples_per_condition=min_samples_per_condition,
    )

    return SubisoformPipelineResult(
        grouped_em=grouped,
        subisoform_model=subisoform_model,
        subisoform_tpm=subisoform_tpm,
        differential_usage=differential_usage,
    )
