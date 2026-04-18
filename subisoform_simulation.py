from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import scipy.sparse as sp

from differential_usage import dirichlet_lrt
from subisoform_pipeline import SubisoformPipelineResult, run_subisoform_pipeline
from subisoforms import SubisoformModel, build_subisoform_model


@dataclass
class SimulationResult:
    values: np.ndarray
    condition_labels: np.ndarray
    sample_ids: np.ndarray


def simulate_dirichlet_multinomial(
    alpha_by_condition: dict[object, Sequence[float]],
    n_samples_per_condition: int | dict[object, int],
    total_count: int | Sequence[int] = 500,
    random_state: int = 0,
) -> SimulationResult:
    if isinstance(n_samples_per_condition, int):
        n_samples_lookup = {condition: n_samples_per_condition for condition in alpha_by_condition}
    else:
        n_samples_lookup = {condition: int(value) for condition, value in n_samples_per_condition.items()}

    rng = np.random.default_rng(random_state)
    values: list[np.ndarray] = []
    condition_labels: list[object] = []
    sample_ids: list[str] = []

    for condition, alpha in alpha_by_condition.items():
        alpha_arr = np.asarray(alpha, dtype=float)
        if alpha_arr.ndim != 1 or np.any(alpha_arr <= 0):
            raise ValueError("Each alpha vector must be one-dimensional and strictly positive.")
        n_samples = n_samples_lookup[condition]
        if isinstance(total_count, int):
            totals = np.full(n_samples, total_count, dtype=int)
        else:
            totals = np.asarray(total_count, dtype=int)
            if totals.size != n_samples:
                raise ValueError("total_count sequence must match the number of samples per condition.")
        compositions = rng.dirichlet(alpha_arr, size=n_samples)
        counts = np.vstack([rng.multinomial(int(total), comp) for total, comp in zip(totals, compositions)])
        values.append(counts)
        condition_labels.extend([condition] * n_samples)
        sample_ids.extend([f"{condition}_{idx:03d}" for idx in range(n_samples)])

    return SimulationResult(
        values=np.vstack(values),
        condition_labels=np.asarray(condition_labels, dtype=object),
        sample_ids=np.asarray(sample_ids, dtype=object),
    )


@dataclass
class LrtOperatingCharacteristics:
    n_reps: int
    alpha_threshold: float
    null_p_values: np.ndarray
    alt_p_values: np.ndarray

    @property
    def type1_error(self) -> float:
        return float(np.mean(self.null_p_values < self.alpha_threshold))

    @property
    def power(self) -> float:
        return float(np.mean(self.alt_p_values < self.alpha_threshold))


@dataclass
class GtfEventSimulationDesign:
    gtf_filename: Path
    subisoform_model: SubisoformModel
    target_event_id: str
    event_subisoform_ids: np.ndarray
    event_subisoform_to_transcript: np.ndarray
    ec_transcript_mat: sp.csr_matrix
    transcript_ec_probs: np.ndarray

    @property
    def transcript_ids(self) -> np.ndarray:
        return self.subisoform_model.transcript_ids


@dataclass
class EndToEndSimulationResult:
    design: GtfEventSimulationDesign
    cell_ec_matrix: sp.csr_matrix
    group_labels: np.ndarray
    transcript_counts: np.ndarray
    transcript_fractions: np.ndarray
    subisoform_fractions: np.ndarray


def _build_feature_compatibility(
    transcript_segment_orders: dict[str, tuple[int, ...]],
    transcript_ids: np.ndarray,
) -> tuple[sp.csr_matrix, np.ndarray]:
    ordered_ids = [str(tid) for tid in transcript_ids]
    transcript_features: dict[str, list[tuple[object, ...]]] = {}
    feature_members: dict[tuple[object, ...], set[str]] = {}

    for transcript_id in ordered_ids:
        path = transcript_segment_orders[transcript_id]
        features = [("segment", int(segment_order)) for segment_order in path]
        features.extend(("junction", int(left), int(right)) for left, right in zip(path[:-1], path[1:]))
        transcript_features[transcript_id] = features
        for feature in features:
            feature_members.setdefault(feature, set()).add(transcript_id)

    ec_lookup: dict[tuple[int, ...], int] = {}
    transcript_ec_probs = np.zeros((len(transcript_ids), 0), dtype=float)
    feature_ec_index: dict[tuple[object, ...], int] = {}

    for feature, members in feature_members.items():
        ec_key = tuple(idx for idx, transcript_id in enumerate(ordered_ids) if transcript_id in members)
        if ec_key not in ec_lookup:
            ec_lookup[ec_key] = len(ec_lookup)
            transcript_ec_probs = np.pad(transcript_ec_probs, ((0, 0), (0, 1)))
        feature_ec_index[feature] = ec_lookup[ec_key]

    for transcript_index, transcript_id in enumerate(ordered_ids):
        features = transcript_features[transcript_id]
        if not features:
            raise ValueError(f"Transcript {transcript_id} has no observable features.")
        weight = 1.0 / len(features)
        for feature in features:
            transcript_ec_probs[transcript_index, feature_ec_index[feature]] += weight

    rows: list[int] = []
    cols: list[int] = []
    for ec_key, ec_index in ec_lookup.items():
        rows.extend([ec_index] * len(ec_key))
        cols.extend(ec_key)
    ec_transcript_mat = sp.csr_matrix(
        (np.ones(len(rows), dtype=float), (rows, cols)),
        shape=(len(ec_lookup), len(transcript_ids)),
    )
    return ec_transcript_mat, transcript_ec_probs


def _select_event(
    model: SubisoformModel,
    event_id: Optional[str] = None,
) -> str:
    event_groups = list(model.subisoform_table.groupby("event_id", sort=False))
    if not event_groups:
        raise ValueError("The selected transcripts did not produce any subisoform events.")

    if event_id is not None:
        if event_id not in model.subisoform_table["event_id"].values:
            raise ValueError(f"Requested event_id {event_id} was not found in the subisoform model.")
        return str(event_id)

    n_transcripts = len(model.transcript_ids)
    candidates: list[tuple[int, int, str]] = []
    for event_name, event_table in event_groups:
        n_subisoforms = len(event_table)
        covered = sorted({tid for members in event_table["transcript_ids"] for tid in members})
        if n_subisoforms < 2:
            continue
        candidates.append((int(len(covered) == n_transcripts), n_subisoforms, str(event_name)))
    if not candidates:
        raise ValueError("No multi-subisoform events were found for the selected transcripts.")
    candidates.sort(key=lambda item: (-item[0], -item[1], item[2]))
    return candidates[0][2]


def build_gtf_event_simulation_design(
    gtf_filename: str | Path,
    transcript_ids: Sequence[object],
    event_id: Optional[str] = None,
) -> GtfEventSimulationDesign:
    gtf_path = Path(gtf_filename)
    model = build_subisoform_model(gtf_path, transcript_ids=transcript_ids)
    selected_event_id = _select_event(model, event_id=event_id)
    event_table = model.subisoform_table.loc[model.subisoform_table["event_id"] == selected_event_id].copy()
    if len(event_table) < 2:
        raise ValueError("Simulation requires an event with at least two subisoforms.")

    transcript_index = {str(tid): idx for idx, tid in enumerate(model.transcript_ids)}
    event_subisoform_to_transcript = np.zeros((len(event_table), len(model.transcript_ids)), dtype=float)
    for row_index, row in enumerate(event_table.itertuples(index=False)):
        members = [str(tid) for tid in row.transcript_ids]
        weight = 1.0 / len(members)
        for transcript_id in members:
            event_subisoform_to_transcript[row_index, transcript_index[transcript_id]] = weight

    ec_transcript_mat, transcript_ec_probs = _build_feature_compatibility(
        transcript_segment_orders=model.transcript_segment_orders,
        transcript_ids=model.transcript_ids,
    )
    return GtfEventSimulationDesign(
        gtf_filename=gtf_path,
        subisoform_model=model,
        target_event_id=selected_event_id,
        event_subisoform_ids=event_table["subisoform_id"].to_numpy(),
        event_subisoform_to_transcript=event_subisoform_to_transcript,
        ec_transcript_mat=ec_transcript_mat,
        transcript_ec_probs=transcript_ec_probs,
    )


def simulate_gtf_event_counts(
    design: GtfEventSimulationDesign,
    alpha_by_condition: dict[object, Sequence[float]],
    n_cells_per_condition: int | dict[object, int],
    mean_gene_count: float = 150.0,
    random_state: int = 0,
) -> EndToEndSimulationResult:
    if isinstance(n_cells_per_condition, int):
        n_cells_lookup = {condition: n_cells_per_condition for condition in alpha_by_condition}
    else:
        n_cells_lookup = {condition: int(value) for condition, value in n_cells_per_condition.items()}

    rng = np.random.default_rng(random_state)
    group_labels: list[object] = []
    transcript_counts: list[np.ndarray] = []
    transcript_fractions: list[np.ndarray] = []
    subisoform_fractions: list[np.ndarray] = []
    cell_ec_rows: list[np.ndarray] = []

    for condition, alpha in alpha_by_condition.items():
        alpha_arr = np.asarray(alpha, dtype=float)
        if alpha_arr.ndim != 1 or alpha_arr.size != len(design.event_subisoform_ids):
            raise ValueError("Each alpha vector must match the number of simulated event subisoforms.")
        if np.any(alpha_arr <= 0):
            raise ValueError("Each alpha vector must be strictly positive.")

        for _ in range(n_cells_lookup[condition]):
            subisoform_fraction = rng.dirichlet(alpha_arr)
            isoform_fraction = subisoform_fraction @ design.event_subisoform_to_transcript
            isoform_fraction /= isoform_fraction.sum()
            total_count = max(1, int(rng.poisson(mean_gene_count)))
            isoform_counts = rng.multinomial(total_count, isoform_fraction)
            ec_counts = np.zeros(design.ec_transcript_mat.shape[0], dtype=float)
            for transcript_index, transcript_count in enumerate(isoform_counts):
                if transcript_count <= 0:
                    continue
                ec_counts += rng.multinomial(int(transcript_count), design.transcript_ec_probs[transcript_index])
            group_labels.append(condition)
            transcript_counts.append(isoform_counts.astype(float))
            transcript_fractions.append(isoform_fraction)
            subisoform_fractions.append(subisoform_fraction)
            cell_ec_rows.append(ec_counts)

    return EndToEndSimulationResult(
        design=design,
        cell_ec_matrix=sp.csr_matrix(np.vstack(cell_ec_rows)),
        group_labels=np.asarray(group_labels, dtype=object),
        transcript_counts=np.vstack(transcript_counts),
        transcript_fractions=np.vstack(transcript_fractions),
        subisoform_fractions=np.vstack(subisoform_fractions),
    )


def run_end_to_end_simulation(
    design: GtfEventSimulationDesign,
    alpha_by_condition: dict[object, Sequence[float]],
    n_cells_per_condition: int | dict[object, int],
    mean_gene_count: float = 150.0,
    metacell_size: Optional[int] = 3,
    random_state: int = 0,
    dirichlet_pseudocount: float = 1e-6,
) -> tuple[EndToEndSimulationResult, SubisoformPipelineResult]:
    simulation = simulate_gtf_event_counts(
        design=design,
        alpha_by_condition=alpha_by_condition,
        n_cells_per_condition=n_cells_per_condition,
        mean_gene_count=mean_gene_count,
        random_state=random_state,
    )
    pipeline = run_subisoform_pipeline(
        cell_ec_matrix=simulation.cell_ec_matrix,
        ec_transcript_mat=design.ec_transcript_mat,
        transcript_ids=design.transcript_ids,
        transcript_weights=np.ones(len(design.transcript_ids), dtype=float),
        group_labels=simulation.group_labels,
        gtf_filename=design.gtf_filename,
        metacell_size=metacell_size,
        random_state=random_state,
        dirichlet_pseudocount=dirichlet_pseudocount,
        min_samples_per_condition=2,
    )
    return simulation, pipeline


def _event_p_value(pipeline: SubisoformPipelineResult, event_id: str) -> float:
    if pipeline.differential_usage.empty:
        return 1.0
    matches = pipeline.differential_usage.loc[pipeline.differential_usage["event_id"] == event_id, "p_value"]
    if matches.empty:
        return 1.0
    return float(matches.iloc[0])


def estimate_lrt_operating_characteristics(
    alpha_null: Sequence[float],
    alpha_alt_a: Sequence[float],
    alpha_alt_b: Sequence[float],
    n_samples_per_condition: int = 6,
    total_count: int = 500,
    n_reps: int = 200,
    alpha_threshold: float = 0.05,
    pseudocount: float = 1e-6,
    random_state: int = 0,
) -> LrtOperatingCharacteristics:
    rng = np.random.default_rng(random_state)
    null_p_values = np.zeros(n_reps, dtype=float)
    alt_p_values = np.zeros(n_reps, dtype=float)

    for idx in range(n_reps):
        null_seed = int(rng.integers(2**31 - 1))
        alt_seed = int(rng.integers(2**31 - 1))

        null_data = simulate_dirichlet_multinomial(
            alpha_by_condition={"a": alpha_null, "b": alpha_null},
            n_samples_per_condition=n_samples_per_condition,
            total_count=total_count,
            random_state=null_seed,
        )
        alt_data = simulate_dirichlet_multinomial(
            alpha_by_condition={"a": alpha_alt_a, "b": alpha_alt_b},
            n_samples_per_condition=n_samples_per_condition,
            total_count=total_count,
            random_state=alt_seed,
        )

        null_p_values[idx] = dirichlet_lrt(
            values=null_data.values,
            condition_labels=null_data.condition_labels,
            pseudocount=pseudocount,
        ).p_value
        alt_p_values[idx] = dirichlet_lrt(
            values=alt_data.values,
            condition_labels=alt_data.condition_labels,
            pseudocount=pseudocount,
        ).p_value

    return LrtOperatingCharacteristics(
        n_reps=n_reps,
        alpha_threshold=alpha_threshold,
        null_p_values=null_p_values,
        alt_p_values=alt_p_values,
    )


def estimate_end_to_end_operating_characteristics(
    design: GtfEventSimulationDesign,
    alpha_null: Sequence[float],
    alpha_alt_a: Sequence[float],
    alpha_alt_b: Sequence[float],
    n_cells_per_condition: int = 12,
    mean_gene_count: float = 150.0,
    metacell_size: Optional[int] = 3,
    n_reps: int = 50,
    alpha_threshold: float = 0.05,
    random_state: int = 0,
) -> LrtOperatingCharacteristics:
    rng = np.random.default_rng(random_state)
    null_p_values = np.zeros(n_reps, dtype=float)
    alt_p_values = np.zeros(n_reps, dtype=float)

    for idx in range(n_reps):
        null_seed = int(rng.integers(2**31 - 1))
        alt_seed = int(rng.integers(2**31 - 1))

        _, null_pipeline = run_end_to_end_simulation(
            design=design,
            alpha_by_condition={"a": alpha_null, "b": alpha_null},
            n_cells_per_condition=n_cells_per_condition,
            mean_gene_count=mean_gene_count,
            metacell_size=metacell_size,
            random_state=null_seed,
        )
        _, alt_pipeline = run_end_to_end_simulation(
            design=design,
            alpha_by_condition={"a": alpha_alt_a, "b": alpha_alt_b},
            n_cells_per_condition=n_cells_per_condition,
            mean_gene_count=mean_gene_count,
            metacell_size=metacell_size,
            random_state=alt_seed,
        )

        null_p_values[idx] = _event_p_value(null_pipeline, design.target_event_id)
        alt_p_values[idx] = _event_p_value(alt_pipeline, design.target_event_id)

    return LrtOperatingCharacteristics(
        n_reps=n_reps,
        alpha_threshold=alpha_threshold,
        null_p_values=null_p_values,
        alt_p_values=alt_p_values,
    )
