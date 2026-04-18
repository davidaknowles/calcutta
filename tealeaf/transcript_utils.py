"""Transcript-length and alias-resolution helpers for transcript EM."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from . import calcutta


def load_alias_map(alias_path: str | Path) -> dict[str, str]:
    """Load transcript identifier aliases from a two-column TSV file."""
    alias_df = pd.read_csv(alias_path, sep="\t", names=["source", "target"])
    return dict(zip(alias_df["source"].astype(str), alias_df["target"].astype(str)))


def resolve_transcript_lengths(
    transcript_ids: Sequence[object],
    transcript_lengths: Mapping[str, int],
    alias_map: Mapping[str, str] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Resolve transcript lengths, optionally falling back through aliases."""
    transcript_ids = np.asarray(transcript_ids, dtype=object)
    lengths = np.zeros(len(transcript_ids), dtype=float)
    missing: list[object] = []
    for idx, transcript_id in enumerate(transcript_ids):
        transcript_key = str(transcript_id)
        resolved = transcript_key
        if resolved not in transcript_lengths and alias_map is not None:
            resolved = alias_map.get(resolved, resolved)
        if resolved not in transcript_lengths:
            missing.append(transcript_id)
            continue
        lengths[idx] = float(transcript_lengths[resolved])
    return lengths, np.asarray(missing, dtype=object)


def get_effective_length_weights(
    transcript_ids: Sequence[object],
    fasta_path: str | Path,
    alias_path: str | Path | None = None,
    fragment_size: int = 300,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert transcript lengths into effective-length EM weights."""
    transcript_lengths = calcutta.get_transcript_lengths(Path(fasta_path))
    alias_map = load_alias_map(alias_path) if alias_path is not None else None
    feature_lengths, missing = resolve_transcript_lengths(
        transcript_ids=transcript_ids,
        transcript_lengths=transcript_lengths,
        alias_map=alias_map,
    )
    if missing.size:
        raise KeyError(f"Transcript lengths not found for: {missing.tolist()}")
    effective_lengths = np.where(feature_lengths > fragment_size, feature_lengths - fragment_size, feature_lengths)
    return feature_lengths, 1.0 / effective_lengths
