"""Build GTF-derived subisoform models from transcript exon structures."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import gzip
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp


def _open_text(path: Path):
    """Open plain-text or gzipped files for text iteration."""
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "r")


def _parse_gtf_attributes(attribute_field: str) -> dict[str, str]:
    """Parse a GTF attributes column into a key-value dictionary."""
    attributes: dict[str, str] = {}
    for item in attribute_field.strip().split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        attributes[key] = value.strip().strip('"')
    return attributes


def load_transcript_exons(
    gtf_filename: str | Path,
    feature_types: Iterable[str] = ("exon", "five_prime_utr", "three_prime_utr"),
) -> pd.DataFrame:
    """Load exon-like transcript features from a GTF into a tidy table."""
    feature_types = set(feature_types)
    records: list[dict[str, object]] = []
    with _open_text(Path(gtf_filename)) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            chrom, _, feature, start, end, _, strand, _, attributes = line.rstrip().split("\t")
            if feature not in feature_types:
                continue
            metadata = _parse_gtf_attributes(attributes)
            transcript_id = metadata.get("transcript_id")
            gene_id = metadata.get("gene_id")
            if transcript_id is None or gene_id is None:
                continue
            start_i = int(start)
            end_i = int(end)
            if end_i < start_i:
                continue
            records.append(
                {
                    "gene_id": gene_id,
                    "transcript_id": transcript_id,
                    "chrom": chrom,
                    "strand": strand,
                    "start": start_i,
                    "end": end_i,
                }
            )
    return pd.DataFrame.from_records(records)


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping or directly adjacent genomic intervals."""
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def _build_atomic_segments(gene_df: pd.DataFrame) -> pd.DataFrame:
    """Split a gene's exonic span into disjoint atomic segments."""
    breakpoints = sorted(set(gene_df["start"]).union(gene_df["end"] + 1))
    segments: list[dict[str, object]] = []
    chrom = gene_df["chrom"].iloc[0]
    strand = gene_df["strand"].iloc[0]
    gene_id = gene_df["gene_id"].iloc[0]
    for order, (left, right) in enumerate(zip(breakpoints[:-1], breakpoints[1:])):
        start = int(left)
        end = int(right) - 1
        overlaps = ((gene_df["start"] <= end) & (gene_df["end"] >= start)).any()
        if not overlaps:
            continue
        segment_order = len(segments)
        segments.append(
            {
                "gene_id": gene_id,
                "segment_id": f"{gene_id}:seg{segment_order:03d}",
                "segment_order": segment_order,
                "chrom": chrom,
                "strand": strand,
                "start": start,
                "end": end,
            }
        )
    return pd.DataFrame.from_records(segments)


def _segments_in_interval(segments: pd.DataFrame, start: int, end: int, strand: str) -> list[int]:
    """Return ordered atomic segment IDs fully contained in an exon interval."""
    segment_orders = segments.loc[(segments["start"] >= start) & (segments["end"] <= end), "segment_order"].tolist()
    if strand == "-":
        segment_orders.reverse()
    return segment_orders


@dataclass
class SubisoformModel:
    """Sparse mapping from transcript isoforms to local subisoform paths."""

    transcript_ids: np.ndarray
    isoform_to_subisoform: sp.csr_matrix
    segment_table: pd.DataFrame
    subisoform_table: pd.DataFrame
    unmatched_transcripts: np.ndarray
    transcript_segment_orders: dict[str, tuple[int, ...]]

    @property
    def subisoform_ids(self) -> np.ndarray:
        """Return the ordered subisoform identifiers in this model."""
        return self.subisoform_table["subisoform_id"].to_numpy()

    def collapse(self, values: np.ndarray, transcript_ids: Optional[Sequence[object]] = None) -> np.ndarray:
        """Collapse transcript-level values into subisoform-level values."""
        matrix = np.asarray(values, dtype=float)
        if matrix.ndim == 1:
            matrix = matrix[None, :]
        if transcript_ids is not None:
            transcript_index = pd.Index(np.asarray(transcript_ids, dtype=object))
            loc = transcript_index.get_indexer(self.transcript_ids)
            if np.any(loc < 0):
                missing = self.transcript_ids[loc < 0]
                raise ValueError(f"Missing transcript values for {missing.tolist()}")
            matrix = matrix[:, loc]
        elif matrix.shape[1] != len(self.transcript_ids):
            raise ValueError("values must be aligned to model.transcript_ids when transcript_ids is omitted.")
        collapsed = sp.csr_matrix(matrix) @ self.isoform_to_subisoform
        return collapsed.toarray()


def build_subisoform_model(
    gtf_filename: str | Path,
    transcript_ids: Optional[Sequence[object]] = None,
) -> SubisoformModel:
    """Construct an anchor-based subisoform model from a GTF transcript set."""
    exons = load_transcript_exons(gtf_filename)
    if exons.empty:
        raise ValueError("No transcript exons were loaded from the GTF.")

    provided_transcripts: Optional[np.ndarray]
    transcript_order_lookup: dict[str, int]
    if transcript_ids is None:
        provided_transcripts = None
        transcript_order_lookup = {}
    else:
        provided_transcripts = np.asarray(transcript_ids, dtype=object)
        keep = exons["transcript_id"].isin(provided_transcripts)
        exons = exons.loc[keep].copy()
        transcript_order_lookup = {str(tid): idx for idx, tid in enumerate(provided_transcripts)}

    if exons.empty:
        raise ValueError("None of the requested transcripts were found in the GTF.")

    transcript_exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    transcript_gene: dict[str, str] = {}
    transcript_strand: dict[str, str] = {}
    transcript_chrom: dict[str, str] = {}
    for row in exons.itertuples(index=False):
        transcript_exons[row.transcript_id].append((int(row.start), int(row.end)))
        transcript_gene[row.transcript_id] = row.gene_id
        transcript_strand[row.transcript_id] = row.strand
        transcript_chrom[row.transcript_id] = row.chrom

    if provided_transcripts is None:
        transcript_ids_matched = np.asarray(list(transcript_exons.keys()), dtype=object)
        unmatched_transcripts = np.array([], dtype=object)
    else:
        matched_set = set(transcript_exons)
        transcript_ids_matched = np.asarray([tid for tid in provided_transcripts if tid in matched_set], dtype=object)
        unmatched_transcripts = np.asarray([tid for tid in provided_transcripts if tid not in matched_set], dtype=object)

    global_transcript_index = {str(tid): idx for idx, tid in enumerate(transcript_ids_matched)}
    transcript_segment_orders: dict[str, tuple[int, ...]] = {}

    segment_tables: list[pd.DataFrame] = []
    subisoform_records: list[dict[str, object]] = []
    event_id_lookup: dict[tuple[str, str, str], str] = {}
    subisoform_lookup: dict[tuple[str, str, tuple[int, ...]], int] = {}
    subisoform_members: dict[int, set[str]] = defaultdict(set)
    rows: list[int] = []
    cols: list[int] = []

    grouped = exons.groupby("gene_id", sort=False)
    for gene_id, gene_df in grouped:
        gene_segments = _build_atomic_segments(gene_df)
        if gene_segments.empty:
            continue

        transcript_paths: dict[str, tuple[int, ...]] = {}
        segment_transcripts: dict[int, set[str]] = defaultdict(set)
        prev_neighbors: dict[int, set[int]] = defaultdict(set)
        next_neighbors: dict[int, set[int]] = defaultdict(set)
        first_segments: set[int] = set()
        last_segments: set[int] = set()

        gene_transcripts = [tid for tid in transcript_ids_matched if transcript_gene.get(str(tid)) == gene_id]
        if not gene_transcripts:
            continue

        for transcript_id in gene_transcripts:
            merged_intervals = _merge_intervals(transcript_exons[str(transcript_id)])
            strand = transcript_strand[str(transcript_id)]
            ordered_intervals = merged_intervals if strand == "+" else list(reversed(merged_intervals))
            path: list[int] = []
            for start, end in ordered_intervals:
                path.extend(_segments_in_interval(gene_segments, start, end, strand))
            if not path:
                continue
            transcript_paths[str(transcript_id)] = tuple(path)
            transcript_segment_orders[str(transcript_id)] = tuple(path)
            first_segments.add(path[0])
            last_segments.add(path[-1])
            for segment_order in path:
                segment_transcripts[segment_order].add(str(transcript_id))
            for left, right in zip(path[:-1], path[1:]):
                prev_neighbors[right].add(left)
                next_neighbors[left].add(right)

        if not transcript_paths:
            continue

        n_transcripts = len(transcript_paths)
        anchor_mask: dict[int, bool] = {}
        for segment_order in gene_segments["segment_order"].tolist():
            if n_transcripts == 1:
                is_anchor = segment_order in first_segments or segment_order in last_segments
            else:
                is_anchor = (
                    len(segment_transcripts.get(segment_order, set())) == n_transcripts
                    or len(prev_neighbors.get(segment_order, set())) != 1
                    or len(next_neighbors.get(segment_order, set())) != 1
                    or segment_order in first_segments
                    or segment_order in last_segments
                )
            anchor_mask[segment_order] = is_anchor

        segment_tables.append(
            gene_segments.assign(
                transcript_count=[len(segment_transcripts.get(order, set())) for order in gene_segments["segment_order"]],
                is_anchor=[anchor_mask[order] for order in gene_segments["segment_order"]],
                is_terminal=[
                    (order in first_segments) or (order in last_segments) for order in gene_segments["segment_order"]
                ],
            )
        )

        event_counter = 0
        subisoform_counter = 0
        for transcript_id, path in transcript_paths.items():
            anchor_positions = sorted({0, len(path) - 1, *[i for i, order in enumerate(path) if anchor_mask[order]]})
            if len(anchor_positions) == 1:
                subpaths = [tuple(path)]
            else:
                subpaths = [tuple(path[left : right + 1]) for left, right in zip(anchor_positions[:-1], anchor_positions[1:])]

            for subpath in subpaths:
                start_order = subpath[0]
                end_order = subpath[-1]
                start_segment_id = gene_segments.loc[
                    gene_segments["segment_order"] == start_order, "segment_id"
                ].iloc[0]
                end_segment_id = gene_segments.loc[
                    gene_segments["segment_order"] == end_order, "segment_id"
                ].iloc[0]
                event_key = (gene_id, start_segment_id, end_segment_id)
                if event_key not in event_id_lookup:
                    event_id_lookup[event_key] = f"{gene_id}:evt{event_counter:03d}"
                    event_counter += 1
                event_id = event_id_lookup[event_key]

                subisoform_key = (gene_id, event_id, subpath)
                if subisoform_key not in subisoform_lookup:
                    subisoform_id = f"{gene_id}:sub{subisoform_counter:03d}"
                    subisoform_lookup[subisoform_key] = len(subisoform_records)
                    subsegments = gene_segments.set_index("segment_order").loc[list(subpath)]
                    subisoform_records.append(
                        {
                            "subisoform_id": subisoform_id,
                            "event_id": event_id,
                            "gene_id": gene_id,
                            "chrom": transcript_chrom[transcript_id],
                            "strand": transcript_strand[transcript_id],
                            "start": int(subsegments["start"].min()),
                            "end": int(subsegments["end"].max()),
                            "anchor_start_segment_id": start_segment_id,
                            "anchor_end_segment_id": end_segment_id,
                            "segment_orders": subpath,
                            "segment_ids": tuple(subsegments["segment_id"].tolist()),
                            "transcript_ids": tuple(),
                        }
                    )
                    subisoform_counter += 1
                subisoform_index = subisoform_lookup[subisoform_key]
                subisoform_members[subisoform_index].add(transcript_id)
                rows.append(global_transcript_index[transcript_id])
                cols.append(subisoform_index)

    if not subisoform_records:
        raise ValueError("No subisoforms could be constructed from the provided transcripts.")

    subisoform_table = pd.DataFrame.from_records(subisoform_records)
    subisoform_table["transcript_ids"] = [
        tuple(sorted(subisoform_members[idx], key=lambda tid: transcript_order_lookup.get(tid, 0)))
        for idx in range(len(subisoform_table))
    ]

    isoform_to_subisoform = sp.csr_matrix(
        (np.ones(len(rows), dtype=float), (rows, cols)),
        shape=(len(transcript_ids_matched), len(subisoform_table)),
    )

    return SubisoformModel(
        transcript_ids=transcript_ids_matched,
        isoform_to_subisoform=isoform_to_subisoform,
        segment_table=pd.concat(segment_tables, axis=0, ignore_index=True),
        subisoform_table=subisoform_table,
        unmatched_transcripts=unmatched_transcripts,
        transcript_segment_orders=transcript_segment_orders,
    )
