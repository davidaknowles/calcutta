from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

import calcutta


DEFAULT_DATA_DIR = Path("/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice")
DEFAULT_PARSE_META_DIR = Path("/gpfs/commons/groups/knowles_lab/data/sc/splitpool/Parse_meta")
DEFAULT_INDEX_DIR = Path("/gpfs/commons/groups/knowles_lab/index/salmon/mus_spliceu")


@dataclass
class MicroglialessInputs:
    cell_ec_matrix: sp.csr_matrix
    ec_transcript_mat: sp.csr_matrix
    transcript_ids: np.ndarray
    cell_metadata: pd.DataFrame
    selected_gene_ids: np.ndarray
    selected_sublibraries: np.ndarray


def _load_sparse_matrix(mtx_path: Path, cache_path: Optional[Path] = None) -> sp.csr_matrix:
    if cache_path is None:
        cache_path = mtx_path.with_suffix(".npz")
    if cache_path.exists():
        return scipy.sparse.load_npz(cache_path).tocsr()
    matrix = scipy.io.mmread(mtx_path).tocsr()
    scipy.sparse.save_npz(cache_path, matrix)
    return matrix


def load_preprocessed_cell_metadata(data_dir: Path = DEFAULT_DATA_DIR) -> pd.DataFrame:
    barcode_ids = [line.strip() for line in open(data_dir / "GSM5693472_barcodes.tsv", "r")]
    barcode_prefix = [value.split("_")[0] for value in barcode_ids]
    metadata = pd.read_csv(data_dir / "GSM5693472_cell_metadata.txt.gz", sep="\t")
    metadata["rnd1_well"] = [int(value.split("_")[1]) for value in barcode_ids]
    metadata["rnd3_BC"] = [value[:8] for value in barcode_prefix]
    metadata["rnd2_BC"] = [value[8:] for value in barcode_prefix]
    return metadata


def load_transcript_gene_map(index_dir: Path = DEFAULT_INDEX_DIR) -> pd.DataFrame:
    mapping = pd.read_csv(index_dir / "spliceu_t2g.tsv", sep="\t", names=["transcript_id", "gene_id"])
    mapping["gene_stable_id"] = mapping["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    return mapping


def select_transcripts_for_genes(
    transcript_ids: Sequence[object],
    gene_ids: Sequence[object],
    index_dir: Path = DEFAULT_INDEX_DIR,
) -> np.ndarray:
    transcript_ids = np.asarray(transcript_ids, dtype=object)
    gene_ids = np.asarray(gene_ids, dtype=object)
    mapping = load_transcript_gene_map(index_dir=index_dir)
    selected = mapping.loc[mapping["gene_id"].isin(gene_ids), "transcript_id"]
    return np.asarray([tid for tid in transcript_ids if tid in set(selected)], dtype=object)


def _prepare_parse_barcode_tables(parse_meta_dir: Path = DEFAULT_PARSE_META_DIR) -> tuple[pd.DataFrame, dict[str, int]]:
    ligation = pd.read_csv(parse_meta_dir / "Parse_ligation_barcodes.txt", sep="\t", header=None)[0].astype(str).tolist()
    rt = pd.read_csv(parse_meta_dir / "Parse_RT_barcodes.txt", sep="\t", header=None)[0].astype(str).tolist()
    polydt_hex_pairs = pd.DataFrame(
        {"rnd1_BC": rt, "rnd1_well": np.tile(np.arange(48), 2), "polydT": [True] * 48 + [False] * 48}
    )
    ligation_lookup = dict(zip(ligation, np.arange(96)))
    return polydt_hex_pairs, ligation_lookup


def match_feature_dump_to_preprocessed(
    feature_dump: pd.DataFrame,
    sublibrary: str,
    preprocessed_metadata: pd.DataFrame,
    polydt_hex_pairs: pd.DataFrame,
    ligation_lookup: dict[str, int],
) -> pd.DataFrame:
    cells = feature_dump.copy()
    cells["local_cell_index"] = np.arange(len(cells))
    cells["sublibrary"] = sublibrary
    cells["rnd3_BC"] = cells["CB"].str[:8]
    cells["rnd2_BC"] = cells["CB"].str[8:16]
    cells["rnd1_BC"] = cells["CB"].str[16:]
    cells["rnd2_well"] = cells["rnd2_BC"].map(ligation_lookup)
    cells["rnd3_well"] = cells["rnd3_BC"].map(ligation_lookup)
    cells = pd.merge(cells, polydt_hex_pairs, on="rnd1_BC")
    return pd.merge(
        cells,
        preprocessed_metadata,
        on=["sublibrary", "rnd1_well", "rnd2_well", "rnd3_well"],
        how="inner",
    )


def _load_eqclass_subset(
    fry_dir: Path,
    selected_transcript_ids: np.ndarray,
    selected_row_indices: np.ndarray,
) -> tuple[sp.csr_matrix, list[tuple[int, ...]]]:
    feature_ids = np.asarray([line.strip() for line in open(fry_dir / "alevin" / "quants_mat_cols.txt", "r")], dtype=object)
    selected_lookup = {tid: idx for idx, tid in enumerate(selected_transcript_ids)}
    feature_to_selected = np.full(len(feature_ids), -1, dtype=int)
    for feature_index, transcript_id in enumerate(feature_ids):
        if transcript_id in selected_lookup:
            feature_to_selected[feature_index] = selected_lookup[transcript_id]

    _, _, ecs = calcutta.read_alevin_ec(fry_dir / "alevin" / "gene_eqclass.txt.gz")
    local_cols: list[int] = []
    local_keys: list[tuple[int, ...]] = []
    for ec_index in range(len(ecs)):
        transcript_subset = feature_to_selected[np.asarray(ecs[ec_index], dtype=int)]
        transcript_subset = transcript_subset[transcript_subset >= 0]
        if transcript_subset.size == 0:
            continue
        local_cols.append(ec_index)
        local_keys.append(tuple(sorted(set(int(value) for value in transcript_subset))))

    if not local_cols:
        return sp.csr_matrix((len(selected_row_indices), 0)), []

    geqc = _load_sparse_matrix(fry_dir / "alevin" / "geqc_counts.mtx")
    subset = geqc[selected_row_indices, :][:, local_cols].tocsr()
    return subset, local_keys


def load_all_sublib_inputs(
    selected_gene_ids: Sequence[object],
    classes_to_keep: Optional[Sequence[object]] = None,
    data_dir: Path = DEFAULT_DATA_DIR,
    parse_meta_dir: Path = DEFAULT_PARSE_META_DIR,
    index_dir: Path = DEFAULT_INDEX_DIR,
    quant_name: str = "quant_spliceu_t2t",
    sublibraries: Optional[Sequence[str]] = None,
) -> MicroglialessInputs:
    polydt_hex_pairs, ligation_lookup = _prepare_parse_barcode_tables(parse_meta_dir=parse_meta_dir)
    preprocessed_metadata = load_preprocessed_cell_metadata(data_dir=data_dir)

    all_sublib_dir = data_dir / "salmon_per_sublib"
    if sublibraries is None:
        sublibraries = sorted(
            [
                path.name
                for path in all_sublib_dir.iterdir()
                if path.is_dir() and path.name != "merge" and (path / "out_permit_known" / quant_name).exists()
            ]
        )

    first_dir = all_sublib_dir / sublibraries[0] / "out_permit_known" / quant_name
    transcript_ids_all = np.asarray(
        [line.strip() for line in open(first_dir / "alevin" / "quants_mat_cols.txt", "r")],
        dtype=object,
    )
    selected_gene_ids = np.asarray(selected_gene_ids, dtype=object)
    selected_transcript_ids = select_transcripts_for_genes(
        transcript_ids=transcript_ids_all,
        gene_ids=selected_gene_ids,
        index_dir=index_dir,
    )
    if selected_transcript_ids.size == 0:
        raise ValueError("No transcripts from the selected genes were found in the transcriptome index.")

    global_ec_lookup: dict[tuple[int, ...], int] = {}
    combined_matrices: list[sp.csr_matrix] = []
    combined_metadata: list[pd.DataFrame] = []

    for sublibrary in sublibraries:
        fry_dir = all_sublib_dir / sublibrary / "out_permit_known" / quant_name
        feature_dump = pd.read_csv(fry_dir / "featureDump.txt", sep="\t")
        matched = match_feature_dump_to_preprocessed(
            feature_dump=feature_dump,
            sublibrary=sublibrary,
            preprocessed_metadata=preprocessed_metadata,
            polydt_hex_pairs=polydt_hex_pairs,
            ligation_lookup=ligation_lookup,
        )
        if classes_to_keep is not None:
            matched = matched.loc[matched["class"].isin(classes_to_keep)].copy()
        if matched.empty:
            continue

        matched["analysis_cell_id"] = [f"{sublibrary}_{idx}" for idx in matched["local_cell_index"]]
        local_matrix, local_keys = _load_eqclass_subset(
            fry_dir=fry_dir,
            selected_transcript_ids=selected_transcript_ids,
            selected_row_indices=matched["local_cell_index"].to_numpy(),
        )
        if local_matrix.shape[1] == 0:
            continue

        mapped_cols = []
        for key in local_keys:
            if key not in global_ec_lookup:
                global_ec_lookup[key] = len(global_ec_lookup)
            mapped_cols.append(global_ec_lookup[key])
        local_coo = local_matrix.tocoo()
        remapped = sp.csr_matrix(
            (local_coo.data, (local_coo.row, np.asarray(mapped_cols, dtype=int)[local_coo.col])),
            shape=(local_matrix.shape[0], len(global_ec_lookup)),
        )
        combined_matrices.append(remapped)
        combined_metadata.append(matched.reset_index(drop=True))

    if not combined_matrices:
        raise ValueError("No matched cells were found for the requested classes and genes.")

    global_width = len(global_ec_lookup)
    padded_matrices = [
        matrix if matrix.shape[1] == global_width else sp.hstack([matrix, sp.csr_matrix((matrix.shape[0], global_width - matrix.shape[1]))]).tocsr()
        for matrix in combined_matrices
    ]
    cell_ec_matrix = sp.vstack(padded_matrices).tocsr()
    cell_metadata = pd.concat(combined_metadata, axis=0, ignore_index=True)

    ec_transcript_mat = sp.csr_matrix(
        (
            np.ones(sum(len(key) for key in global_ec_lookup), dtype=float),
            (
                np.concatenate([[global_ec_lookup[key]] * len(key) for key in global_ec_lookup]),
                np.concatenate([np.asarray(key, dtype=int) for key in global_ec_lookup]),
            ),
        ),
        shape=(len(global_ec_lookup), len(selected_transcript_ids)),
    )

    return MicroglialessInputs(
        cell_ec_matrix=cell_ec_matrix,
        ec_transcript_mat=ec_transcript_mat,
        transcript_ids=selected_transcript_ids,
        cell_metadata=cell_metadata,
        selected_gene_ids=selected_gene_ids,
        selected_sublibraries=np.asarray(sublibraries, dtype=object),
    )
