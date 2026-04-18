# Copilot instructions

## Cluster config

This cluster uses the module system, so `module load python` is needed.

## Repository purpose

This repository analyzes Split-seq/Parse single-cell RNA-seq with an emphasis on splicing-aware quantification and downstream equivalence-class analysis.

## Build, test, and lint commands

- Build the spliceu index with `bash alevin_pipeline/1_index.sh`.
- Run the original SPLIT-seq alignment with `bash alevin_pipeline/2_align.sh`.
- Re-run quantification against an existing permit list with `bash alevin_pipeline/4_quant.sh`.
- Run the checked-in test script with `python test_sparse_cholesky.py`.
- There is no dedicated lint configuration in this repository.
- There is no test runner wrapper; the single-test command is the same `python test_sparse_cholesky.py`.

## High-level architecture

- The overall workflow is: build a splice-aware reference, quantify Split-seq/Parse reads with `salmon alevin` plus `alevin-fry`, then analyze splicing-relevant equivalence classes and cell-level matrices in Python notebooks.
- `alevin_pipeline/` contains the cluster-oriented data generation workflow. The shell scripts build a spliceu index, align Parse/SPLIT-seq reads with `salmon alevin`, then run `alevin-fry generate-permit-list`, `collate`, and `quant` to produce sparse count matrices and equivalence classes. The `README.md` in this directory explains the spliceu/splici distinction and why the microglialess dataset is often processed per SRR or per sublibrary.
- `calcutta.py` is the shared analysis library that downstream notebooks rely on. It reads FASTA, and `alevin-fry --dump-eqclasses` outputs, collapses barcode+UMI evidence into transcript or gene equivalence classes, converts per-cell results into SciPy sparse matrices, pairs Split-seq half-cells, and runs the EM routine used for transcript-level pseudobulk estimates.
- `load_fry.py` is the bridge from `alevin-fry` output into Scanpy/AnnData. It understands both standard quantification and USA mode, and the notebooks use it to build matrices for exploratory analysis.
- The notebooks are the main analysis entry points. In particular, `scanpy_alevin.ipynb`, `alevin_og.ipynb`, and `em_alevin.ipynb` compose `load_fry.py` and `calcutta.py` rather than reimplementing parsing logic.
- `get_precursors.py`, `alevin_pipeline/make_t2t.py`, and `alevin_pipeline/dedup_t2g.py` generate the precursor FASTA and custom transcript mappings that make transcript-level EC analysis possible after `alevin-fry` quantification.

## Key conventions

- Many scripts hard-code GPFS paths, Knowles Lab directories, and conda environments. Treat the repository as cluster-specific unless a task explicitly asks for portability.
- Repository code assumes `alevin-fry` outputs are central artifacts: `quant.json`/`meta_info.json`, `quants_mat.mtx`, `quants_mat_rows.txt`, `quants_mat_cols.txt`, and `gene_eqclass.txt.gz`. Reuse `load_fry.load_fry()` and `calcutta.read_alevin_ec()` before writing new parsers.
- Transcript-level analysis is encoded by swapping the usual transcript-to-gene mapping for transcript-to-transcript mappings such as `t2t.tsv` or deduplicated variants. Gene-level analysis uses `spliceu_t2g.tsv`, while USA mode uses a 3-column mapping and `load_fry(..., which_counts=...)`.
- Split-seq/Parse barcode handling is repository-specific. The top-level `README.md` and pipeline README both note that V1 and V2 barcode sets differ, the microglialess Parse data should often be processed per SRR or per sublibrary, and half-cell pairing is done downstream with `calcutta.make_cell_halfcell_matrix()` rather than by rewriting FASTQs.
- `calcutta.process_bus()` and the collapse helpers expect BUS-like text that is already grouped by barcode and UMI. Preserve that ordering assumption when generating intermediate files.
- Not every `.py` file is a pure library module. Files such as `gtf.py`, `get_precursors.py`, and `test_sparse_cholesky.py` execute work at import time, so check for top-level side effects before importing or refactoring them into shared code.
