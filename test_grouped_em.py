from __future__ import annotations

import unittest

import numpy as np
import scipy.sparse as sp

from tealeaf.grouped_em import aggregate_rows, make_metacell_assignments, run_grouped_em
from tealeaf.transcript_utils import resolve_transcript_lengths


class GroupedEmTest(unittest.TestCase):
    def test_aggregate_rows_preserves_label_order(self) -> None:
        matrix = sp.csr_matrix(np.array([[1, 0], [2, 1], [0, 3], [4, 0]], dtype=float))
        labels = np.array(["b", "a", "b", "a"], dtype=object)
        uniques, aggregated, counts = aggregate_rows(matrix, labels)
        self.assertEqual(uniques.tolist(), ["b", "a"])
        np.testing.assert_array_equal(counts, np.array([2, 2]))
        np.testing.assert_array_equal(aggregated.toarray(), np.array([[1, 3], [6, 1]], dtype=float))

    def test_make_metacell_assignments_respects_group_and_batch(self) -> None:
        groups = np.array(["x"] * 5 + ["y"] * 4, dtype=object)
        batches = np.array(["b1", "b1", "b2", "b2", "b2", "b1", "b1", "b2", "b2"], dtype=object)
        assignments = make_metacell_assignments(groups, batch_labels=batches, target_size=2, random_state=0)
        self.assertEqual(len(assignments), len(groups))
        grouped = assignments.groupby("sample_id").agg(condition=("condition", "nunique"), batch=("batch", "nunique"))
        self.assertTrue((grouped["condition"] == 1).all())
        self.assertTrue((grouped["batch"] == 1).all())
        self.assertTrue((assignments.groupby("sample_id").size() <= 2).all())

    def test_run_grouped_em_separates_groups(self) -> None:
        compatibility = sp.csr_matrix(
            np.array(
                [
                    [1, 0],
                    [0, 1],
                    [1, 1],
                ],
                dtype=float,
            )
        )
        cell_ec = sp.csr_matrix(
            np.array(
                [
                    [40, 2, 10],
                    [36, 4, 10],
                    [3, 35, 10],
                    [2, 38, 11],
                ],
                dtype=float,
            )
        )
        result = run_grouped_em(
            cell_ec_matrix=cell_ec,
            ec_transcript_mat=compatibility,
            transcript_weights=np.ones(2, dtype=float),
            transcript_ids=np.array(["t1", "t2"], dtype=object),
            group_labels=np.array(["a", "a", "b", "b"], dtype=object),
            iterations=200,
            tol=1e-10,
        )
        tpm = result.tpm_frame()
        self.assertGreater(tpm.loc["a", "t1"], tpm.loc["a", "t2"])
        self.assertGreater(tpm.loc["b", "t2"], tpm.loc["b", "t1"])

    def test_resolve_transcript_lengths_uses_alias_map(self) -> None:
        lengths, missing = resolve_transcript_lengths(
            transcript_ids=["tx1", "tx2"],
            transcript_lengths={"tx1": 1000, "alias_tx2": 750},
            alias_map={"tx2": "alias_tx2"},
        )
        np.testing.assert_array_equal(lengths, np.array([1000.0, 750.0]))
        self.assertEqual(missing.size, 0)


if __name__ == "__main__":
    unittest.main()
