from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from tealeaf.subisoform_pipeline import run_subisoform_pipeline
from tealeaf.subisoform_simulation import build_gtf_event_simulation_design, run_end_to_end_simulation
from tealeaf.subisoforms import build_subisoform_model


SYNTHETIC_GTF = """\
chr1\ttest\texon\t100\t199\t.\t+\t.\tgene_id "gene1"; transcript_id "t1";
chr1\ttest\texon\t300\t349\t.\t+\t.\tgene_id "gene1"; transcript_id "t1";
chr1\ttest\texon\t400\t499\t.\t+\t.\tgene_id "gene1"; transcript_id "t1";
chr1\ttest\texon\t100\t199\t.\t+\t.\tgene_id "gene1"; transcript_id "t2";
chr1\ttest\texon\t400\t499\t.\t+\t.\tgene_id "gene1"; transcript_id "t2";
chr1\ttest\texon\t50\t199\t.\t+\t.\tgene_id "gene1"; transcript_id "t3";
chr1\ttest\texon\t400\t499\t.\t+\t.\tgene_id "gene1"; transcript_id "t3";
"""

REAL_GTF_EVENT = """\
chr6\tHAVANA\texon\t128803122\t128803215\t.\t-\t.\tgene_id "ENSMUSG00000107872.3"; transcript_id "ENSMUST00000204394.3";
chr6\tHAVANA\texon\t128797905\t128798003\t.\t-\t.\tgene_id "ENSMUSG00000107872.3"; transcript_id "ENSMUST00000204394.3";
chr6\tHAVANA\texon\t128797228\t128797299\t.\t-\t.\tgene_id "ENSMUSG00000107872.3"; transcript_id "ENSMUST00000204394.3";
chr6\tHAVANA\texon\t128761068\t128761222\t.\t-\t.\tgene_id "ENSMUSG00000107872.3"; transcript_id "ENSMUST00000204394.3";
chr6\tHAVANA\texon\t128760627\t128760742\t.\t-\t.\tgene_id "ENSMUSG00000107872.3"; transcript_id "ENSMUST00000204394.3";
chr6\tHAVANA\texon\t128757248\t128757389\t.\t-\t.\tgene_id "ENSMUSG00000107872.3"; transcript_id "ENSMUST00000204394.3";
chr6\tHAVANA\texon\t128803122\t128803215\t.\t-\t.\tgene_id "ENSMUSG00000107872.3"; transcript_id "ENSMUST00000204423.3";
chr6\tHAVANA\texon\t128797905\t128798003\t.\t-\t.\tgene_id "ENSMUSG00000107872.3"; transcript_id "ENSMUST00000204423.3";
chr6\tHAVANA\texon\t128797228\t128797299\t.\t-\t.\tgene_id "ENSMUSG00000107872.3"; transcript_id "ENSMUST00000204423.3";
chr6\tHAVANA\texon\t128795876\t128796030\t.\t-\t.\tgene_id "ENSMUSG00000107872.3"; transcript_id "ENSMUST00000204423.3";
chr6\tHAVANA\texon\t128795436\t128795551\t.\t-\t.\tgene_id "ENSMUSG00000107872.3"; transcript_id "ENSMUST00000204423.3";
chr6\tHAVANA\texon\t128757248\t128757389\t.\t-\t.\tgene_id "ENSMUSG00000107872.3"; transcript_id "ENSMUST00000204423.3";
"""


class SubisoformPipelineTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.gtf_path = Path(self.tempdir.name) / "synthetic.gtf"
        self.gtf_path.write_text(SYNTHETIC_GTF)
        self.transcript_ids = np.array(["t1", "t2", "t3"], dtype=object)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_build_subisoform_model_finds_alternative_event_groups(self) -> None:
        model = build_subisoform_model(self.gtf_path, transcript_ids=self.transcript_ids)
        self.assertEqual(model.unmatched_transcripts.size, 0)
        event_sizes = model.subisoform_table.groupby("event_id").size()
        self.assertTrue((event_sizes >= 2).any())
        self.assertGreater(model.isoform_to_subisoform.nnz, len(self.transcript_ids))

    def test_end_to_end_pipeline_detects_differential_usage(self) -> None:
        compatibility = sp.csr_matrix(
            np.array(
                [
                    [1, 0, 0],  # t1 unique
                    [0, 1, 0],  # t2 unique
                    [1, 1, 0],  # t1/t2 ambiguous
                    [0, 0, 1],  # t3 unique
                    [1, 0, 1],  # t1/t3 ambiguous
                ],
                dtype=float,
            )
        )

        cell_ec = sp.csr_matrix(
            np.array(
                [
                    [40, 2, 18, 5, 2],
                    [36, 4, 17, 6, 3],
                    [42, 1, 20, 4, 3],
                    [38, 3, 19, 5, 2],
                    [41, 2, 16, 6, 2],
                    [39, 3, 18, 5, 3],
                    [2, 40, 18, 5, 2],
                    [4, 36, 17, 6, 3],
                    [1, 42, 20, 4, 3],
                    [3, 38, 19, 5, 2],
                    [2, 41, 16, 6, 2],
                    [3, 39, 18, 5, 3],
                ],
                dtype=float,
            )
        )
        group_labels = np.array(["type_a"] * 6 + ["type_b"] * 6, dtype=object)

        result = run_subisoform_pipeline(
            cell_ec_matrix=cell_ec,
            ec_transcript_mat=compatibility,
            transcript_ids=self.transcript_ids,
            transcript_weights=np.ones(len(self.transcript_ids), dtype=float),
            group_labels=group_labels,
            gtf_filename=self.gtf_path,
            metacell_size=2,
            random_state=0,
            em_iterations=200,
            em_tol=1e-10,
            min_samples_per_condition=2,
        )

        self.assertEqual(result.grouped_em.sample_metadata["condition"].nunique(), 2)
        self.assertEqual(result.subisoform_tpm.shape[0], len(result.grouped_em.sample_ids))
        self.assertFalse(result.differential_usage.empty)
        self.assertLess(result.differential_usage["p_value"].min(), 0.05)

    def test_real_gtf_event_simulation_runs_full_pipeline(self) -> None:
        real_gtf_path = Path(self.tempdir.name) / "real_event.gtf"
        real_gtf_path.write_text(REAL_GTF_EVENT)
        transcript_ids = np.array(["ENSMUST00000204394.3", "ENSMUST00000204423.3"], dtype=object)

        design = build_gtf_event_simulation_design(real_gtf_path, transcript_ids=transcript_ids)
        self.assertEqual(design.target_event_id, "ENSMUSG00000107872.3:evt002")
        self.assertEqual(design.ec_transcript_mat.shape[1], len(transcript_ids))

        _, pipeline = run_end_to_end_simulation(
            design=design,
            alpha_by_condition={"a": [18, 2], "b": [2, 18]},
            n_cells_per_condition=18,
            mean_gene_count=120,
            metacell_size=3,
            random_state=7,
        )

        self.assertFalse(pipeline.differential_usage.empty)
        detected = pipeline.differential_usage.set_index("event_id")
        self.assertIn(design.target_event_id, detected.index)
        self.assertLess(detected.loc[design.target_event_id, "p_value"], 0.05)


if __name__ == "__main__":
    unittest.main()
