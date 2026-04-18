from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from tealeaf.differential_usage import benjamini_hochberg, dirichlet_lrt, test_differential_subisoform_usage
from tealeaf.subisoform_simulation import estimate_lrt_operating_characteristics, simulate_dirichlet_multinomial


class DifferentialUsageTest(unittest.TestCase):
    def test_benjamini_hochberg_is_monotone(self) -> None:
        adjusted = benjamini_hochberg([0.01, 0.04, 0.03, 0.2])
        self.assertTrue(np.all((adjusted >= 0) & (adjusted <= 1)))
        self.assertLessEqual(adjusted[0], adjusted[1])

    def test_dirichlet_lrt_detects_large_shift(self) -> None:
        data = simulate_dirichlet_multinomial(
            alpha_by_condition={"a": [12, 3], "b": [3, 12]},
            n_samples_per_condition=8,
            total_count=400,
            random_state=0,
        )
        result = dirichlet_lrt(data.values, data.condition_labels)
        self.assertLess(result.p_value, 0.01)

    def test_subisoform_usage_wrapper_returns_event(self) -> None:
        values = np.array(
            [
                [80, 20, 0],
                [75, 25, 0],
                [20, 80, 0],
                [25, 75, 0],
            ],
            dtype=float,
        )
        table = pd.DataFrame(
            {
                "subisoform_id": ["s1", "s2", "s3"],
                "event_id": ["e1", "e1", "e2"],
                "gene_id": ["g1", "g1", "g2"],
            }
        )
        result = test_differential_subisoform_usage(
            subisoform_values=values,
            subisoform_table=table,
            condition_labels=np.array(["a", "a", "b", "b"], dtype=object),
        )
        self.assertEqual(result["event_id"].tolist(), ["e1"])
        self.assertLess(result["p_value"].iloc[0], 0.05)

    def test_lrt_operating_characteristics_show_power_over_type1(self) -> None:
        metrics = estimate_lrt_operating_characteristics(
            alpha_null=[8, 8],
            alpha_alt_a=[14, 2],
            alpha_alt_b=[2, 14],
            n_samples_per_condition=6,
            total_count=300,
            n_reps=200,
            random_state=0,
        )
        self.assertLess(metrics.type1_error, 0.13)
        self.assertGreater(metrics.power, 0.75)
        self.assertGreater(metrics.power, metrics.type1_error + 0.5)


if __name__ == "__main__":
    unittest.main()
