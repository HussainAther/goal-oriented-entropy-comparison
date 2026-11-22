"""Regression tests for the exploratory trajectory-error analysis."""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "code" / "compartment"),
)

from run_benchmark import N_CONFIRMATION, N_DISCOVERY
from run_state_sensitivity import (
    THRESHOLD,
    prepare_cases,
    paired_bootstrap,
)
from run_benchmark import evaluate


class TestStateSensitivity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.discovery, cls.confirmation = prepare_cases()
        cls.metrics, cls.predictions = evaluate(
            cls.discovery, cls.confirmation
        )

    def test_sample_sizes(self):
        self.assertEqual(len(self.discovery), N_DISCOVERY)
        self.assertEqual(len(self.confirmation), N_CONFIRMATION)

    def test_threshold(self):
        self.assertEqual(THRESHOLD, 0.001)

    def test_warning_prevalence(self):
        self.assertEqual(int(self.discovery.label.sum()), 112)
        self.assertEqual(int(self.confirmation.label.sum()), 86)

    def test_original_warning_preserved(self):
        for frame in (self.discovery, self.confirmation):
            self.assertIn("divergence_label", frame.columns)
            self.assertTrue(
                np.array_equal(
                    frame.label.to_numpy(),
                    (frame.state_error1 > THRESHOLD).astype(int).to_numpy(),
                )
            )

    def test_auc_ranking_reversal(self):
        results = self.metrics.set_index("method")
        self.assertGreater(
            results.loc["fisher_norm", "auc"],
            results.loc["entropy_exact", "auc"],
        )
        self.assertGreater(
            results.loc["state_norm", "auc"],
            results.loc["entropy_exact", "auc"],
        )

    def test_bootstrap_reproducibility(self):
        first = paired_bootstrap(self.predictions)
        second = paired_bootstrap(self.predictions)
        self.assertTrue(first.equals(second))
        self.assertEqual(len(first), 1500)

    def test_bootstrap_favors_fisher_norm(self):
        boot = paired_bootstrap(self.predictions)
        lo, hi = np.nanpercentile(
            boot["auc_entropy_minus_fisher_norm"],
            [2.5, 97.5],
        )
        self.assertLess(hi, 0)
        self.assertAlmostEqual(lo, -0.226206, places=5)
        self.assertAlmostEqual(hi, -0.085048, places=5)


if __name__ == "__main__":
    unittest.main()
