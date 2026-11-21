"""Numerical and design checks for the synthetic compartment benchmark."""

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code" / "compartment"))

from model import augmented_rhs, hessian_direction, jacobian, rhs
from run_benchmark import (
    FEATURES,
    N_CONFIRMATION,
    N_DISCOVERY,
    N_BOOT,
    SEEDS,
    divergence,
    generate,
    integrate,
    one_case,
)


class TestCompartmentModel(unittest.TestCase):
    def setUp(self):
        self.p = (0.15, 0.20, 1.2, 2.0, 0.15)
        self.y0 = np.array([5.0, 3.0])

    def test_jacobian_against_finite_differences(self):
        h = 1e-6
        numerical = np.column_stack([
            (
                rhs(0, self.y0 + h * np.eye(2)[j], self.p)
                - rhs(0, self.y0 - h * np.eye(2)[j], self.p)
            ) / (2 * h)
            for j in range(2)
        ])
        np.testing.assert_allclose(
            jacobian(self.y0, self.p),
            numerical,
            rtol=1e-6,
            atol=1e-8,
        )

    def test_directional_hessian(self):
        u = np.array([0.8, -0.3])
        h = 1e-4
        numerical = (
            rhs(0, self.y0 + h * u, self.p)
            - 2 * rhs(0, self.y0, self.p)
            + rhs(0, self.y0 - h * u, self.p)
        ) / h**2
        np.testing.assert_allclose(
            hessian_direction(self.y0, u, self.p),
            numerical,
            rtol=1e-5,
            atol=1e-7,
        )

    def test_second_order_approximation_improves_accuracy(self):
        z0 = np.r_[self.y0, 1.0, 0.0, 0.0, 0.0]
        augmented = integrate(augmented_rhs, z0, self.p)
        X = augmented[:, :2]
        u = augmented[:, 2:4]
        w = augmented[:, 4:6]

        self.assertGreater(np.max(np.abs(w)), 1e-6)

        for eps in (0.05, 0.10, 0.20):
            actual = integrate(
                rhs, self.y0 + np.array([eps, 0.0]), self.p
            )
            first = X + eps * u
            second = first + eps**2 * w

            error_first = np.max(np.abs(actual - first))
            error_second = np.max(np.abs(actual - second))

            self.assertLess(error_second, error_first)
            self.assertTrue(np.all(actual > 0))
            self.assertTrue(np.all(second > 0))

    def test_bregman_identity(self):
        X = np.array([[5.0, 3.0], [4.0, 2.5]])
        p1 = X + np.array([[0.2, 0.1], [0.15, 0.08]])
        h = np.array([[0.01, -0.005], [0.008, -0.003]])
        p2 = p1 + h

        lhs = (
            divergence(p2, X, self.y0)
            - divergence(p1, X, self.y0)
        )
        S = np.sum(h * np.log(p1 / X) / self.y0, axis=1)
        C = divergence(p2, p1, self.y0)

        np.testing.assert_allclose(lhs, S + C, atol=1e-12)


class TestBenchmarkDesign(unittest.TestCase):
    def test_prespecified_design(self):
        self.assertEqual(N_DISCOVERY, 192)
        self.assertEqual(N_CONFIRMATION, 144)
        self.assertEqual(N_BOOT, 1500)
        self.assertNotEqual(SEEDS["discovery"], SEEDS["confirmation"])

    def test_case_reproducibility(self):
        rng1 = np.random.default_rng(SEEDS["discovery"])
        rng2 = np.random.default_rng(SEEDS["discovery"])

        a = one_case("discovery", 0, rng1)
        b = one_case("discovery", 0, rng2)

        self.assertEqual(a, b)
        self.assertIn(a["label"], (0, 1))
        self.assertLess(a["max_identity_error"], 1e-9)

    def test_distinct_case_splits(self):
        discovery = generate("discovery", 3)
        confirmation = generate("confirmation", 3)

        self.assertEqual(set(discovery["split"]), {"discovery"})
        self.assertEqual(set(confirmation["split"]), {"confirmation"})
        self.assertNotEqual(
            discovery["alpha"].tolist(),
            confirmation["alpha"].tolist(),
        )

    def test_no_duplicate_curvature_method(self):
        self.assertIn("entropy_exact", FEATURES)
        self.assertNotIn("goal_plus_curvature", FEATURES)

        case = one_case(
            "discovery",
            0,
            np.random.default_rng(SEEDS["discovery"]),
        )
        self.assertIn("entropy_exact", case)
        self.assertNotIn("goal_plus_curvature", case)


if __name__ == "__main__":
    unittest.main()
