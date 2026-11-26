"""Verify the auxiliary Poisson information geometry."""

import numpy as np
import pytest

from scipy.stats import poisson


def weighted_divergence(z, x, omega):
    return np.sum(
        (z * np.log(z / x) - z + x) / omega
    )


@pytest.mark.parametrize("seed", [123, 456, 789])
def test_poisson_kl_identity(seed):
    rng = np.random.default_rng(seed)

    for _ in range(30):
        x = rng.uniform(0.5, 5.0, size=4)
        z = rng.uniform(0.5, 5.0, size=4)
        omega = rng.uniform(0.5, 2.0, size=4)

        expected = weighted_divergence(z, x, omega)

        numerical = 0.0
        support = np.arange(150)

        for zi, xi, wi in zip(z, x, omega):
            mu_z = zi / wi
            mu_x = xi / wi

            log_p = poisson.logpmf(support, mu_z)
            log_q = poisson.logpmf(support, mu_x)

            numerical += np.sum(
                np.exp(log_p) * (log_p - log_q)
            )

        np.testing.assert_allclose(
            numerical,
            expected,
            rtol=1e-10,
            atol=1e-12,
        )


@pytest.mark.parametrize("seed", [321, 654, 987])
def test_poisson_fisher_identity(seed):
    rng = np.random.default_rng(seed)

    for _ in range(30):
        x = rng.uniform(0.5, 5.0, size=4)
        omega = rng.uniform(0.5, 2.0, size=4)

        support = np.arange(150)

        expected = 1.0 / (omega * x)

        numerical = []

        for xi, wi in zip(x, omega):
            probabilities = poisson.pmf(
                support, xi / wi
            )

            score = support / xi - 1.0 / wi

            np.testing.assert_allclose(
                np.sum(probabilities * score),
                0.0,
                atol=1e-12,
            )

            numerical.append(
                np.sum(probabilities * score**2)
            )

        np.testing.assert_allclose(
            numerical,
            expected,
            rtol=1e-10,
            atol=1e-12,
        )


def test_weighted_divergence_local_metric():
    x = np.array([1.0, 2.0, 3.0])
    omega = np.array([0.5, 1.5, 2.0])
    v = np.array([0.2, -0.3, 0.4])

    eps = 1e-3

    exact = weighted_divergence(
        x + eps * v, x, omega
    )

    quadratic = (
        0.5 * eps**2
        * np.sum(v**2 / (omega * x))
    )

    np.testing.assert_allclose(
        exact,
        quadratic,
        rtol=1e-3,
        atol=1e-12,
    )
