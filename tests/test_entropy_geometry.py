"""Verification of local weighted I-divergence geometry."""

import numpy as np
import pytest


def divergence(z, x, omega):
    return np.sum(
        (z * np.log(z / x) - z + x) / omega
    )


def geometry(p, q, h, omega):
    x = p - q
    p2 = p + h

    g = 1.0 / (omega * p)

    Q = np.sum(g * q**2)
    H = np.sum(g * h**2)
    L = np.sum(g * q * h)

    S = np.sum(h * np.log(p / x) / omega)
    C = divergence(p2, p, omega)
    D = divergence(p, x, omega)

    a = np.max(np.abs(q / p))
    b = np.max(np.abs(h / p))

    bound_s = np.sum(
        np.abs(h) * q**2 / (omega * p**2)
    ) / (2 * (1 - a)**2)

    bound_c = np.sum(
        np.abs(h)**3 / (omega * p**2)
    ) / (6 * (1 - b)**2)

    bound_d = np.sum(
        np.abs(q)**3 / (omega * p**2)
    ) / (3 * (1 - a)**3)

    return {
        "Q": Q,
        "H": H,
        "L": L,
        "S": S,
        "C": C,
        "D": D,
        "bound_s": bound_s,
        "bound_c": bound_c,
        "bound_d": bound_d,
    }


def test_exact_bregman_identity():
    p = np.array([2., 3., 4.])
    q = np.array([0.2, -0.3, 0.4])
    h = np.array([-0.04, 0.06, 0.02])
    omega = np.array([1., 2., 3.])

    g = geometry(p, q, h, omega)

    exact = (
        divergence(p + h, p - q, omega)
        - divergence(p, p - q, omega)
    )

    assert g["S"] + g["C"] == pytest.approx(
        exact, abs=1e-12
    )


@pytest.mark.parametrize("seed", [123, 456, 789])
def test_taylor_remainder_bounds(seed):
    rng = np.random.default_rng(seed)

    for _ in range(100):
        p = rng.uniform(0.5, 10., size=5)
        omega = rng.uniform(0.5, 5., size=5)

        q = p * rng.uniform(-0.5, 0.5, size=5)
        h = p * rng.uniform(-0.3, 0.3, size=5)

        g = geometry(p, q, h, omega)

        rs = abs(g["S"] - g["L"])
        rc = abs(g["C"] - 0.5 * g["H"])
        rd = abs(g["D"] - 0.5 * g["Q"])

        assert rs <= g["bound_s"] + 1e-12
        assert rc <= g["bound_c"] + 1e-12
        assert rd <= g["bound_d"] + 1e-12


def test_normalized_remainder_identity():
    p = np.array([2., 3.])
    q = np.array([0.2, -0.1])
    h = np.array([0.03, 0.02])
    omega = np.array([1., 2.])

    g = geometry(p, q, h, omega)

    rs = g["S"] - g["L"]
    rc = g["C"] - 0.5 * g["H"]
    rd = g["D"] - 0.5 * g["Q"]

    exact = abs(g["S"] + g["C"]) / g["D"]

    reconstructed = abs(
        2 * g["L"] + g["H"] + 2 * rs + 2 * rc
    ) / (g["Q"] + 2 * rd)

    assert reconstructed == pytest.approx(
        exact, abs=1e-12
    )


def test_alignment_difference_bound():
    for r in np.linspace(0., 2., 101):
        for c in np.linspace(-1., 1., 101):
            entropy = abs(2 * r * c + r**2)
            alignment = r * abs(c)

            assert abs(
                entropy - 2 * alignment
            ) <= r**2 + 1e-12


def test_quadratic_ranking_counterexample():
    r_a, c_a = 0.20, 0.0
    r_b, c_b = 0.10, 1.0

    entropy_a = abs(2 * r_a * c_a + r_a**2)
    entropy_b = abs(2 * r_b * c_b + r_b**2)

    assert r_a > r_b
    assert entropy_a < entropy_b

    assert entropy_a == pytest.approx(0.04)
    assert entropy_b == pytest.approx(0.21)


def test_exact_entropy_ranking_counterexample():
    p = np.ones(2)
    omega = np.ones(2)

    q = np.array([0.01, 0.0])

    h_a = np.array([0.0, 0.002])
    h_b = np.array([0.001, 0.0])

    a = geometry(p, q, h_a, omega)
    b = geometry(p, q, h_b, omega)

    norm_a = np.sqrt(a["H"] / a["Q"])
    norm_b = np.sqrt(b["H"] / b["Q"])

    entropy_a = abs(a["S"] + a["C"]) / a["D"]
    entropy_b = abs(b["S"] + b["C"]) / b["D"]

    assert norm_a > norm_b
    assert entropy_a < entropy_b


def test_local_approximation_convergence():
    p = np.array([2., 3.])
    omega = np.array([1., 2.])

    q0 = np.array([0.4, -0.2])
    h0 = np.array([0.06, 0.03])

    errors = []

    for scale in (1., 0.5, 0.25, 0.125):
        g = geometry(
            p, scale * q0, scale * h0, omega
        )

        exact = abs(g["S"] + g["C"]) / g["D"]

        quadratic = abs(
            2 * g["L"] + g["H"]
        ) / g["Q"]

        errors.append(abs(exact - quadratic))

    assert all(
        later < earlier
        for earlier, later in zip(errors, errors[1:])
    )


@pytest.mark.parametrize("seed", [321, 654, 987])
def test_normalized_approximation_bound(seed):
    rng = np.random.default_rng(seed)

    for _ in range(100):
        p = rng.uniform(0.5, 10., size=5)
        omega = rng.uniform(0.5, 5., size=5)

        q = p * rng.uniform(-0.1, 0.1, size=5)
        h = p * rng.uniform(-0.05, 0.05, size=5)

        g = geometry(p, q, h, omega)

        Q = g["Q"]
        H = g["H"]
        L = g["L"]

        bs = g["bound_s"]
        bc = g["bound_c"]
        bd = g["bound_d"]

        assert Q > 2 * bd

        exact = abs(g["S"] + g["C"]) / g["D"]
        quadratic = abs(2 * L + H) / Q

        bound = (
            2 * (bs + bc) + 2 * quadratic * bd
        ) / (Q - 2 * bd)

        alignment = abs(L) / Q
        r_squared = H / Q

        assert abs(
            exact - quadratic
        ) <= bound + 1e-12

        assert abs(
            exact - 2 * alignment
        ) <= r_squared + bound + 1e-12
