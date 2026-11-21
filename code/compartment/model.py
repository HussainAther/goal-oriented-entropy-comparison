"""Positive-state nonlinear two-compartment benchmark.

The perturbation is applied to the first initial-state component.
"""

import numpy as np


def rhs(t, y, p):
    """Two-compartment saturating transfer and decay."""
    x1, x2 = y
    alpha, beta, k, K, rho = p

    transfer = k * x1 / (K + x1)

    return np.array([
        -alpha * x1 - transfer + rho * x2,
        transfer - (beta + rho) * x2,
    ])


def jacobian(y, p):
    """Jacobian of the vector field with respect to the state."""
    x1, _ = y
    alpha, beta, k, K, rho = p

    transfer_prime = k * K / (K + x1) ** 2

    return np.array([
        [-alpha - transfer_prime, rho],
        [transfer_prime, -beta - rho],
    ])


def hessian_direction(y, u, p):
    """Directional second derivative D²f(y)[u,u]."""
    x1, _ = y
    _, _, k, K, _ = p

    transfer_second = -2.0 * k * K / (K + x1) ** 3
    contribution = transfer_second * u[0] ** 2

    return np.array([-contribution, contribution])


def augmented_rhs(t, z, p):
    """Evolve the reference state and first/second variations.

    x(eps) = X + eps*u + eps**2*w + O(eps**3).
    """
    X = z[:2]
    u = z[2:4]
    w = z[4:6]

    J = jacobian(X, p)

    return np.concatenate([
        rhs(t, X, p),
        J @ u,
        J @ w + 0.5 * hessian_direction(X, u, p),
    ])
