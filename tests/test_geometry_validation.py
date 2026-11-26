"""Regression tests for retrospective geometry assessment."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(
    0,
    str(
        Path(__file__).resolve().parents[1]
        / "code" / "compartment"
    ),
)

from run_geometry_validation import (
    pointwise_geometry,
    safe_divide,
    bound_violation,
)


def example(h):
    return pointwise_geometry(
        np.array([[0.99, 1.0]]),
        np.array([[1.0, 1.0]]),
        np.array([h]),
        np.ones(2),
    ).iloc[0]


def test_orthogonal_correction():
    row = example([0.0, 0.002])

    assert row.r == pytest.approx(0.2)
    assert row.c == pytest.approx(0.0)
    assert row.E_quad == pytest.approx(0.04)
    assert row.twice_alignment == pytest.approx(0.0)


def test_parallel_correction():
    row = example([0.001, 0.0])

    assert row.r == pytest.approx(0.1)
    assert row.c == pytest.approx(1.0)
    assert row.E_quad == pytest.approx(0.21)
    assert row.twice_alignment == pytest.approx(0.2)


def test_antiparallel_correction():
    row = example([-0.001, 0.0])

    assert row.c == pytest.approx(-1.0)
    assert row.E_quad == pytest.approx(0.19)


def test_exact_ranking_reversal():
    a = example([0.0, 0.002])
    b = example([0.001, 0.0])

    assert a.r > b.r
    assert a.E_loc < b.E_loc


def test_zero_correction():
    row = example([0.0, 0.0])

    assert row.r == 0
    assert np.isnan(row.c)
    assert row.twice_alignment == 0
    assert row.E_quad == 0
    assert row.E_loc == pytest.approx(0)


def test_zero_displacement():
    row = pointwise_geometry(
        np.ones((1, 2)),
        np.ones((1, 2)),
        np.array([[0.001, 0.0]]),
        np.ones(2),
    ).iloc[0]

    assert not row.condition_Q
    assert not row.applicable
    assert np.isnan(row.r)
    assert np.isnan(row.E_quad)


def test_failed_relative_displacement_condition():
    row = pointwise_geometry(
        np.array([[3.0, 1.0]]),
        np.ones((1, 2)),
        np.array([[0.001, 0.0]]),
        np.ones(2),
    ).iloc[0]

    assert not row.condition_a
    assert not row.applicable
    assert np.isnan(row.B_D)


def test_inactive_floor_requirement():
    row = pointwise_geometry(
        np.array([[1.0 - 1e-8, 1.0]]),
        np.ones((1, 2)),
        np.array([[1e-9, 0.0]]),
        np.ones(2),
    ).iloc[0]

    assert not row.condition_alignment_floor
    assert not row.applicable


def test_nonpositive_state():
    with pytest.raises(ValueError):
        pointwise_geometry(
            np.array([[0.0, 1.0]]),
            np.ones((1, 2)),
            np.zeros((1, 2)),
            np.ones(2),
        )


def test_nonfinite_input():
    with pytest.raises(ValueError):
        pointwise_geometry(
            np.ones((1, 2)),
            np.ones((1, 2)),
            np.array([[np.nan, 0.0]]),
            np.ones(2),
        )


@pytest.mark.parametrize(
    "h",
    [
        [0.0, 0.002],
        [0.001, 0.0],
        [-0.001, 0.0],
        [0.001, -0.002],
    ],
)
def test_analytical_bounds(h):
    row = example(h)

    assert row.applicable

    for bound in (
        "B_S", "B_C", "B_D",
        "B_quad", "B_r2", "B_alignment",
    ):
        assert not row[f"violation_{bound}"]


def test_safe_division():
    result = safe_divide([2, 1, 1], [4, 0, -1])

    assert result[0] == pytest.approx(0.5)
    assert np.isnan(result[1:]).all()


def test_bound_violation_detection():
    result = bound_violation(
        [0.5, 1.5, 1.0, np.nan],
        [1.0, 1.0, 1.0, 1.0],
    )

    np.testing.assert_array_equal(
        result, [False, True, False, False]
    )
