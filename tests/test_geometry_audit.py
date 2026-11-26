"""Regression and corruption-detection tests for geometry auditing."""

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from audit_geometry_results import audit_geometry, violations


ROOT = Path(__file__).resolve().parents[1]

FILES = (
    "geometry_validation_pointwise.csv",
    "geometry_validation_cases.csv",
    "geometry_validation_summary.csv",
    "target_confirmation_cases.csv",
)


@pytest.fixture
def evidence_copy(tmp_path):
    source = ROOT / "code/compartment"
    destination = tmp_path / "code/compartment"
    destination.mkdir(parents=True)

    for name in FILES:
        shutil.copy2(source / name, destination / name)

    return tmp_path, destination


def test_complete_geometry_audit():
    result = audit_geometry(ROOT)

    assert result == {
        "pointwise_records": 40000,
        "configurations": 1000,
        "summary_records": 268,
        "pointwise_inequalities": 240000,
        "temporal_inequalities": 1000,
    }


@pytest.mark.parametrize(
    "filename,column,replacement",
    [
        (
            "geometry_validation_pointwise.csv",
            "time",
            0.123,
        ),
        (
            "geometry_validation_pointwise.csv",
            "E_loc",
            0.987,
        ),
        (
            "geometry_validation_pointwise.csv",
            "B_quad",
            0.987,
        ),
        (
            "geometry_validation_pointwise.csv",
            "condition_a",
            False,
        ),
        (
            "geometry_validation_pointwise.csv",
            "ratio_B_C",
            99.0,
        ),
        (
            "geometry_validation_pointwise.csv",
            "violation_B_S",
            True,
        ),
        (
            "geometry_validation_cases.csv",
            "mean_r",
            0.987,
        ),
        (
            "geometry_validation_cases.csv",
            "temporal_bound",
            0.987,
        ),
        (
            "geometry_validation_cases.csv",
            "divergence_label",
            -1,
        ),
        (
            "geometry_validation_cases.csv",
            "reconstruction_error_state_norm",
            1.0,
        ),
        (
            "geometry_validation_summary.csv",
            "value",
            -999.0,
        ),
        (
            "target_confirmation_cases.csv",
            "entropy_exact",
            0.987,
        ),
    ],
)
def test_corrupted_evidence_rejected(
    evidence_copy,
    filename,
    column,
    replacement,
):
    root, directory = evidence_copy
    path = directory / filename

    frame = pd.read_csv(path)
    frame.loc[0, column] = replacement
    frame.to_csv(path, index=False)

    with pytest.raises(AssertionError):
        audit_geometry(root)


def test_duplicate_pointwise_record_rejected(evidence_copy):
    root, directory = evidence_copy

    path = directory / "geometry_validation_pointwise.csv"
    frame = pd.read_csv(path)

    frame.loc[1, ["case_id", "time"]] = (
        frame.loc[0, ["case_id", "time"]].to_numpy()
    )

    frame.to_csv(path, index=False)

    with pytest.raises(
        AssertionError,
        match="Duplicate pointwise records",
    ):
        audit_geometry(root)


def test_duplicate_summary_key_rejected(evidence_copy):
    root, directory = evidence_copy

    path = directory / "geometry_validation_summary.csv"
    frame = pd.read_csv(path)

    frame.loc[1] = frame.loc[0]
    frame.to_csv(path, index=False)

    with pytest.raises(
        AssertionError,
        match="Duplicate summary keys",
    ):
        audit_geometry(root)


def test_bound_tolerance():
    bound = np.array([1.0, 1.0, 1.0])

    error = np.array([
        1.0,
        1.0 + 1e-10,
        1.0 + 1e-5,
    ])

    np.testing.assert_array_equal(
        violations(error, bound),
        [False, False, True],
    )
