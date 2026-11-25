"""Tests for the prospective target-reversal protocol."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "code/compartment"),
)

import run_target_confirmation as tc


def artificial_predictions():
    rows = []

    for target in tc.TARGETS:
        for i in range(8):
            label = i % 2
            good = 0.9 if label else 0.1
            bad = 1.0 - good

            if target == "divergence":
                entropy, fisher = good, bad
            else:
                entropy, fisher = bad, good

            for method, probability in (
                ("entropy_exact", entropy),
                ("fisher_norm", fisher),
            ):
                rows.append({
                    "target": target,
                    "case_id": i,
                    "label": label,
                    "method": method,
                    "probability": probability,
                })

    return pd.DataFrame(rows)


def test_protocol_constants():
    assert tc.N_DISCOVERY == 192
    assert tc.N_CONFIRMATION == 1000
    assert tc.CONFIRMATION_SEED == 2026112301
    assert tc.BOOT_SEED == 2026112302
    assert tc.N_BOOT == 5000
    assert tc.TARGETS["divergence"] == ("D_error1", 0.01)
    assert tc.TARGETS["trajectory"] == ("state_error1", 0.001)


@pytest.mark.parametrize("target", list(tc.TARGETS))
def test_strict_warning_threshold(target):
    column, threshold = tc.TARGETS[target]

    frame = pd.DataFrame({
        column: [threshold / 2, threshold, threshold * 2]
    })

    result = tc.apply_target(frame, target)

    assert result.label.tolist() == [0, 0, 1]
    assert "label" not in frame.columns


def test_single_class_auc_is_undefined():
    assert np.isnan(tc.safe_auc([0, 0], [0.1, 0.2]))
    assert np.isnan(tc.safe_auc([1, 1], [0.8, 0.9]))


def test_perfect_auc():
    assert tc.safe_auc([0, 1], [0.1, 0.9]) == 1.0


def test_bootstrap_reproducibility():
    predictions = artificial_predictions()

    first = tc.paired_bootstrap(predictions, n_boot=50)
    second = tc.paired_bootstrap(predictions, n_boot=50)

    pd.testing.assert_frame_equal(first, second)


def test_paired_reversal_contrast():
    boot = tc.paired_bootstrap(
        artificial_predictions(), n_boot=100
    )

    valid = boot.dropna(subset=["auc_reversal_contrast"])

    assert len(valid) > 0

    np.testing.assert_allclose(
        valid.auc_delta_divergence, 1.0
    )
    np.testing.assert_allclose(
        valid.auc_delta_trajectory, -1.0
    )
    np.testing.assert_allclose(
        valid.auc_reversal_contrast, 2.0
    )


def test_mismatched_case_ids_rejected():
    predictions = artificial_predictions()

    predictions = predictions[
        ~(
            (predictions.target == "trajectory")
            & (predictions.case_id == 0)
        )
    ]

    with pytest.raises(ValueError, match="case IDs"):
        tc.paired_bootstrap(predictions, n_boot=5)


def test_single_class_confirmation_preserves_brier():
    discovery = pd.DataFrame({
        "case_id": [0, 1, 2, 3],
        "D_error1": [0.005, 0.02, 0.006, 0.03],
    })

    confirmation = pd.DataFrame({
        "case_id": [0, 1],
        "D_error1": [0.002, 0.004],
    })

    for feature in tc.FEATURES:
        discovery[feature] = [0.1, 0.8, 0.2, 0.9]
        confirmation[feature] = [0.1, 0.2]

    metrics, predictions, counts = tc.evaluate_target(
        discovery, confirmation, "divergence"
    )

    assert metrics.auc.isna().all()
    assert np.isfinite(metrics.brier).all()
    assert len(predictions) == 10
    assert counts["confirmation_positive"] == 0
    assert counts["confirmation_negative"] == 2
