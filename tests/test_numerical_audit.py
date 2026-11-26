import numpy as np
import pytest

from audit_numerical_results import (
    independent_auc,
    independent_brier,
    independent_log_loss,
    percentile_interval,
)


def test_perfect_auc():
    assert independent_auc(
        [0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9]
    ) == 1.0


def test_reversed_auc():
    assert independent_auc(
        [0, 0, 1, 1], [0.9, 0.8, 0.2, 0.1]
    ) == 0.0


def test_all_tied_auc():
    assert independent_auc(
        [0, 1, 0, 1], [0.5] * 4
    ) == 0.5


def test_partial_ties():
    assert independent_auc(
        [0, 0, 1, 1], [0.1, 0.5, 0.5, 0.9]
    ) == 0.875


def test_perfect_brier():
    assert independent_brier(
        [0, 1], [0.0, 1.0]
    ) == 0.0


def test_constant_brier():
    assert independent_brier(
        [0, 1], [0.5, 0.5]
    ) == 0.25


def test_log_loss():
    assert independent_log_loss(
        [0, 1], [0.5, 0.5]
    ) == pytest.approx(np.log(2))


def test_percentile_interval():
    np.testing.assert_allclose(
        percentile_interval(np.arange(101)),
        [2.5, 97.5],
    )


@pytest.mark.parametrize(
    "probabilities",
    [[-0.1, 0.5], [0.5, 1.1], [np.nan, 0.5]],
)
def test_invalid_probabilities(probabilities):
    with pytest.raises(AssertionError):
        independent_auc([0, 1], probabilities)


def test_single_class_auc():
    with pytest.raises(AssertionError):
        independent_auc([1, 1], [0.3, 0.8])


def test_invalid_labels():
    with pytest.raises(AssertionError):
        independent_brier([0, 2], [0.3, 0.8])
