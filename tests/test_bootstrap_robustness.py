"""Regression tests for ecological bootstrap sensitivity."""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

sys.path.insert(
    0, str(ROOT / "code" / "comparator")
)

import run_bootstrap_robustness as br


@pytest.fixture
def predictions():
    return pd.read_csv(
        ROOT / "code/comparator/matched_predictions.csv",
        dtype={"chemostat": str},
    )


@pytest.fixture
def design(predictions):
    return br.prepare_predictions(predictions)


def test_protocol_constants():
    assert br.N_BOOT == 5000
    assert list(br.SCHEMES.values()) == [
        2026112601,
        2026112602,
        2026112603,
        2026112604,
    ]
    assert len(br.COMPARATORS) == 6


def test_complete_design(design):
    wide, grids = design

    assert len(wide) == 144
    assert len(grids) == 2

    for grid in grids.values():
        assert grid.shape == (12, 6)

    indices = np.concatenate([
        grid.ravel() for grid in grids.values()
    ])

    assert np.array_equal(
        np.sort(indices), np.arange(144)
    )


def test_duplicate_rejected(predictions):
    bad = pd.concat([
        predictions, predictions.iloc[[0]]
    ])

    with pytest.raises(ValueError, match="Duplicate"):
        br.prepare_predictions(bad)


def test_missing_prediction_rejected(predictions):
    with pytest.raises(ValueError):
        br.prepare_predictions(predictions.iloc[1:])


@pytest.mark.parametrize("scheme", list(br.SCHEMES))
def test_resample_size_and_determinism(design, scheme):
    _, grids = design

    a = br.sample_indices(
        grids, scheme, np.random.default_rng(123)
    )
    b = br.sample_indices(
        grids, scheme, np.random.default_rng(123)
    )

    assert len(a) == 144
    assert np.array_equal(a, b)
    assert np.all((a >= 0) & (a < 144))


@pytest.mark.parametrize("scheme", [
    "vessel_stratified",
    "draw_stratified",
    "crossed_stratified",
])
def test_stratified_group_counts(design, scheme):
    wide, grids = design

    idx = br.sample_indices(
        grids, scheme, np.random.default_rng(123)
    )

    counts = wide.iloc[idx].group.value_counts()

    assert sorted(counts.tolist()) == [72, 72]


def test_vessel_multiplicities(design):
    _, grids = design

    idx = br.sample_indices(
        grids,
        "vessel_stratified",
        np.random.default_rng(123),
    )

    counts = np.bincount(idx, minlength=144)

    for grid in grids.values():
        m = counts[grid]

        assert np.all(m == m[:, [0]])
        assert m[:, 0].sum() == 12


def test_draw_multiplicities(design):
    _, grids = design

    idx = br.sample_indices(
        grids,
        "draw_stratified",
        np.random.default_rng(123),
    )

    counts = np.bincount(idx, minlength=144)

    for grid in grids.values():
        m = counts[grid]

        assert np.all(m == m[[0], :])
        assert m[0, :].sum() == 6


def test_crossed_multiplicities(design):
    _, grids = design

    idx = br.sample_indices(
        grids,
        "crossed_stratified",
        np.random.default_rng(123),
    )

    counts = np.bincount(idx, minlength=144)

    for grid in grids.values():
        m = counts[grid]

        assert m.sum() == 72
        assert np.linalg.matrix_rank(m) == 1

        assert np.array_equal(
            m * 72,
            np.outer(
                m.sum(axis=1),
                m.sum(axis=0),
            ),
        )


def test_original_point_estimates(design):
    wide, _ = design

    contrasts = br.compute_contrasts(
        wide.label.to_numpy(),
        wide[list(br.METHODS)].to_numpy(),
    )

    original = pd.read_csv(
        ROOT / "code/comparator/matched_metrics.csv"
    ).set_index("method")

    for comparator in br.COMPARATORS:
        for metric in ("auc", "brier"):
            expected = (
                original.loc[br.REFERENCE, metric]
                - original.loc[comparator, metric]
            )

            actual = contrasts[
                f"{metric}_entropy_minus_{comparator}"
            ]

            assert actual == pytest.approx(
                expected, abs=1e-12
            )


def test_single_class_resample():
    y = np.zeros(8, dtype=int)

    probabilities = np.full(
        (8, len(br.METHODS)), 0.2
    )

    result = br.compute_contrasts(y, probabilities)

    for comparator in br.COMPARATORS:
        assert np.isnan(
            result[f"auc_entropy_minus_{comparator}"]
        )

        assert result[
            f"brier_entropy_minus_{comparator}"
        ] == pytest.approx(0.0)


def test_identical_predictions():
    y = np.array([0, 1, 0, 1])

    probabilities = np.tile(
        np.array([0.1, 0.9, 0.2, 0.8])[:, None],
        (1, len(br.METHODS)),
    )

    result = br.compute_contrasts(y, probabilities)

    assert all(
        value == pytest.approx(0.0)
        for value in result.values()
    )


def test_small_bootstrap_and_summary(design):
    wide, grids = design

    boot = br.run_bootstrap(wide, grids, n_boot=5)
    summary = br.summarize(wide, boot)

    assert len(boot) == 20
    assert len(summary) == 48

    assert (
        summary.valid_replicates
        + summary.undefined_replicates
    ).eq(5).all()

    assert summary[
        summary.metric == "brier"
    ].undefined_replicates.eq(0).all()

    valid = summary.ci_width.dropna()
    assert valid.ge(0).all()
