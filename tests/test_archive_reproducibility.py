"""Regression tests for archive integrity and numerical comparison."""

import hashlib
from pathlib import Path

import pandas as pd
import pytest

from verify_archive import MANIFESTS, verify_archives
from reproduce_ecology import compare_predictions


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def archive(tmp_path):
    payload = tmp_path / "evidence.csv"
    payload.write_text("value\n1\n")

    digest = hashlib.sha256(
        payload.read_bytes()
    ).hexdigest()

    for name in MANIFESTS:
        manifest = tmp_path / name
        manifest.parent.mkdir(
            parents=True, exist_ok=True
        )
        manifest.write_text(
            f"{digest}  evidence.csv\n"
        )

    return tmp_path, payload


@pytest.fixture
def predictions(tmp_path):
    frame = pd.DataFrame({
        "group": ["mono", "mono"],
        "chemostat": ["001", "002"],
        "draw": [1, 1],
        "method": ["entropy_exact", "entropy_exact"],
        "label": [0, 1],
        "probability": [0.2, 0.8],
    })

    a = tmp_path / "archived.csv"
    b = tmp_path / "regenerated.csv"

    frame.to_csv(a, index=False)
    frame.to_csv(b, index=False)

    return frame, a, b


def test_repository_archives():
    assert len(verify_archives(ROOT)) == 12


def test_geometry_archives_registered():
    verified = verify_archives(ROOT)

    expected = {
        "code/compartment/geometry_validation_cases.csv",
        "code/compartment/geometry_validation_pointwise.csv",
        "code/compartment/geometry_validation_summary.csv",
    }

    assert expected <= set(verified)


def test_duplicate_manifest_entries_counted_once(archive):
    root, _ = archive
    assert len(verify_archives(root)) == 1


def test_valid_archive(archive):
    root, _ = archive
    assert "evidence.csv" in verify_archives(root)


def test_corrupted_archive(archive):
    root, payload = archive
    payload.write_text("value\n2\n")

    with pytest.raises(ValueError, match="Checksum mismatch"):
        verify_archives(root)


def test_missing_archive(archive):
    root, payload = archive
    payload.unlink()

    with pytest.raises(FileNotFoundError):
        verify_archives(root)


def test_identical_predictions(predictions):
    _, a, b = predictions
    result = compare_predictions(a, b)

    assert result["prediction_rows"] == 2
    assert result["max_absolute_probability_difference"] == 0


def test_prediction_row_order(predictions):
    frame, a, b = predictions
    frame.iloc[::-1].to_csv(b, index=False)

    compare_predictions(a, b)


def test_small_numerical_difference(predictions):
    frame, a, b = predictions
    frame["probability"] += 1e-9
    frame.to_csv(b, index=False)

    result = compare_predictions(a, b)

    assert result["max_absolute_probability_difference"] < 1e-8


def test_large_numerical_difference(predictions):
    frame, a, b = predictions
    frame["probability"] += 0.01
    frame.to_csv(b, index=False)

    with pytest.raises(AssertionError):
        compare_predictions(a, b)


def test_changed_labels(predictions):
    frame, a, b = predictions
    frame["label"] = 1 - frame["label"]
    frame.to_csv(b, index=False)

    with pytest.raises(AssertionError):
        compare_predictions(a, b)


def test_duplicate_predictions(predictions):
    frame, a, b = predictions
    pd.concat([frame, frame]).to_csv(b, index=False)

    with pytest.raises(ValueError, match="Duplicate"):
        compare_predictions(a, b)


def test_missing_predictions(predictions):
    frame, a, b = predictions
    frame.iloc[:1].to_csv(b, index=False)

    with pytest.raises(AssertionError):
        compare_predictions(a, b)
