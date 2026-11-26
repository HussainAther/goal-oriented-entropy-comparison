"""Regenerate ecological results without modifying archived evidence."""

from pathlib import Path
from tempfile import TemporaryDirectory
import json
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

from verify_archive import verify_archives


STAGES = (
    "perturbation",
    "features",
    "confirmation",
    "comparator",
)

DIRECTORIES = (
    "code",
    "data",
    "reports",
    "figures",
)

FILES = (
    "run_stages.sh",
    "build_release.py",
)

KEYS = ["group", "chemostat", "draw", "method"]

PROBABILITY_ATOL = 1e-7
PROBABILITY_RTOL = 1e-6


def compare_predictions(archived, regenerated):
    def load(path):
        frame = pd.read_csv(
            path,
            dtype={"chemostat": str},
        )

        if frame.duplicated(KEYS).any():
            raise ValueError("Duplicate prediction keys.")

        return frame.sort_values(KEYS).reset_index(drop=True)

    a = load(archived)
    b = load(regenerated)

    pd.testing.assert_frame_equal(
        a[KEYS + ["label"]],
        b[KEYS + ["label"]],
    )

    np.testing.assert_allclose(
        a["probability"].to_numpy(),
        b["probability"].to_numpy(),
        atol=PROBABILITY_ATOL,
        rtol=PROBABILITY_RTOL,
    )

    return {
        "prediction_rows": len(a),
        "max_absolute_probability_difference": float(
            np.max(np.abs(
                a["probability"].to_numpy()
                - b["probability"].to_numpy()
            ))
        ),
        "probability_atol": PROBABILITY_ATOL,
        "probability_rtol": PROBABILITY_RTOL,
    }


def reproduce_ecology(root):
    root = Path(root).resolve()
    before = verify_archives(root)

    with TemporaryDirectory(
        prefix="entropy-ecology-"
    ) as directory:

        work = Path(directory)

        for name in DIRECTORIES:
            shutil.copytree(
                root / name,
                work / name,
                ignore=shutil.ignore_patterns(
                    "__pycache__",
                    "*.pyc",
                ),
            )

        for name in FILES:
            shutil.copy2(root / name, work / name)

        for stage in STAGES:
            print(f"Running isolated {stage} stage...", flush=True)

            subprocess.run(
                ["sh", "run_stages.sh", stage],
                cwd=work,
                check=True,
            )

        subprocess.run(
            [sys.executable, "build_release.py"],
            cwd=work,
            check=True,
        )

        perturbation = pd.read_csv(
            work / "code/perturbation/summary.csv"
        ).sort_values("delta_N")

        metrics = pd.read_csv(
            work / "code/comparator/matched_metrics.csv"
        ).set_index("method")

        bootstrap = pd.read_csv(
            work / "code/comparator/cluster_bootstrap.csv"
        )

        verification = json.loads(
            (work / "reports/verification.json").read_text()
        )

        comparison = compare_predictions(
            root / "code/comparator/matched_predictions.csv",
            work / "code/comparator/matched_predictions.csv",
        )

        assert verification["all_checks_passed"]

    after = verify_archives(root)

    if before != after:
        raise RuntimeError("Archived evidence changed.")

    verification["isolated_regeneration"] = True
    verification["verified_archived_files"] = len(after)
    verification["prediction_comparison"] = comparison

    return perturbation, metrics, bootstrap, verification


if __name__ == "__main__":
    root = Path(__file__).resolve().parent

    _, _, _, verification = reproduce_ecology(root)

    print(json.dumps(verification, indent=2))
