"""Exploratory sensitivity analysis using trajectory-error warnings.

Run from repository root:
    python code/compartment/run_state_sensitivity.py

The 0.1% threshold was selected after inspecting error distributions,
including confirmation cases. Results are exploratory, not confirmatory.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from run_benchmark import (
    N_CONFIRMATION,
    N_DISCOVERY,
    evaluate,
    generate,
)

OUT = Path(__file__).resolve().parent
THRESHOLD = 0.001
N_BOOT = 1500
BOOT_SEED = 2026112104


def prepare_cases():
    discovery = generate("discovery", N_DISCOVERY)
    confirmation = generate("confirmation", N_CONFIRMATION)

    for frame in (discovery, confirmation):
        frame["divergence_label"] = frame["label"].copy()
        frame["label"] = (frame["state_error1"] > THRESHOLD).astype(int)

    return discovery, confirmation


def paired_bootstrap(predictions):
    wide = predictions.pivot(
        index=["case_id", "label"],
        columns="method",
        values="probability",
    ).reset_index()

    from sklearn.metrics import brier_score_loss, roc_auc_score

    rng = np.random.default_rng(BOOT_SEED)
    rows = []

    for iteration in range(N_BOOT):
        sample = wide.iloc[
            rng.integers(0, len(wide), size=len(wide))
        ]

        y = sample["label"].to_numpy()
        row = {"iteration": iteration}

        for comparator in ("fisher_norm", "fisher_alignment", "state_norm"):
            entropy = sample["entropy_exact"].to_numpy()
            other = sample[comparator].to_numpy()

            row[f"brier_entropy_minus_{comparator}"] = (
                brier_score_loss(y, entropy)
                - brier_score_loss(y, other)
            )

            row[f"auc_entropy_minus_{comparator}"] = (
                roc_auc_score(y, entropy)
                - roc_auc_score(y, other)
                if len(np.unique(y)) == 2
                else np.nan
            )

        rows.append(row)

    return pd.DataFrame(rows)


def main():
    discovery, confirmation = prepare_cases()

    for name, frame in (
        ("discovery", discovery),
        ("confirmation", confirmation),
    ):
        counts = frame["label"].value_counts().reindex(
            [0, 1], fill_value=0
        )
        print(
            f"{name}: negative={counts[0]}, positive={counts[1]}"
        )
        if frame["label"].nunique() < 2:
            print("STOP: warning class absent; metrics undefined")
            return 2

    metrics, predictions = evaluate(discovery, confirmation)
    boot = paired_bootstrap(predictions)

    pd.concat([discovery, confirmation], ignore_index=True).to_csv(
        OUT / "state_sensitivity_cases.csv", index=False
    )
    metrics.to_csv(OUT / "state_sensitivity_metrics.csv", index=False)
    predictions.to_csv(
        OUT / "state_sensitivity_predictions.csv", index=False
    )
    boot.to_csv(
        OUT / "state_sensitivity_bootstrap.csv", index=False
    )

    print("\nTRAJECTORY-ERROR WARNING RESULTS")
    print(metrics.to_string(index=False))

    for comparator in ("fisher_norm", "fisher_alignment", "state_norm"):
        for metric in ("auc", "brier"):
            column = f"{metric}_entropy_minus_{comparator}"
            lo, hi = np.nanpercentile(boot[column], [2.5, 97.5])
            print(
                f"{column}: 95% bootstrap CI "
                f"[{lo:+.6f}, {hi:+.6f}]"
            )

    print("\nDONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
