"""Prospective confirmation of warning-target-dependent rankings.

Protocol: protocols/target_reversal_confirmation.md
Protocol commit: ef3352f

Run from repository root:
    python code/compartment/run_target_confirmation.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score

from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from run_benchmark import FEATURES, design, generate, one_case

OUT = Path(__file__).resolve().parent

N_DISCOVERY = 192
N_CONFIRMATION = 1000
CONFIRMATION_SEED = 2026112301
BOOT_SEED = 2026112302
N_BOOT = 5000

TARGETS = {
    "divergence": ("D_error1", 0.01),
    "trajectory": ("state_error1", 0.001),
}


def generate_fresh_cases():
    rng = np.random.default_rng(CONFIRMATION_SEED)
    rows = []

    for case_id in range(N_CONFIRMATION):
        try:
            rows.append(
                one_case("prospective_confirmation", case_id, rng)
            )
        except Exception as exc:
            pd.DataFrame([{
                "case_id": case_id,
                "seed": CONFIRMATION_SEED,
                "error": repr(exc),
            }]).to_csv(
                OUT / "target_confirmation_failure.csv",
                index=False,
            )
            raise RuntimeError(
                f"Confirmation stopped at case {case_id}; "
                "no replacement cases generated."
            ) from exc

    return pd.DataFrame(rows)


def apply_target(frame, target):
    column, threshold = TARGETS[target]
    result = frame.copy()
    result["label"] = (result[column] > threshold).astype(int)
    return result


def safe_auc(y, probability):
    if np.unique(y).size < 2:
        return np.nan
    return float(roc_auc_score(y, probability))


def evaluate_target(discovery, confirmation, target):
    train = apply_target(discovery, target)
    test = apply_target(confirmation, target)

    counts = {
        "target": target,
        "discovery_negative": int((train.label == 0).sum()),
        "discovery_positive": int((train.label == 1).sum()),
        "confirmation_negative": int((test.label == 0).sum()),
        "confirmation_positive": int((test.label == 1).sum()),
        "confirmation_prevalence": float(test.label.mean()),
    }

    print(counts, flush=True)

    if train.label.nunique() < 2:
        raise RuntimeError(
            f"{target}: discovery has only one warning class."
        )

    metric_rows = []
    prediction_rows = []
    y = test.label.to_numpy()

    for method, columns in FEATURES.items():
        clf = make_pipeline(
            StandardScaler(),
            LogisticRegression(C=1., max_iter=2000),
        )
        clf.fit(design(train, columns), train.label)
        probs = clf.predict_proba(design(test, columns))[:, 1]

        metric_rows.append({
            "method": method,
            "auc": safe_auc(y, probs),
            "brier": float(brier_score_loss(y, probs)),
        })

        prediction_rows.extend(
            {
                "case_id": int(case_id),
                "method": method,
                "label": int(label),
                "probability": float(probability),
            }
            for case_id, label, probability in zip(
                test.case_id, y, probs
            )
        )

    metrics = pd.DataFrame(metric_rows)
    predictions = pd.DataFrame(prediction_rows)
    metrics.insert(0, "target", target)
    predictions.insert(0, "target", target)

    return metrics, predictions, counts


def paired_bootstrap(predictions, n_boot=N_BOOT):
    frames = {}

    for target in TARGETS:
        subset = predictions[predictions.target == target]

        frames[target] = subset.pivot(
            index=["case_id", "label"],
            columns="method",
            values="probability",
        ).reset_index().sort_values("case_id").reset_index(drop=True)

    divergence = frames["divergence"]
    trajectory = frames["trajectory"]

    if not np.array_equal(
        divergence.case_id.to_numpy(),
        trajectory.case_id.to_numpy(),
    ):
        raise ValueError("Confirmation case IDs do not match.")

    if not divergence.case_id.is_unique:
        raise ValueError("Duplicate confirmation case IDs.")

    arrays = {}

    for target, frame in frames.items():
        arrays[target] = (
            frame.label.to_numpy(),
            frame.entropy_exact.to_numpy(),
            frame.fisher_norm.to_numpy(),
        )

    rng = np.random.default_rng(BOOT_SEED)
    rows = []
    n = len(divergence)

    for iteration in range(n_boot):
        idx = rng.integers(0, n, size=n)
        row = {"iteration": iteration}

        for target, (y, entropy, fisher) in arrays.items():
            yy = y[idx]
            ee = entropy[idx]
            ff = fisher[idx]

            row[f"auc_delta_{target}"] = (
                safe_auc(yy, ee) - safe_auc(yy, ff)
            )

            row[f"brier_delta_{target}"] = (
                brier_score_loss(yy, ee)
                - brier_score_loss(yy, ff)
            )

        row["auc_reversal_contrast"] = (
            row["auc_delta_divergence"]
            - row["auc_delta_trajectory"]
        )

        rows.append(row)

    return pd.DataFrame(rows)


def summarize(metrics, boot):
    indexed = metrics.set_index(["target", "method"])
    estimates = {}

    for target in TARGETS:
        for metric in ("auc", "brier"):
            estimates[f"{metric}_delta_{target}"] = float(
                indexed.loc[(target, "entropy_exact"), metric]
                - indexed.loc[(target, "fisher_norm"), metric]
            )

    estimates["auc_reversal_contrast"] = (
        estimates["auc_delta_divergence"]
        - estimates["auc_delta_trajectory"]
    )

    rows = []

    for column, estimate in estimates.items():
        values = boot[column].dropna()

        lo, hi = (
            np.percentile(values, [2.5, 97.5])
            if len(values)
            else (np.nan, np.nan)
        )

        rows.append({
            "comparison": column,
            "estimate": estimate,
            "ci_low": lo,
            "ci_high": hi,
            "valid_replicates": len(values),
            "undefined_replicates": len(boot) - len(values),
        })

    return pd.DataFrame(rows)


def main():
    discovery = generate("discovery", N_DISCOVERY)
    confirmation = generate_fresh_cases()

    confirmation.to_csv(
        OUT / "target_confirmation_cases.csv", index=False
    )

    all_metrics = []
    all_predictions = []
    all_counts = []

    for target in TARGETS:
        metrics, predictions, counts = evaluate_target(
            discovery, confirmation, target
        )

        all_metrics.append(metrics)
        all_predictions.append(predictions)
        all_counts.append(counts)

    metrics = pd.concat(all_metrics, ignore_index=True)
    predictions = pd.concat(all_predictions, ignore_index=True)
    counts = pd.DataFrame(all_counts)

    boot = paired_bootstrap(predictions)
    summary = summarize(metrics, boot)

    for name, frame in (
        ("metrics", metrics),
        ("predictions", predictions),
        ("counts", counts),
        ("bootstrap", boot),
        ("summary", summary),
    ):
        frame.to_csv(
            OUT / f"target_confirmation_{name}.csv",
            index=False,
        )

    print("\nCONFIRMATION METRICS")
    print(metrics.to_string(index=False))

    print("\nPAIRED BOOTSTRAP RESULTS")
    print(summary.to_string(index=False))

    intervals = summary.set_index("comparison")

    supported = (
        intervals.loc["auc_delta_divergence", "ci_low"] > 0
        and intervals.loc["auc_delta_trajectory", "ci_high"] < 0
    )

    print(f"\nPrespecified reversal supported: {supported}")
    print("DONE")


if __name__ == "__main__":
    main()
