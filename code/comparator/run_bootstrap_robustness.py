"""Retrospective ecological bootstrap sensitivity analysis."""

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

OUT = Path(__file__).resolve().parent

N_BOOT = 5000

SCHEMES = {
    "vessel_unstratified": 2026112601,
    "vessel_stratified": 2026112602,
    "draw_stratified": 2026112603,
    "crossed_stratified": 2026112604,
}

REFERENCE = "entropy_exact"

COMPARATORS = (
    "goal_linear",
    "goal_plus_curvature",
    "fisher_norm",
    "fisher_alignment",
    "state_norm",
    "goal_linear_plus_curvature",
)

METHODS = (REFERENCE,) + COMPARATORS
KEYS = ["group", "chemostat", "draw"]


def prepare_predictions(frame):
    required = set(KEYS + [
        "label", "method", "probability"
    ])

    if not required.issubset(frame.columns):
        raise ValueError("Missing required columns.")

    if frame[list(required)].isna().any().any():
        raise ValueError("Missing prediction values.")

    if frame.duplicated(KEYS + ["method"]).any():
        raise ValueError("Duplicate predictions.")

    if set(frame.method) != set(METHODS):
        raise ValueError("Unexpected diagnostic methods.")

    if not frame.label.isin([0, 1]).all():
        raise ValueError("Invalid warning labels.")

    if not frame.probability.between(0, 1).all():
        raise ValueError("Invalid probabilities.")

    if not frame.groupby(KEYS).label.nunique().eq(1).all():
        raise ValueError("Inconsistent warning labels.")

    wide = frame.pivot(
        index=KEYS + ["label"],
        columns="method",
        values="probability",
    ).reset_index()

    wide = wide.sort_values(KEYS).reset_index(drop=True)

    if wide[list(METHODS)].isna().any().any():
        raise ValueError("Incomplete diagnostic predictions.")

    if len(wide) != 144 or wide.group.nunique() != 2:
        raise ValueError("Unexpected confirmation design.")

    grids = {}

    for group, g in wide.groupby("group", sort=True):
        grid = (
            g.assign(row_index=g.index)
            .pivot(
                index="chemostat",
                columns="draw",
                values="row_index",
            )
            .sort_index()
            .sort_index(axis=1)
        )

        if grid.shape != (12, 6):
            raise ValueError("Unexpected crossed design.")

        if grid.isna().any().any():
            raise ValueError("Incomplete crossed design.")

        grids[group] = grid.to_numpy(dtype=int)

    return wide, grids


def sample_indices(grids, scheme, rng):
    if scheme not in SCHEMES:
        raise ValueError(f"Unknown scheme: {scheme}")

    if scheme == "vessel_unstratified":
        vessels = np.concatenate(
            list(grids.values()), axis=0
        )
        selected = rng.integers(
            0, len(vessels), size=len(vessels)
        )
        return vessels[selected].ravel()

    blocks = []

    for grid in grids.values():
        nv, nd = grid.shape

        vessels = np.arange(nv)
        draws = np.arange(nd)

        if scheme in (
            "vessel_stratified",
            "crossed_stratified",
        ):
            vessels = rng.integers(0, nv, size=nv)

        if scheme in (
            "draw_stratified",
            "crossed_stratified",
        ):
            draws = rng.integers(0, nd, size=nd)

        blocks.append(
            grid[np.ix_(vessels, draws)].ravel()
        )

    return np.concatenate(blocks)


def compute_contrasts(labels, probabilities):
    labels = np.asarray(labels)
    probabilities = np.asarray(probabilities)

    brier = np.mean(
        (probabilities - labels[:, None]) ** 2,
        axis=0,
    )

    if np.unique(labels).size == 2:
        auc = np.array([
            roc_auc_score(labels, probabilities[:, j])
            for j in range(len(METHODS))
        ])
    else:
        auc = np.full(len(METHODS), np.nan)

    result = {}

    for j, comparator in enumerate(COMPARATORS, start=1):
        result[f"auc_entropy_minus_{comparator}"] = (
            auc[0] - auc[j]
        )
        result[f"brier_entropy_minus_{comparator}"] = (
            brier[0] - brier[j]
        )

    return result


def run_bootstrap(wide, grids, n_boot=N_BOOT):
    labels = wide.label.to_numpy(dtype=int)
    probabilities = wide[list(METHODS)].to_numpy(float)

    rows = []

    for scheme, seed in SCHEMES.items():
        rng = np.random.default_rng(seed)

        print(
            f"Running {scheme}: {n_boot} replicates",
            flush=True,
        )

        for iteration in range(n_boot):
            idx = sample_indices(grids, scheme, rng)

            if len(idx) != len(wide):
                raise RuntimeError("Invalid resample size.")

            row = {
                "scheme": scheme,
                "iteration": iteration,
            }

            row.update(
                compute_contrasts(
                    labels[idx],
                    probabilities[idx],
                )
            )

            rows.append(row)

        print(f"Completed {scheme}", flush=True)

    return pd.DataFrame(rows)


def summarize(wide, boot):
    estimates = compute_contrasts(
        wide.label.to_numpy(dtype=int),
        wide[list(METHODS)].to_numpy(float),
    )

    rows = []

    for scheme in SCHEMES:
        subset = boot[boot.scheme == scheme]

        for comparator in COMPARATORS:
            for metric in ("auc", "brier"):
                column = (
                    f"{metric}_entropy_minus_{comparator}"
                )

                values = subset[column].dropna().to_numpy()

                if len(values):
                    lo, hi = np.percentile(
                        values, [2.5, 97.5]
                    )
                    mean = np.mean(values)
                else:
                    lo = hi = mean = np.nan

                sd = (
                    np.std(values, ddof=1)
                    if len(values) > 1
                    else np.nan
                )

                rows.append({
                    "scheme": scheme,
                    "comparator": comparator,
                    "metric": metric,
                    "estimate": estimates[column],
                    "bootstrap_mean": mean,
                    "bootstrap_sd": sd,
                    "ci_low": lo,
                    "ci_high": hi,
                    "ci_width": hi - lo,
                    "valid_replicates": len(values),
                    "undefined_replicates": (
                        len(subset) - len(values)
                    ),
                })

    return pd.DataFrame(rows)


def main():
    frame = pd.read_csv(
        OUT / "matched_predictions.csv",
        dtype={"chemostat": str},
    )

    wide, grids = prepare_predictions(frame)

    print("Validated configurations:", len(wide))
    print("Ecological groups:", list(grids))

    boot = run_bootstrap(wide, grids)
    summary = summarize(wide, boot)

    boot.to_csv(
        OUT / "bootstrap_robustness_replicates.csv",
        index=False,
    )

    summary.to_csv(
        OUT / "bootstrap_robustness_summary.csv",
        index=False,
    )

    print("\nPRIMARY COMPARISON: ENTROPY VS FISHER NORM")

    print(
        summary[
            summary.comparator == "fisher_norm"
        ].to_string(index=False)
    )

    print("\nALL COMPARISONS")

    print(summary.to_string(index=False))
    print("\nDONE")


if __name__ == "__main__":
    main()
