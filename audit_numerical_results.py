"""Independent, read-only audit of archived numerical results."""

from pathlib import Path
from itertools import product

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
COMP = ROOT / "code" / "compartment"
ECO = ROOT / "code" / "comparator"

ATOL = 1e-12

METHODS = {
    "entropy_exact",
    "goal_linear",
    "fisher_norm",
    "fisher_alignment",
    "state_norm",
}

ECO_METHODS = METHODS | {
    "goal_plus_curvature",
    "goal_linear_plus_curvature",
}

TARGETS = {
    "divergence": ("D_error1", 0.01),
    "trajectory": ("state_error1", 0.001),
}

SCHEMES = {
    "vessel_unstratified",
    "vessel_stratified",
    "draw_stratified",
    "crossed_stratified",
}


def read_csv(path, expected_rows):
    df = pd.read_csv(path)

    assert len(df) == expected_rows, path
    assert not df.isna().any().any(), path

    numeric = df.select_dtypes(include="number")
    assert np.isfinite(numeric.to_numpy()).all(), path

    return df


def close(actual, expected):
    np.testing.assert_allclose(
        actual,
        expected,
        atol=ATOL,
        rtol=0,
        equal_nan=False,
    )


def validate_predictions(labels, probabilities):
    y = np.asarray(labels)
    p = np.asarray(probabilities, dtype=float)

    assert y.ndim == p.ndim == 1
    assert len(y) == len(p) > 0
    assert np.isin(y, [0, 1]).all()
    assert np.isfinite(p).all()
    assert ((p >= 0) & (p <= 1)).all()

    return y, p


def independent_auc(labels, probabilities):
    """Mann-Whitney AUC with average ranks for ties."""

    y, p = validate_predictions(labels, probabilities)

    n_pos = int(np.sum(y == 1))
    n_neg = int(np.sum(y == 0))

    assert n_pos > 0 and n_neg > 0

    ranks = pd.Series(p).rank(method="average").to_numpy()

    rank_sum = ranks[y == 1].sum()

    return float(
        (rank_sum - n_pos * (n_pos + 1) / 2)
        / (n_pos * n_neg)
    )


def independent_brier(labels, probabilities):
    y, p = validate_predictions(labels, probabilities)
    return float(np.mean((p - y) ** 2))


def independent_log_loss(labels, probabilities):
    y, p = validate_predictions(labels, probabilities)

    eps = np.finfo(np.float64).eps
    p = np.clip(p, eps, 1 - eps)

    return float(
        -np.mean(y * np.log(p) + (1 - y) * np.log1p(-p))
    )


def percentile_interval(values):
    values = np.asarray(values, dtype=float)

    assert np.isfinite(values).all()
    assert len(values) > 0

    return np.quantile(values, [0.025, 0.975])


def audit_original_cases():
    cases = read_csv(COMP / "cases.csv", 336)

    assert not cases.duplicated(["split", "case_id"]).any()

    expected = (cases["D_error1"] > 0.01).astype(int)

    np.testing.assert_array_equal(cases["label"], expected)

    discovery = cases[cases["split"] == "discovery"]
    confirmation = cases[cases["split"] == "confirmation"]

    assert len(discovery) == 192
    assert len(confirmation) == 144

    assert int(discovery["label"].sum()) == 124

    print("PASS: original synthetic divergence labels")

    return discovery


def audit_prospective(discovery):
    cases = read_csv(
        COMP / "target_confirmation_cases.csv", 1000
    )
    counts = read_csv(
        COMP / "target_confirmation_counts.csv", 2
    )
    predictions = read_csv(
        COMP / "target_confirmation_predictions.csv", 10000
    )
    metrics = read_csv(
        COMP / "target_confirmation_metrics.csv", 10
    )
    boot = read_csv(
        COMP / "target_confirmation_bootstrap.csv", 5000
    )
    summary = read_csv(
        COMP / "target_confirmation_summary.csv", 5
    )

    assert set(cases["case_id"]) == set(range(1000))
    assert cases["case_id"].is_unique
    assert set(cases["split"]) == {"prospective_confirmation"}

    assert not predictions.duplicated(
        ["target", "case_id", "method"]
    ).any()

    assert not metrics.duplicated(
        ["target", "method"]
    ).any()

    assert counts["target"].is_unique
    assert summary["comparison"].is_unique

    assert set(predictions["target"]) == set(TARGETS)
    assert set(counts["target"]) == set(TARGETS)
    assert set(predictions["method"]) == METHODS

    assert set(boot["iteration"]) == set(range(5000))

    cases = cases.set_index("case_id")
    counts = counts.set_index("target")
    metrics = metrics.set_index(["target", "method"])
    summary = summary.set_index("comparison")

    close(cases["epsilon"], 0.2 * cases["x10"])

    assert (cases["max_identity_error"] >= 0).all()
    assert (cases["max_identity_error"] <= 1e-9).all()

    np.testing.assert_array_equal(
        cases["label"],
        (cases["D_error1"] > 0.01).astype(int),
    )

    computed = {}

    for target, (column, threshold) in TARGETS.items():
        expected = (cases[column] > threshold).astype(int)

        count = counts.loc[target]

        assert int(expected.sum()) == count["confirmation_positive"]
        assert int((expected == 0).sum()) == count["confirmation_negative"]

        close(expected.mean(), count["confirmation_prevalence"])

        assert (
            count["discovery_positive"]
            + count["discovery_negative"]
        ) == 192

        if target == "divergence":
            assert discovery["label"].sum() == count["discovery_positive"]

        for method in sorted(METHODS):
            block = predictions[
                (predictions["target"] == target)
                & (predictions["method"] == method)
            ].sort_values("case_id")

            assert len(block) == 1000

            np.testing.assert_array_equal(
                block["case_id"], np.arange(1000)
            )

            np.testing.assert_array_equal(
                block["label"],
                expected.reindex(block["case_id"]).to_numpy(),
            )

            auc = independent_auc(
                block["label"], block["probability"]
            )
            brier = independent_brier(
                block["label"], block["probability"]
            )

            close(auc, metrics.loc[(target, method), "auc"])
            close(brier, metrics.loc[(target, method), "brier"])

            computed[target, method] = {
                "auc": auc,
                "brier": brier,
            }

    estimates = {}

    for target in TARGETS:
        for metric in ("auc", "brier"):
            estimates[f"{metric}_delta_{target}"] = (
                computed[target, "entropy_exact"][metric]
                - computed[target, "fisher_norm"][metric]
            )

    estimates["auc_reversal_contrast"] = (
        estimates["auc_delta_divergence"]
        - estimates["auc_delta_trajectory"]
    )

    close(
        boot["auc_reversal_contrast"],
        boot["auc_delta_divergence"]
        - boot["auc_delta_trajectory"],
    )

    assert set(summary.index) == set(estimates)

    for name, estimate in estimates.items():
        row = summary.loc[name]
        low, high = percentile_interval(boot[name])

        close(estimate, row["estimate"])
        close(low, row["ci_low"])
        close(high, row["ci_high"])

        assert row["valid_replicates"] == 5000
        assert row["undefined_replicates"] == 0

        print(
            f"PASS: {name}: {estimate:.6f} "
            f"[{low:.6f}, {high:.6f}]"
        )

    assert summary.loc["auc_delta_divergence", "ci_low"] > 0
    assert summary.loc["auc_delta_trajectory", "ci_high"] < 0

    print("PASS: all 10 prospective metric records")
    print("PASS: all 5 prospective bootstrap summaries")
    print("PASS: all 5,000 replicate-level reversal identities")
    print("PASS: prespecified ranking reversal")

    print(
        "LIMITATION: original discovery trajectory labels "
        "are not independently verified from cases.csv."
    )


def audit_ecology():
    predictions = read_csv(
        ECO / "matched_predictions.csv", 1008
    )
    metrics = read_csv(
        ECO / "matched_metrics.csv", 7
    )
    boot = read_csv(
        ECO / "bootstrap_robustness_replicates.csv", 20000
    )
    summary = read_csv(
        ECO / "bootstrap_robustness_summary.csv", 48
    )

    keys = ["group", "chemostat", "draw"]

    assert not predictions.duplicated(keys + ["method"]).any()

    assert set(predictions["method"]) == ECO_METHODS
    assert set(metrics["method"]) == ECO_METHODS
    assert metrics["method"].is_unique

    grouped = predictions.groupby(keys)

    assert grouped.ngroups == 144
    assert (grouped["method"].nunique() == 7).all()
    assert (grouped["label"].nunique() == 1).all()

    configurations = predictions[keys].drop_duplicates()

    assert set(configurations["group"]) == {
        "monocultures", "polycultures"
    }

    for group, block in configurations.groupby("group"):
        assert len(block) == 72
        assert block["chemostat"].nunique() == 12
        assert block["draw"].nunique() == 6

    metrics = metrics.set_index("method")
    computed = {}

    for method in sorted(ECO_METHODS):
        block = predictions[predictions["method"] == method]

        assert len(block) == 144

        y = block["label"]
        p = block["probability"]

        values = {
            "auc": independent_auc(y, p),
            "brier": independent_brier(y, p),
            "log_loss": independent_log_loss(y, p),
        }

        for metric, value in values.items():
            close(value, metrics.loc[method, metric])

        computed[method] = values

    print("PASS: all 7 ecological metric records")

    assert set(boot["scheme"]) == SCHEMES
    assert not boot.duplicated(["scheme", "iteration"]).any()

    comparators = ECO_METHODS - {"entropy_exact"}

    expected_keys = set(
        product(SCHEMES, comparators, ["auc", "brier"])
    )

    assert not summary.duplicated(
        ["scheme", "comparator", "metric"]
    ).any()

    summary = summary.set_index(
        ["scheme", "comparator", "metric"]
    )

    assert set(summary.index) == expected_keys

    for scheme in sorted(SCHEMES):
        block = boot[boot["scheme"] == scheme]

        assert len(block) == 5000
        assert set(block["iteration"]) == set(range(5000))

        for comparator in sorted(comparators):
            for metric in ("auc", "brier"):
                column = f"{metric}_entropy_minus_{comparator}"
                values = block[column].to_numpy()

                row = summary.loc[(scheme, comparator, metric)]

                estimate = (
                    computed["entropy_exact"][metric]
                    - computed[comparator][metric]
                )

                low, high = percentile_interval(values)

                close(estimate, row["estimate"])
                close(values.mean(), row["bootstrap_mean"])
                close(values.std(ddof=1), row["bootstrap_sd"])
                close(low, row["ci_low"])
                close(high, row["ci_high"])
                close(high - low, row["ci_width"])

                assert row["valid_replicates"] == 5000
                assert row["undefined_replicates"] == 0

                if comparator == "goal_plus_curvature":
                    close(values, 0)

    print("PASS: all 48 ecological bootstrap summaries")

    for metric in ("auc", "brier"):
        baseline = summary.loc[
            ("vessel_unstratified", "fisher_norm", metric)
        ]

        crossed = summary.loc[
            ("crossed_stratified", "fisher_norm", metric)
        ]

        ratio = crossed["ci_width"] / baseline["ci_width"]

        assert ratio > 2

        print(
            f"PASS: {metric.upper()} crossed interval "
            f"width ratio = {ratio:.6f}"
        )

        for scheme in SCHEMES:
            row = summary.loc[(scheme, "fisher_norm", metric)]
            assert row["ci_low"] <= 0 <= row["ci_high"]

    print("PASS: ecological near-tie remains unresolved")


def main():
    print("INDEPENDENT NUMERICAL CONSISTENCY AUDIT\n")

    discovery = audit_original_cases()
    audit_prospective(discovery)
    audit_ecology()

    print("\nALL NUMERICAL AUDIT CHECKS PASSED")


if __name__ == "__main__":
    main()
