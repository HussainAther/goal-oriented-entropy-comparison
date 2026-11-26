"""Independent consistency audit of archived geometry evidence.

Does not import the geometry experiment or integrate the ODE.
"""

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


ROOT = Path(__file__).resolve().parent

PAIRS = {
    "B_S": "error_S",
    "B_C": "error_C",
    "B_D": "error_D",
    "B_quad": "error_quad",
    "B_r2": "error_alignment_quad",
    "B_alignment": "error_alignment_exact",
}

RECONSTRUCTION_FIELDS = (
    "entropy_exact",
    "goal_linear",
    "fisher_norm",
    "fisher_alignment",
    "state_norm",
    "D_error1",
    "D_error2",
    "state_error1",
)


def check(actual, expected, name, atol=1e-12, rtol=1e-9):
    np.testing.assert_allclose(
        actual,
        expected,
        atol=atol,
        rtol=rtol,
        equal_nan=False,
        err_msg=name,
    )


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def violations(error, bound):
    return error > bound + 1e-10 + 1e-8 * np.abs(bound)


def audit_geometry(root=ROOT):
    directory = Path(root) / "code/compartment"

    p = pd.read_csv(
        directory / "geometry_validation_pointwise.csv"
    )
    c = pd.read_csv(
        directory / "geometry_validation_cases.csv"
    )
    s = pd.read_csv(
        directory / "geometry_validation_summary.csv"
    )
    original = pd.read_csv(
        directory / "target_confirmation_cases.csv"
    )

    require(p.shape == (40000, 58), "Pointwise shape")
    require(c.shape == (1000, 30), "Configuration shape")
    require(s.shape == (268, 4), "Summary shape")

    for name, frame in (
        ("pointwise", p),
        ("configuration", c),
        ("summary", s),
        ("original", original),
    ):
        require(
            not frame.isna().any().any(),
            f"Missing values: {name}",
        )
        require(
            np.isfinite(
                frame.select_dtypes(include="number")
            ).all().all(),
            f"Nonfinite values: {name}",
        )

    require(
        not p.duplicated(["case_id", "time"]).any(),
        "Duplicate pointwise records",
    )
    require(
        not c.case_id.duplicated().any(),
        "Duplicate configuration identifiers",
    )
    require(
        not original.case_id.duplicated().any(),
        "Duplicate original identifiers",
    )

    p = p.sort_values(
        ["case_id", "time"]
    ).reset_index(drop=True)

    c = c.set_index("case_id").sort_index()
    original = original.set_index("case_id").sort_index()

    ids = np.arange(1000)

    np.testing.assert_array_equal(c.index, ids)
    np.testing.assert_array_equal(original.index, ids)
    np.testing.assert_array_equal(
        p.case_id, np.repeat(ids, 40)
    )

    check(
        p.time,
        np.tile(np.arange(1, 41) / 8, 1000),
        "Positive-time grid",
        atol=0,
        rtol=0,
    )

    print("PASS: identifiers, time grid and finite records")

    # Labels reconstructed from archived warning targets.
    labels = {
        "divergence_label":
            (original.D_error1 > 0.01).astype(int),
        "trajectory_label":
            (original.state_error1 > 0.001).astype(int),
    }

    np.testing.assert_array_equal(
        original.label, labels["divergence_label"]
    )

    for name, expected in labels.items():
        np.testing.assert_array_equal(c[name], expected)
        np.testing.assert_array_equal(
            p[name], expected.loc[p.case_id]
        )

    print("PASS: both warning targets and archived labels")

    # Independent algebraic checks using archived primitives.
    Q, H, L = p.Q, p.H, p.L
    d1, S, C = p.d1, p.S, p.C

    require((Q > 0).all(), "Nonpositive Q")
    require((H > 0).all(), "Nonpositive H")
    require((d1 > 0).all(), "Nonpositive divergence")

    r = np.sqrt(H / Q)
    cosine = L / np.sqrt(Q * H)
    alignment = 2 * np.abs(L) / Q
    exact = np.abs(S + C) / d1
    quadratic = np.abs(2 * L + H) / Q

    require(
        (np.abs(cosine) <= 1 + 1e-12).all(),
        "Invalid alignment cosine",
    )

    identities = {
        "r": r,
        "c": cosine,
        "twice_alignment": alignment,
        "E_loc": exact,
        "E_quad": quadratic,
        "entropy_regularized":
            np.abs(S + C) / np.maximum(d1, 1e-12),
        "goal_linear_regularized":
            np.abs(S) / np.maximum(d1, 1e-12),
        "alignment_regularized":
            np.abs(L) / np.maximum(Q, 1e-12),
        "fisher_norm_regularized":
            np.sqrt(H) / np.maximum(np.sqrt(Q), 1e-12),
    }

    for name, expected in identities.items():
        check(p[name], expected, name)

    conditions = pd.DataFrame({
        "condition_a": p.a < 1,
        "condition_b": p.b < 1,
        "condition_Q": Q > 0,
        "condition_entropy_floor": d1 >= 1e-12,
        "condition_alignment_floor": Q >= 1e-12,
        "condition_remainder": 2 * p.B_D < Q,
    })

    conditions["applicable"] = conditions.all(axis=1)

    for name in conditions:
        require(
            pd.api.types.is_bool_dtype(p[name]),
            f"Invalid condition type: {name}",
        )
        np.testing.assert_array_equal(
            p[name], conditions[name]
        )

    require(
        conditions.applicable.all(),
        "Local sufficient conditions not satisfied",
    )

    print("PASS: diagnostic identities and sufficient conditions")

    # Elementary bounds are archived inputs.
    # Their component-level construction is not re-audited.
    for name in ("B_S", "B_C", "B_D"):
        require((p[name] > 0).all(), name)

    combined = (
        2 * (p.B_S + p.B_C + quadratic * p.B_D)
        / (Q - 2 * p.B_D)
    )

    check(p.B_quad, combined, "Combined quadratic bound")
    check(p.B_r2, r**2, "Quadratic alignment bound")
    check(
        p.B_alignment,
        combined + r**2,
        "Combined alignment bound",
    )

    errors = {
        "error_S": np.abs(S - L),
        "error_C": np.abs(C - H / 2),
        "error_D": np.abs(d1 - Q / 2),
        "error_quad": np.abs(exact - quadratic),
        "error_alignment_quad":
            np.abs(quadratic - alignment),
        "error_alignment_exact":
            np.abs(exact - alignment),
    }

    for name, expected in errors.items():
        check(p[name], expected, name)

    for bound, error in PAIRS.items():
        require((p[bound] > 0).all(), bound)

        # Use archived errors for ratio verification to avoid
        # amplifying CSV rounding in very small remainders.
        check(
            p[f"ratio_{bound}"],
            p[error] / p[bound],
            f"Ratio: {bound}",
        )

        check(
            p[f"excess_{bound}"],
            p[error] - p[bound],
            f"Excess: {bound}",
        )

        expected = violations(errors[error], p[bound])

        np.testing.assert_array_equal(
            p[f"violation_{bound}"], expected
        )

        require(
            not expected.any(),
            f"Detected inequality violation: {bound}",
        )

    print("PASS: 240,000 pointwise inequality comparisons")

    # Independently aggregate the pointwise records.
    grouped = p.groupby("case_id", sort=True)

    for name in conditions:
        check(
            c[f"fraction_{name}"],
            grouped[name].mean(),
            f"Configuration fraction: {name}",
        )

    for name in (
        "E_loc", "E_quad", "twice_alignment", "r", "c"
    ):
        check(
            c[f"mean_{name}"],
            grouped[name].mean(),
            f"Configuration mean: {name}",
        )

    check(
        c.mean_abs_c,
        p.c.abs().groupby(p.case_id).mean(),
        "Mean absolute cosine",
    )

    np.testing.assert_array_equal(
        c.defined_c_count, grouped.c.count()
    )

    np.testing.assert_array_equal(
        c.all_times_applicable,
        grouped.applicable.all(),
    )

    temporal_error = np.abs(
        grouped.entropy_regularized.mean()
        - 2 * grouped.alignment_regularized.mean()
    )

    temporal_bound = grouped.B_alignment.mean()

    check(c.temporal_error, temporal_error, "Temporal error")
    check(c.temporal_bound, temporal_bound, "Temporal bound")
    check(
        c.temporal_ratio,
        temporal_error / temporal_bound,
        "Temporal ratio",
    )

    temporal_flags = violations(
        temporal_error, temporal_bound
    )

    np.testing.assert_array_equal(
        c.temporal_violation, temporal_flags
    )

    require(
        not temporal_flags.any(),
        "Detected temporal inequality violation",
    )

    print("PASS: configuration aggregation and temporal bounds")

    # Four diagnostic features can be reconstructed directly
    # from the pointwise archive.
    diagnostic_columns = {
        "entropy_exact": "entropy_regularized",
        "goal_linear": "goal_linear_regularized",
        "fisher_norm": "fisher_norm_regularized",
        "fisher_alignment": "alignment_regularized",
    }

    for feature, column in diagnostic_columns.items():
        check(
            grouped[column].mean(),
            original[feature],
            f"Archived diagnostic: {feature}",
            atol=1e-10,
            rtol=1e-7,
        )

    for feature in RECONSTRUCTION_FIELDS:
        discrepancy = c[f"reconstruction_error_{feature}"]

        require(
            (discrepancy >= 0).all(),
            f"Negative reconstruction error: {feature}",
        )

        require(
            (
                discrepancy
                <= 1e-10 + 1e-7 * original[feature].abs()
            ).all(),
            f"Reconstruction tolerance: {feature}",
        )

    print("PASS: archived diagnostics and reconstruction tolerances")

    # Independently recompute the complete summary table.
    expected = {}

    def add(scope, metric, statistic, value):
        key = (scope, metric, statistic)
        require(key not in expected, f"Duplicate key: {key}")
        expected[key] = float(value)

    def describe(scope, metric, values):
        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]

        add(scope, metric, "count", len(values))

        for statistic, quantile in (
            ("p05", 0.05),
            ("median", 0.5),
            ("p95", 0.95),
            ("max", 1.0),
        ):
            add(
                scope, metric, statistic,
                np.quantile(values, quantile),
            )

    for name in conditions:
        add("pointwise", name, "count", p[name].sum())
        add("pointwise", name, "fraction", p[name].mean())

        describe(
            "configuration",
            f"fraction_{name}",
            c[f"fraction_{name}"],
        )

    add(
        "configuration",
        "all_times_applicable",
        "count",
        c.all_times_applicable.sum(),
    )

    for bound, error in PAIRS.items():
        describe("pointwise", error, p[error])

        describe(
            "pointwise",
            f"ratio_{bound}",
            p[f"ratio_{bound}"],
        )

        add(
            "pointwise", bound, "violations",
            p[f"violation_{bound}"].sum(),
        )

        add(
            "pointwise", bound, "max_positive_excess",
            max(0.0, p[f"excess_{bound}"].max()),
        )

    for name in (
        "temporal_error",
        "temporal_bound",
        "temporal_ratio",
    ):
        describe("configuration", name, c[name])

    add(
        "configuration", "temporal_bound", "violations",
        c.temporal_violation.sum(),
    )

    correlation_fields = (
        "mean_E_loc",
        "mean_E_quad",
        "mean_twice_alignment",
        "mean_r",
    )

    for left, right in combinations(correlation_fields, 2):
        x = c[left].to_numpy()
        y = c[right].to_numpy()

        metric = f"{left}_vs_{right}"

        add("correlation", metric, "count", len(x))

        add(
            "correlation", metric, "pearson",
            np.corrcoef(x, y)[0, 1],
        )

        add(
            "correlation", metric, "spearman",
            np.corrcoef(rankdata(x), rankdata(y))[0, 1],
        )

    for target in ("divergence_label", "trajectory_label"):
        for label in (0, 1):
            subset = c.loc[c[target] == label]
            scope = f"{target}_{label}"

            add(
                scope, "configurations", "count", len(subset)
            )

            for name in (
                "mean_r",
                "mean_c",
                "mean_abs_c",
                "mean_E_loc",
                "mean_twice_alignment",
            ):
                describe(scope, name, subset[name])

    for feature in RECONSTRUCTION_FIELDS:
        name = f"reconstruction_error_{feature}"

        add("reconstruction", name, "max", c[name].max())

    keys = ["scope", "metric", "statistic"]

    require(
        not s.duplicated(keys).any(),
        "Duplicate summary keys",
    )

    observed = s.set_index(keys).value.to_dict()

    require(
        len(expected) == 268,
        f"Unexpected recomputed summary count: {len(expected)}",
    )

    require(
        expected.keys() == observed.keys(),
        "Summary key mismatch",
    )

    for key, value in expected.items():
        check(observed[key], value, str(key))

    print("PASS: all 268 geometry summary records")
    print("PASS: all 12 correlation coefficients")
    print("PASS: all 20 warning-class medians")

    print(
        "LIMITATION: no independent ODE integration or "
        "component-level reconstruction of elementary bounds."
    )

    print("ALL GEOMETRY AUDIT CHECKS PASSED")

    return {
        "pointwise_records": len(p),
        "configurations": len(c),
        "summary_records": len(expected),
        "pointwise_inequalities": len(p) * len(PAIRS),
        "temporal_inequalities": len(c),
    }


if __name__ == "__main__":
    audit_geometry()
