"""Retrospective assessment of local information geometry.

Protocol: protocols/geometry_validation.md
Protocol commit: e2be115

Run from repository root:
    python code/compartment/run_geometry_validation.py
"""

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from model import rhs, augmented_rhs
from run_benchmark import TIMES, integrate, divergence


OUT = Path(__file__).resolve().parent
N_CASES = 1000
ATOL = 1e-10
RTOL = 1e-7
FLOOR = 1e-12

CONDITIONS = [
    "condition_a",
    "condition_b",
    "condition_Q",
    "condition_entropy_floor",
    "condition_alignment_floor",
    "condition_remainder",
    "applicable",
]

GEOMETRY_COLUMNS = [
    "E_loc",
    "E_quad",
    "twice_alignment",
    "r",
    "c",
]

BOUND_PAIRS = [
    ("error_quad", "B_quad"),
    ("error_alignment_quad", "B_r2"),
    ("error_alignment_exact", "B_alignment"),
]


def safe_divide(numerator, denominator):
    """Return NaN when the denominator is not positive."""
    numerator, denominator = np.broadcast_arrays(
        np.asarray(numerator, dtype=float),
        np.asarray(denominator, dtype=float),
    )
    result = np.full(numerator.shape, np.nan)
    np.divide(
        numerator,
        denominator,
        out=result,
        where=denominator > 0,
    )
    return result


def bound_violation(error, bound):
    """Evaluate a bound with the protocol's numerical tolerance."""
    error = np.asarray(error, dtype=float)
    bound = np.asarray(bound, dtype=float)

    valid = np.isfinite(error) & np.isfinite(bound)
    tolerance = 1e-10 + 1e-8 * np.abs(bound)

    return valid & (error > bound + tolerance)


def pointwise_geometry(X, p1, h, omega):
    """Evaluate the manuscript's local geometry."""
    X = np.asarray(X, dtype=float)
    p1 = np.asarray(p1, dtype=float)
    h = np.asarray(h, dtype=float)
    omega = np.asarray(omega, dtype=float)

    p2 = p1 + h

    for array in (X, p1, h, p2, omega):
        if not np.isfinite(array).all():
            raise ValueError("Nonfinite geometry input")

    for array in (X, p1, p2, omega):
        if not np.all(array > 0):
            raise ValueError("Nonpositive geometry input")

    q = p1 - X
    G = 1.0 / (omega * p1)

    Q = np.sum(q**2 * G, axis=1)
    H = np.sum(h**2 * G, axis=1)
    L = np.sum(q * h * G, axis=1)

    r = np.sqrt(safe_divide(H, Q))
    c = safe_divide(L, np.sqrt(Q * H))
    alignment = safe_divide(np.abs(L), Q)

    a = np.max(np.abs(q) / p1, axis=1)
    b = np.max(np.abs(h) / p1, axis=1)

    d1 = divergence(p1, X, omega)
    S = np.sum(h * np.log(p1 / X) / omega, axis=1)
    C = divergence(p2, p1, omega)

    E_loc = safe_divide(np.abs(S + C), d1)
    E_quad = safe_divide(np.abs(2 * L + H), Q)

    condition_a = a < 1
    condition_b = b < 1
    condition_Q = Q > 0

    condition_entropy_floor = d1 >= FLOOR
    condition_alignment_floor = Q >= FLOOR

    B_S = np.full(len(Q), np.nan)
    B_C = np.full(len(Q), np.nan)
    B_D = np.full(len(Q), np.nan)

    m = condition_a
    B_S[m] = (
        np.sum(
            np.abs(h[m]) * q[m]**2 /
            (omega * p1[m]**2),
            axis=1,
        ) / (2 * (1 - a[m])**2)
    )

    B_D[m] = (
        np.sum(
            np.abs(q[m])**3 /
            (omega * p1[m]**2),
            axis=1,
        ) / (3 * (1 - a[m])**3)
    )

    m = condition_b
    B_C[m] = (
        np.sum(
            np.abs(h[m])**3 /
            (omega * p1[m]**2),
            axis=1,
        ) / (6 * (1 - b[m])**2)
    )

    condition_remainder = (
        np.isfinite(B_D) & (2 * B_D < Q)
    )

    applicable = (
        condition_a
        & condition_b
        & condition_Q
        & condition_entropy_floor
        & condition_alignment_floor
        & condition_remainder
    )

    B_quad = np.full(len(Q), np.nan)

    m = applicable

    B_quad[m] = (
        2 * (B_S[m] + B_C[m])
        + 2 * E_quad[m] * B_D[m]
    ) / (Q[m] - 2 * B_D[m])

    B_alignment = r**2 + B_quad

    result = pd.DataFrame({
        "Q": Q,
        "H": H,
        "L": L,
        "r": r,
        "c": c,
        "a": a,
        "b": b,
        "d1": d1,
        "S": S,
        "C": C,
        "E_loc": E_loc,
        "E_quad": E_quad,
        "twice_alignment": 2 * alignment,
        "entropy_regularized":
            np.abs(S + C) / np.maximum(d1, FLOOR),
        "goal_linear_regularized":
            np.abs(S) / np.maximum(d1, FLOOR),
        "alignment_regularized":
            np.abs(L) / np.maximum(Q, FLOOR),
        "fisher_norm_regularized":
            np.sqrt(H) / np.maximum(np.sqrt(Q), FLOOR),
        "B_S": B_S,
        "B_C": B_C,
        "B_D": B_D,
        "B_quad": B_quad,
        "B_r2": r**2,
        "B_alignment": B_alignment,
        "condition_a": condition_a,
        "condition_b": condition_b,
        "condition_Q": condition_Q,
        "condition_entropy_floor": condition_entropy_floor,
        "condition_alignment_floor": condition_alignment_floor,
        "condition_remainder": condition_remainder,
        "applicable": applicable,
        "error_S": np.abs(S - L),
        "error_C": np.abs(C - 0.5 * H),
        "error_D": np.abs(d1 - 0.5 * Q),
        "error_quad": np.abs(E_loc - E_quad),
        "error_alignment_quad":
            np.abs(E_quad - 2 * alignment),
        "error_alignment_exact":
            np.abs(E_loc - 2 * alignment),
    })

    pairs = BOUND_PAIRS + [
        ("error_S", "B_S"),
        ("error_C", "B_C"),
        ("error_D", "B_D"),
    ]

    for error, bound in pairs:
        result[f"ratio_{bound}"] = safe_divide(
            result[error], result[bound]
        )

        result[f"violation_{bound}"] = bound_violation(
            result[error], result[bound]
        )

        result[f"excess_{bound}"] = (
            result[error] - result[bound]
        )

    return result


def reconstruct_case(row):
    """Reconstruct a case from its archived parameters."""
    params = (
        row.alpha, row.beta, row.k, row.K, row.rho
    )

    y0 = np.array([row.x10, row.x20])
    eps = float(row.epsilon)

    aug = integrate(
        augmented_rhs,
        np.r_[y0, 1., 0., 0., 0.],
        params,
    )

    X = aug[:, :2]
    u = aug[:, 2:4]
    w = aug[:, 4:6]

    p1 = X + eps * u
    h = eps**2 * w
    p2 = p1 + h

    actual = integrate(
        rhs,
        y0 + np.array([eps, 0.]),
        params,
    )

    for array in (X, p1, p2, actual):
        if not (
            np.isfinite(array).all()
            and np.all(array > 0)
        ):
            raise ValueError(
                f"Invalid trajectory: {row.case_id}"
            )

    pos = TIMES > 0

    frame = pointwise_geometry(
        X[pos], p1[pos], h[pos], y0
    )

    frame.insert(0, "time", TIMES[pos])
    frame.insert(0, "case_id", int(row.case_id))

    dtrue = divergence(actual, X, y0)
    d1 = divergence(p1, X, y0)
    d2 = divergence(p2, X, y0)

    valid = pos & (dtrue > 1e-10)

    if not valid.any():
        raise ValueError("No usable divergence values")

    state_norm = np.mean(
        np.linalg.norm(h[pos] / y0, axis=1)
        / np.maximum(
            np.linalg.norm(eps * u[pos] / y0, axis=1),
            1e-14,
        )
    )

    reconstructed = {
        "entropy_exact":
            frame.entropy_regularized.mean(),
        "goal_linear":
            frame.goal_linear_regularized.mean(),
        "fisher_norm":
            frame.fisher_norm_regularized.mean(),
        "fisher_alignment":
            frame.alignment_regularized.mean(),
        "state_norm": state_norm,
        "D_error1": np.mean(
            np.abs(d1[valid] - dtrue[valid]) / dtrue[valid]
        ),
        "D_error2": np.mean(
            np.abs(d2[valid] - dtrue[valid]) / dtrue[valid]
        ),
        "state_error1":
            np.linalg.norm(actual - p1)
            / max(np.linalg.norm(actual), FLOOR),
    }

    summary = {"case_id": int(row.case_id)}

    for name, value in reconstructed.items():
        archived = float(getattr(row, name))

        np.testing.assert_allclose(
            value,
            archived,
            atol=ATOL,
            rtol=RTOL,
            err_msg=f"Case {row.case_id}: {name}",
        )

        summary[f"reconstruction_error_{name}"] = (
            abs(value - archived)
        )

    divergence_label = int(
        reconstructed["D_error1"] > 0.01
    )
    trajectory_label = int(
        reconstructed["state_error1"] > 0.001
    )

    assert divergence_label == int(row.label)
    assert trajectory_label == int(row.state_error1 > 0.001)

    summary["divergence_label"] = divergence_label
    summary["trajectory_label"] = trajectory_label

    frame["divergence_label"] = divergence_label
    frame["trajectory_label"] = trajectory_label

    for name in CONDITIONS:
        summary[f"fraction_{name}"] = frame[name].mean()

    for name in GEOMETRY_COLUMNS:
        values = frame[name].to_numpy()

        # Undefined values must not silently change the
        # averaging time grid for normalized diagnostics.
        if name == "c":
            defined = values[np.isfinite(values)]
            summary[f"mean_{name}"] = (
                defined.mean() if len(defined) else np.nan
            )
            summary["defined_c_count"] = len(defined)
        else:
            summary[f"mean_{name}"] = (
                values.mean()
                if np.isfinite(values).all()
                else np.nan
            )

    summary["mean_abs_c"] = (
        frame.c.abs().mean()
    )

    summary["temporal_error"] = abs(
        frame.entropy_regularized.mean()
        - 2 * frame.alignment_regularized.mean()
    )

    all_applicable = bool(frame.applicable.all())
    summary["all_times_applicable"] = all_applicable

    temporal_bound = (
        frame.B_alignment.mean()
        if all_applicable
        else np.nan
    )

    summary["temporal_bound"] = temporal_bound

    summary["temporal_violation"] = bool(
        bound_violation(
            summary["temporal_error"], temporal_bound
        )
    )

    summary["temporal_ratio"] = float(
        safe_divide(
            summary["temporal_error"], temporal_bound
        )
    )

    return frame, summary


def summarize(pointwise, cases):
    """Create descriptive summaries without pointwise inference."""
    rows = []

    def add(scope, metric, statistic, value):
        rows.append({
            "scope": scope,
            "metric": metric,
            "statistic": statistic,
            "value": value,
        })

    def describe(scope, metric, values):
        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]

        add(scope, metric, "count", len(values))

        if not len(values):
            return

        for name, value in (
            ("p05", np.percentile(values, 5)),
            ("median", np.median(values)),
            ("p95", np.percentile(values, 95)),
            ("max", np.max(values)),
        ):
            add(scope, metric, name, value)

    for name in CONDITIONS:
        add(
            "pointwise", name, "count",
            int(pointwise[name].sum()),
        )
        add(
            "pointwise", name, "fraction",
            pointwise[name].mean(),
        )

        describe(
            "configuration",
            f"fraction_{name}",
            cases[f"fraction_{name}"],
        )

    add(
        "configuration",
        "all_times_applicable",
        "count",
        int(cases.all_times_applicable.sum()),
    )

    for error, bound in BOUND_PAIRS + [
        ("error_S", "B_S"),
        ("error_C", "B_C"),
        ("error_D", "B_D"),
    ]:
        valid = (
            np.isfinite(pointwise[error])
            & np.isfinite(pointwise[bound])
        )

        describe(
            "pointwise", error,
            pointwise.loc[valid, error],
        )

        describe(
            "pointwise", f"ratio_{bound}",
            pointwise.loc[valid, f"ratio_{bound}"],
        )

        add(
            "pointwise", bound, "violations",
            int(pointwise[f"violation_{bound}"].sum()),
        )

        if valid.any():
            add(
                "pointwise", bound, "max_positive_excess",
                max(
                    0.0,
                    pointwise.loc[
                        valid, f"excess_{bound}"
                    ].max(),
                ),
            )

    for name in (
        "temporal_error",
        "temporal_bound",
        "temporal_ratio",
    ):
        describe("configuration", name, cases[name])

    add(
        "configuration",
        "temporal_bound",
        "violations",
        int(cases.temporal_violation.sum()),
    )

    columns = [
        "mean_E_loc",
        "mean_E_quad",
        "mean_twice_alignment",
        "mean_r",
    ]

    for left, right in combinations(columns, 2):
        paired = cases[[left, right]].dropna()

        add(
            "correlation",
            f"{left}_vs_{right}",
            "count",
            len(paired),
        )

        for method in ("pearson", "spearman"):
            add(
                "correlation",
                f"{left}_vs_{right}",
                method,
                paired[left].corr(
                    paired[right], method=method
                ),
            )

    for target in ("divergence", "trajectory"):
        for label in (0, 1):
            subset = cases[cases[f"{target}_label"] == label]
            scope = f"{target}_label_{label}"

            add(scope, "configurations", "count", len(subset))

            for name in (
                "mean_r",
                "mean_c",
                "mean_abs_c",
                "mean_E_loc",
                "mean_twice_alignment",
            ):
                describe(scope, name, subset[name])

    for name in cases.columns:
        if name.startswith("reconstruction_error_"):
            add(
                "reconstruction", name, "max",
                cases[name].max(),
            )

    return pd.DataFrame(rows)


def main():
    archived = pd.read_csv(
        OUT / "target_confirmation_cases.csv",
        float_precision="round_trip",
    )

    assert len(archived) == N_CASES
    assert archived.case_id.is_unique
    assert set(archived.case_id) == set(range(N_CASES))
    assert archived.split.eq("prospective_confirmation").all()
    assert not archived.isna().any().any()

    pointwise_rows = []
    case_rows = []

    for row in archived.itertuples(index=False):
        pointwise, case = reconstruct_case(row)
        pointwise_rows.append(pointwise)
        case_rows.append(case)

        if (row.case_id + 1) % 100 == 0:
            print(
                f"Reconstructed {row.case_id + 1}/{N_CASES}",
                flush=True,
            )

    pointwise = pd.concat(
        pointwise_rows, ignore_index=True
    )
    cases = pd.DataFrame(case_rows)

    assert len(pointwise) == 40000
    assert len(cases) == 1000
    assert not pointwise.duplicated(
        ["case_id", "time"]
    ).any()

    summary = summarize(pointwise, cases)

    for name, frame in (
        ("pointwise", pointwise),
        ("cases", cases),
        ("summary", summary),
    ):
        filename = OUT / f"geometry_validation_{name}.csv"
        frame.to_csv(filename, index=False)
        print(f"Wrote {filename.name}: {len(frame)} rows")

    print("\nSUFFICIENT CONDITIONS")

    for name in CONDITIONS:
        print(
            f"{name}: {pointwise[name].sum()}/"
            f"{len(pointwise)} "
            f"({pointwise[name].mean():.2%})"
        )

    print(
        "\nAll-times applicable configurations:",
        int(cases.all_times_applicable.sum()),
    )

    print("\nBOUND VIOLATIONS")

    violations = summary[
        summary.statistic == "violations"
    ]

    print(violations.to_string(index=False))

    print("\nCONFIGURATION-LEVEL CORRELATIONS")

    print(
        summary[
            summary.scope == "correlation"
        ].to_string(index=False)
    )

    print("\nMAXIMUM RECONSTRUCTION DISCREPANCIES")

    print(
        summary[
            summary.scope == "reconstruction"
        ].to_string(index=False)
    )

    if violations.value.sum() > 0:
        raise RuntimeError(
            "Bound violations detected; evidence retained."
        )

    print("\nGEOMETRY ASSESSMENT COMPLETE")


if __name__ == "__main__":
    main()
