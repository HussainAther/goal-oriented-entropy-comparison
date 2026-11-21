"""Reproducible synthetic nonlinear-compartment diagnostic benchmark.

Run from repository root: python code/compartment/run_benchmark.py
Outputs are written beside this script. The design is fixed before inspecting
comparative metrics. A degenerate warning class is reported, not tuned away.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from model import rhs, augmented_rhs

OUT = Path(__file__).resolve().parent
TIMES = np.linspace(0.0, 5.0, 41)
N_DISCOVERY, N_CONFIRMATION = 192, 144
SEEDS = {"discovery": 2026112101, "confirmation": 2026112102}
EPS_FRACTION = 0.20
WARNING_THRESHOLD = 0.01
N_BOOT = 1500
FEATURES = {
    "entropy_exact": ["entropy_exact"],
    "goal_linear": ["goal_linear"],
    "fisher_norm": ["fisher_norm"],
    "fisher_alignment": ["fisher_alignment"],
    "state_norm": ["state_norm"],
}


def integrate(fun, y0, p):
    sol = solve_ivp(lambda t, y: fun(t, y, p), (0.0, 5.0), y0,
                    t_eval=TIMES, method="DOP853", rtol=1e-10, atol=1e-12)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y.T


def divergence(y, ref, y0):
    if not (np.all(y > 0) and np.all(ref > 0)):
        raise ValueError("Generalized I-divergence requires positive states")
    return np.sum((y * np.log(y / ref) - y + ref) / y0, axis=1)


def one_case(split, case_id, rng):
    # Independent synthetic ensemble, not a fitted posterior.
    alpha = rng.uniform(0.05, 0.35)
    beta = rng.uniform(0.05, 0.35)
    k = rng.uniform(0.5, 2.0)
    K = rng.uniform(0.5, 4.0)
    rho = rng.uniform(0.05, 0.5)
    p = (alpha, beta, k, K, rho)
    x10 = rng.uniform(1.0, 10.0)
    x20 = rng.uniform(1.0, 10.0)
    y0 = np.array([x10, x20])
    eps = EPS_FRACTION * x10

    aug = integrate(augmented_rhs, np.r_[y0, 1., 0., 0., 0.], p)
    X, u, w = aug[:, :2], aug[:, 2:4], aug[:, 4:6]
    actual = integrate(rhs, y0 + np.array([eps, 0.]), p)
    p1 = X + eps * u
    h = eps**2 * w
    p2 = p1 + h
    if not all(np.all(arr > 0) for arr in (X, actual, p1, p2)):
        raise ValueError(f"Nonpositive trajectory: {split}/{case_id}")

    dtrue = divergence(actual, X, y0)
    d1 = divergence(p1, X, y0)
    d2 = divergence(p2, X, y0)
    valid = (TIMES > 0) & (dtrue > 1e-10)
    if not valid.any():
        raise ValueError(f"No usable divergence values: {split}/{case_id}")
    err1 = float(np.mean(np.abs(d1[valid] - dtrue[valid]) / dtrue[valid]))
    err2 = float(np.mean(np.abs(d2[valid] - dtrue[valid]) / dtrue[valid]))

    # Exact Bregman decomposition: D(p2,X)-D(p1,X) = S+C.
    S = np.sum(h * np.log(p1 / X) / y0, axis=1)
    C = divergence(p2, p1, y0)
    identity_error = float(np.max(np.abs(S + C - (d2 - d1))))
    if identity_error > 1e-9:
        raise AssertionError(f"Bregman identity error: {identity_error}")
    pos = TIMES > 0
    denom = np.maximum(d1[pos], 1e-12)
    fisher_num = np.sqrt(np.sum(h[pos]**2 / (p1[pos] * y0), axis=1))
    fisher_den = np.maximum(np.sqrt(np.sum((p1[pos] - X[pos])**2 / (p1[pos] * y0), axis=1)), 1e-12)
    align_num = np.abs(np.sum((p1[pos] - X[pos]) * h[pos] / (p1[pos] * y0), axis=1))
    align_den = np.maximum(np.sum((p1[pos] - X[pos])**2 / (p1[pos] * y0), axis=1), 1e-12)
    state_num = np.linalg.norm(h[pos] / y0, axis=1)
    state_den = np.maximum(np.linalg.norm(eps * u[pos] / y0, axis=1), 1e-14)

    return dict(split=split, case_id=case_id, alpha=alpha, beta=beta, k=k,
                K=K, rho=rho, x10=x10, x20=x20, epsilon=eps,
                D_error1=err1, D_error2=err2,
                label=int(err1 > WARNING_THRESHOLD),
                entropy_exact=float(np.mean(np.abs(S[pos] + C[pos]) / denom)),
                goal_linear=float(np.mean(np.abs(S[pos]) / denom)),
                fisher_norm=float(np.mean(fisher_num / fisher_den)),
                fisher_alignment=float(np.mean(align_num / align_den)),
                state_norm=float(np.mean(state_num / state_den)),
                max_identity_error=identity_error)


def generate(split, n):
    rng = np.random.default_rng(SEEDS[split])
    rows = [one_case(split, i, rng) for i in range(n)]
    return pd.DataFrame(rows)


def design(df, columns):
    return np.log10(np.maximum(df[columns].to_numpy(dtype=float), 1e-15))


def evaluate(discovery, confirmation):
    predictions, metrics = [], []
    y = confirmation.label.to_numpy()
    for method, cols in FEATURES.items():
        clf = make_pipeline(StandardScaler(), LogisticRegression(C=1., max_iter=2000))
        clf.fit(design(discovery, cols), discovery.label)
        probs = clf.predict_proba(design(confirmation, cols))[:, 1]
        metrics.append(dict(method=method, auc=float(roc_auc_score(y, probs)),
                            brier=float(brier_score_loss(y, probs))))
        predictions.extend(dict(case_id=int(i), method=method, label=int(label),
                                probability=float(prob))
                           for i, label, prob in zip(confirmation.case_id, y, probs))
    return pd.DataFrame(metrics), pd.DataFrame(predictions)


def bootstrap(predictions):
    # Independent synthetic cases are resampled by case_id; no artificial clusters.
    wide = predictions.pivot(index=["case_id", "label"], columns="method", values="probability").reset_index()
    rng = np.random.default_rng(2026112103)
    rows = []
    for j in range(N_BOOT):
        sample = wide.iloc[rng.integers(0, len(wide), size=len(wide))]
        y = sample.label.to_numpy()
        row = {"iteration": j}
        for comp in ("fisher_norm", "fisher_alignment"):
            row[f"brier_entropy_minus_{comp}"] = (brier_score_loss(y, sample.entropy_exact) -
                                                  brier_score_loss(y, sample[comp]))
            row[f"auc_entropy_minus_{comp}"] = (roc_auc_score(y, sample.entropy_exact) -
                                                roc_auc_score(y, sample[comp])) if len(np.unique(y)) == 2 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    discovery = generate("discovery", N_DISCOVERY)
    confirmation = generate("confirmation", N_CONFIRMATION)
    cases = pd.concat([discovery, confirmation], ignore_index=True)
    cases.to_csv(OUT / "cases.csv", index=False)
    for split, df in (("discovery", discovery), ("confirmation", confirmation)):
        counts = df.label.value_counts().reindex([0, 1], fill_value=0)
        print(f"{split}: negative={counts[0]}, positive={counts[1]}", flush=True)
    if discovery.label.nunique() < 2 or confirmation.label.nunique() < 2:
        print("STOP: one warning class absent; metrics are undefined. Design remains unchanged.")
        return 2
    metrics, predictions = evaluate(discovery, confirmation)
    metrics.to_csv(OUT / "metrics.csv", index=False)
    predictions.to_csv(OUT / "predictions.csv", index=False)
    bootstrap(predictions).to_csv(OUT / "bootstrap.csv", index=False)
    print(metrics.to_string(index=False))
    print("DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
