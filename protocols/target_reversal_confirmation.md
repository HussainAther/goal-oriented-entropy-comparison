# Prospective confirmation of warning-target-dependent diagnostic rankings

## 1. Motivation

An exploratory analysis of the nonlinear compartment benchmark
identified a reversal in the relative discrimination performance
of entropy correction and Fisher-weighted norm when the warning
target changed from divergence error to state-trajectory error.

The trajectory-error threshold was selected after inspecting the
original discovery and confirmation error distributions.
Consequently, the original finding is exploratory.

This protocol specifies a prospective confirmation experiment
using previously unexamined synthetic configurations.

## 2. Research question

Does the entropy-versus-Fisher-norm ranking reversal persist
in fresh configurations drawn from the same synthetic
parameter-generating distribution?

## 3. Hypotheses

Let Delta_D denote entropy-minus-Fisher-norm ROC AUC for
divergence-error warnings.

Let Delta_S denote entropy-minus-Fisher-norm ROC AUC for
state-trajectory-error warnings.

The directional hypotheses are:

- H1: Delta_D > 0.
- H2: Delta_S < 0.

The reversal is supported if both paired 95% bootstrap
intervals exclude zero in their hypothesized directions.

The contrast Delta_D - Delta_S will also be estimated using
paired resampling of the same confirmation configurations.

## 4. Fixed experimental design

- Model: existing nonlinear two-compartment system.
- Parameter-generating distributions: unchanged.
- Discovery sample: original 192 synthetic configurations.
- Fresh confirmation sample: 1000 configurations.
- Confirmation random seed: 2026112301.
- Bootstrap random seed: 2026112302.
- Bootstrap replicates: 5000.
- Initial perturbation: 20% of the first compartment.
- Time interval: [0, 5].
- Evaluation grid: 41 equally spaced time points.
- Integrator: DOP853.
- Relative tolerance: 1e-10.
- Absolute tolerance: 1e-12.

Independent uniform sampling distributions:

| Quantity | Lower bound | Upper bound |
|---|---:|---:|
| alpha | 0.05 | 0.35 |
| beta | 0.05 | 0.35 |
| k | 0.5 | 2.0 |
| K | 0.5 | 4.0 |
| rho | 0.05 | 0.5 |
| Initial x1 | 1.0 | 10.0 |
| Initial x2 | 1.0 | 10.0 |

The existing numerical model and diagnostic definitions
will be reused without modification.

## 5. Warning targets

Divergence-error warning:

    D_error1 > 0.01

Trajectory-error warning:

    state_error1 > 0.001

Both thresholds are inherited from previous analyses.
They are fixed before generating the fresh confirmation
sample but were not selected independently of all prior data.

The original exploratory results will remain unchanged.

## 6. Diagnostic features

All five existing diagnostic features will be evaluated:

- entropy_exact
- goal_linear
- fisher_norm
- fisher_alignment
- state_norm

The primary comparison is entropy_exact versus fisher_norm.
Other comparisons are secondary and descriptive.

## 7. Calibration

Separate logistic classifiers will be fitted for each
diagnostic and warning target.

Calibration will use only the original 192 discovery cases.

The existing calibration procedure will be retained:

- Feature floor: 1e-15.
- Base-10 logarithmic transformation.
- StandardScaler fitted on discovery data.
- LogisticRegression with C=1 and max_iter=2000.

No fresh confirmation labels will be used for fitting,
feature selection, threshold selection, or hyperparameter tuning.

## 8. Statistical analysis

Primary metric: ROC AUC.

Secondary metric: Brier score.

Report all diagnostic metrics for both warning targets,
including warning prevalence and class counts.

Use paired nonparametric bootstrap resampling of the
1000 fresh confirmation configurations.

Each bootstrap replicate will use the same sampled case
indices for both warning targets and all diagnostics.

Report percentile intervals at the 2.5th and 97.5th
percentiles for:

- Delta_D.
- Delta_S.
- Delta_D - Delta_S.
- Entropy-minus-Fisher-norm Brier differences for both targets.

The classifiers remain fixed during resampling.
Intervals are conditional on the original discovery sample
and fitted calibration models.

Replicates lacking either warning class will have undefined
AUC recorded as missing. Their number will be reported.

## 9. Decision and reporting rules

The sample size, seeds, thresholds, generating distributions,
and primary comparison will not be altered in response to
observed confirmation results.

Report results regardless of whether the hypotheses are supported.

Numerical failures and nonpositive trajectories must be
documented. Failed configurations must not be silently
replaced or excluded.

If a warning target has only one class, its AUC is undefined.
The experimental design will not be tuned to eliminate
this outcome.

Any necessary implementation corrections will be documented
with their rationale and associated Git history.

## 10. Scope of inference

This experiment is a prospective confirmation within the
same synthetic model and parameter-generating distribution.

It is not an independent replication by another research
group, validation on observational data, or evidence of
universal diagnostic superiority across dynamical systems.

The protocol is version-controlled before execution.
It is not a claim of formal third-party preregistration.
