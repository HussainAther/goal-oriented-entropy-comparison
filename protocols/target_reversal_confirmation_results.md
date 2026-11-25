# Prospective target-reversal confirmation: results

Execution date: 2026-11-25

## Provenance

- Protocol commit: ef3352f
- Implementation commit: 1215af8
- Regression tests: 27 passed before execution
- Discovery configurations: 192
- Fresh confirmation configurations: 1000
- Confirmation seed: 2026112301
- Bootstrap seed: 2026112302
- Bootstrap replicates: 5000
- Undefined bootstrap replicates: 0

## Primary outcome

Both prespecified directional hypotheses were supported.

| Comparison | Estimate | 95% CI lower | 95% CI upper |
|---|---:|---:|---:|
| Divergence AUC difference | 0.249664 | 0.218876 | 0.281867 |
| Trajectory AUC difference | -0.156068 | -0.183496 | -0.129397 |
| Ranking-reversal contrast | 0.405732 | 0.354269 | 0.459020 |

All differences are entropy minus Fisher norm, except
the reversal contrast, which is the divergence AUC
difference minus the trajectory AUC difference.

## Secondary outcome

| Comparison | Estimate | 95% CI lower | 95% CI upper |
|---|---:|---:|---:|
| Divergence Brier difference | -0.156642 | -0.170074 | -0.143930 |
| Trajectory Brier difference | 0.102448 | 0.088554 | 0.116482 |

Lower Brier scores indicate better probabilistic predictions.

## Diagnostic performance

| Target | Diagnostic | ROC AUC | Brier |
|---|---|---:|---:|
| Divergence | entropy_exact | 0.999572 | 0.032246 |
| Divergence | goal_linear | 0.999446 | 0.032902 |
| Divergence | fisher_norm | 0.749908 | 0.188888 |
| Divergence | fisher_alignment | 0.996564 | 0.037252 |
| Divergence | state_norm | 0.759143 | 0.183283 |
| Trajectory | entropy_exact | 0.810527 | 0.176724 |
| Trajectory | goal_linear | 0.808359 | 0.177675 |
| Trajectory | fisher_norm | 0.966595 | 0.074276 |
| Trajectory | fisher_alignment | 0.807718 | 0.178020 |
| Trajectory | state_norm | 0.966057 | 0.075792 |

## Confirmation class counts

- Divergence: 608 positive, 392 negative.
- Trajectory: 581 positive, 419 negative.

## Protocol adherence

The prescribed sample size, seeds, warning thresholds,
diagnostic features, and calibration procedure were retained.

The experiment completed successfully without reported
numerical failures or replacement configurations.

Both primary confidence intervals excluded zero in their
prespecified directions.

## Interpretation

The warning-target-dependent ranking reversal observed
in the exploratory analysis was prospectively confirmed
in fresh synthetic configurations.

Goal-oriented diagnostics performed better for divergence
error warnings, whereas norm-based diagnostics performed
better for trajectory error warnings.

The close agreement between exploratory and prospective
effect estimates supports the stability of this finding
within the specified synthetic generating distribution.

## Limitations

The warning thresholds were inherited from earlier analyses.
The bootstrap intervals are conditional on the original
discovery sample and fitted calibration models.

Both warning targets share the same underlying trajectories.

The experiment does not establish universal diagnostic
superiority, external biological validity, or independent
replication by another research group.

## Archived evidence

The following files are retained under code/compartment/:

- target_confirmation_cases.csv
- target_confirmation_metrics.csv
- target_confirmation_predictions.csv
- target_confirmation_counts.csv
- target_confirmation_bootstrap.csv
- target_confirmation_summary.csv
