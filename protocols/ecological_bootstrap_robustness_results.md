# Ecological bootstrap robustness results

## Provenance

Protocol commit: 30c8efa
Implementation commit: 66d49f1

The protocol and implementation were committed and
pushed before execution of the full analysis.

This is a retrospective sensitivity analysis of
archived ecological predictions.

## Execution

- Confirmation configurations: 144
- Diagnostic methods: 7
- Archived predictions: 1008
- Resampling schemes: 4
- Replicates per scheme: 5000
- Total replicates: 20000
- Undefined replicates: 0
- Summary comparisons: 48
- Regression tests passed before execution: 45

No trajectories were regenerated and no calibration
models were refitted.

## Principal comparison

All differences are entropy_exact minus fisher_norm.

Full-sample ROC AUC difference: +0.000775
Full-sample Brier score difference: -0.002034

### ROC AUC differences

| Scheme | 95% interval | Width |
|---|---|---:|
| Unstratified vessel | [-0.014764, 0.020252] | 0.035016 |
| Stratified vessel | [-0.014286, 0.019079] | 0.033365 |
| Stratified draw | [-0.022072, 0.023961] | 0.046034 |
| Stratified crossed | [-0.032656, 0.041507] | 0.074163 |

### Brier score differences

| Scheme | 95% interval | Width |
|---|---|---:|
| Unstratified vessel | [-0.021243, 0.016677] | 0.037920 |
| Stratified vessel | [-0.017349, 0.013353] | 0.030701 |
| Stratified draw | [-0.039828, 0.023904] | 0.063732 |
| Stratified crossed | [-0.049630, 0.027649] | 0.077278 |

All principal comparison intervals include zero.

Relative to unstratified vessel resampling, crossed
resampling increases the AUC interval width by
approximately 2.12-fold and the Brier interval width
by approximately 2.04-fold.

## Additional findings

The entropy-minus-goal-linear AUC interval excludes
zero under the three simpler resampling schemes but
includes zero under crossed resampling:

[-0.001064, 0.021310].

The corresponding Brier interval remains negative:

[-0.015257, -0.003399].

Under crossed resampling, entropy retains favorable
AUC and Brier intervals against Fisher alignment
and state norm.

Comparisons against the two-feature goal-linear
plus curvature classifier remain inconclusive.

The goal-plus-curvature diagnostic is algebraically
identical to the exact entropy correction and is
not an independent competing diagnostic.

## Interpretation

The qualitative conclusion of similar observed
entropy and Fisher-norm performance survives the
alternative resampling schemes.

However, uncertainty is sensitive to the crossed
dependence structure of the ecological configurations.

The analysis does not establish superiority or
formal equivalence between these diagnostics.

## Limitations

The analysis is conditional on the fitted posterior,
archived predictions, and fixed calibration models.

Only six posterior draws per ecological group are
available.

The crossed bootstrap is a robustness analysis,
not a guarantee of nominal interval coverage.

Intervals are marginal percentile intervals and
are not adjusted for multiple comparisons.

No independent biological replication or
population-level inference is established.

## Archived evidence

code/comparator/bootstrap_robustness_replicates.csv

code/comparator/bootstrap_robustness_summary.csv

All 48 specified comparisons are retained.
