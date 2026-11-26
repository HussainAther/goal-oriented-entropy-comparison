# Ecological bootstrap robustness analysis

## Objective

Assess the sensitivity of diagnostic performance
comparisons to the resampling structure used for
uncertainty estimation in the ecological benchmark.

This is a retrospective statistical sensitivity analysis,
not an independent confirmation experiment.

## Input evidence

code/comparator/matched_predictions.csv

The archived predictions comprise 144 configurations
and seven diagnostic methods, yielding 1008 predictions.

Each ecological group contains a complete crossed
design of 12 vessels and six posterior parameter draws.

There are no missing values or duplicate predictions.

## Fixed analysis conditions

- Preserve the original warning labels.
- Preserve the original calibrated probabilities.
- Do not refit calibration models.
- Do not regenerate trajectories.
- Do not change diagnostic definitions.
- Preserve the original archived bootstrap results.

## Resampling schemes

1. Unstratified vessel-cluster bootstrap:
   sample 24 vessel clusters with replacement from
   the combined ecological groups.

2. Stratified vessel-cluster bootstrap:
   independently sample 12 vessels with replacement
   within each group, retaining all six draws.

3. Stratified posterior-draw-cluster bootstrap:
   independently sample six posterior draws with
   replacement within each group, retaining all vessels.

4. Stratified crossed bootstrap:
   independently sample 12 vessels and six posterior
   draws with replacement within each group.
   Retain their Cartesian product, including
   replacement multiplicities.

Each scheme produces 144 resampled configurations
per replicate.

## Computational specification

- Bootstrap replicates per scheme: 5000
- Total bootstrap replicates: 20000
- Unstratified vessel seed: 2026112601
- Stratified vessel seed: 2026112602
- Stratified draw seed: 2026112603
- Stratified crossed seed: 2026112604
- Confidence intervals: percentile, 95 percent

All diagnostic methods must use identical resampled
configuration indices within each replicate.

## Performance measures

Primary: ROC AUC differences.

Secondary: Brier score differences.

All contrasts are entropy_exact minus the comparator.

Report comparisons against:

- goal_linear
- goal_plus_curvature
- fisher_norm
- fisher_alignment
- state_norm
- goal_linear_plus_curvature

The final comparator uses two features and is therefore
not a strictly matched-dimensionality comparison.

The principal scientific comparison is entropy_exact
versus fisher_norm.

## Reporting

For every scheme, comparator, and performance measure,
report:

- Original full-sample point estimate
- Bootstrap mean
- Bootstrap standard deviation
- 95 percent percentile confidence interval
- Confidence interval width
- Number of valid replicates
- Number of undefined replicates

Single-class resamples have undefined ROC AUC.
Retain their Brier score results.

Do not replace undefined replicates.

Report all specified comparisons regardless of outcome.

## Interpretation and limitations

The analysis evaluates sensitivity to alternative
resampling assumptions.

The intervals are conditional on the fitted ecological
posterior, archived predictions, and discovery-trained
calibration models.

Only six posterior draws per ecological group are
available. Draw-cluster and crossed bootstrap intervals
may therefore be unstable.

The analysis does not establish independent biological
replication, population-level uncertainty, or formal
equivalence between diagnostics.

Failure to exclude zero does not establish equivalence.

The crossed bootstrap is a sensitivity analysis rather
than a guarantee of nominal confidence interval coverage.
