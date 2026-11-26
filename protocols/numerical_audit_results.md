# Independent numerical consistency audit

## Scope

This audit independently verifies archived prediction
metrics, warning labels, bootstrap summaries and
principal numerical contrasts.

The verification implementation is provided in
audit_numerical_results.py, with regression tests in
tests/test_numerical_audit.py.

The audit does not import the original experimental
programs or scikit-learn scoring functions.

## Methods

- ROC AUC: Mann-Whitney rank statistic with average
  ranks for tied probabilities.
- Brier score: direct mean squared probability error.
- Logarithmic loss: direct logarithmic calculation
  with floating-point probability clipping.
- Bootstrap intervals: empirical 2.5th and 97.5th
  percentiles.
- Bootstrap standard deviations: sample standard
  deviations with one degree-of-freedom correction.
- Numerical agreement: absolute tolerance 1e-12,
  zero relative tolerance.

The audit also verifies configuration identifiers,
method coverage, finite probabilities, probability
bounds, class labels and replicate counts.

## Results

| Verification | Outcome |
|---|---|
| Original synthetic divergence labels | PASS |
| Prospective divergence and trajectory labels | PASS |
| Prospective metric records | 10/10 |
| Prospective bootstrap summaries | 5/5 |
| Replicate-level reversal identities | 5000/5000 |
| Ecological metric records | 7/7 |
| Ecological bootstrap summaries | 48/48 |
| New audit regression tests | 13/13 |
| Complete regression suite | 88/88 |
| Archive checksums before audit | 9/9 |
| Archive checksums after audit | 9/9 |

## Prospective confirmation

| Contrast | Estimate | 95% bootstrap interval |
|---|---|---|
| Divergence AUC | 0.249664 | [0.218876, 0.281867] |
| Divergence Brier | -0.156642 | [-0.170074, -0.143930] |
| Trajectory AUC | -0.156068 | [-0.183496, -0.129397] |
| Trajectory Brier | 0.102448 | [0.088554, 0.116482] |
| AUC reversal | 0.405732 | [0.354269, 0.459020] |

Both prespecified directional AUC comparisons
remain supported.

## Ecological uncertainty

The independently verified crossed-to-unstratified
vessel bootstrap interval-width ratios are:

- AUC: 2.117995.
- Brier score: 2.037948.

All four resampling schemes retain intervals
containing zero for the entropy-minus-Fisher-norm
AUC and Brier contrasts.

The ecological comparison therefore does not
establish superiority or formal equivalence.

## Reproduction commands

    python verify_archive.py
    python audit_numerical_results.py
    python -m pytest tests/test_numerical_audit.py -v
    python -m pytest tests/ -q
    python verify_archive.py

## Limitations

The audit verifies internal numerical consistency
using archived predictions and bootstrap replicates.

It does not independently regenerate every bootstrap
sample, reconstruct the original model trajectories,
refit the ecological posterior or validate physical
forecasting performance.

Original discovery trajectory-error labels cannot
be independently reconstructed from the archived
original cases.csv because that file does not
contain state_error1.

The prospective trajectory-error labels are
independently verified from the prospective
configuration archive.

The audit supplements, rather than replaces,
the clean-worktree reproduction documented in
protocols/clean_worktree_reproduction.md.
