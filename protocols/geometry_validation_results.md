# Retrospective information geometry assessment: results

## Provenance

Protocol commit: e2be115

Implementation commit: 3cc622f

The protocol and implementation were committed before
execution of the geometry assessment.

This analysis is retrospective and uses the existing
1,000-case prospective synthetic confirmation ensemble.
It is not an additional independent confirmation experiment.

## Execution and reconstruction

All 1,000 archived configurations were reconstructed.

All 40,000 positive sampled time points were evaluated.

All 16 new regression tests passed.

The complete regression suite contained 104 passing tests.

Local execution time was 10.208 seconds.

All nine previously archived evidence files retained
their original SHA-256 checksums.

The maximum reconstruction discrepancy across the
eight checked case-level quantities was approximately
1.388e-17.

## Applicability of the mathematical propositions

Every sufficient condition was satisfied at all
40,000 sampled positive times.

All 1,000 configurations satisfied the complete
conditions at every sampled positive time.

No violations of the six pointwise inequalities or
the temporal corollary were detected within the
protocol's numerical comparison tolerance.

## Approximation discrepancies

The following statistics describe absolute discrepancies
over 40,000 sampled positive times.

| Comparison | Median | 95th percentile | Maximum |
|---|---:|---:|---:|
| Exact versus quadratic entropy | 0.000203972 | 0.002089708 | 0.014583498 |
| Quadratic versus twice alignment | 0.000227802 | 0.003131497 | 0.042603782 |
| Exact versus twice alignment | 0.000171293 | 0.004575701 | 0.046209098 |

The median configuration-level temporal discrepancy
was approximately 0.000177, with a 95th percentile
of 0.003330 and a maximum of 0.010870.

## Analytical bound tightness

Median observed-error-to-bound ratios:

| Bound | Median ratio |
|---|---:|
| Goal-linear remainder | 0.77732895 |
| Curvature remainder | 0.63078034 |
| Divergence remainder | 0.69663203 |
| Exact-quadratic discrepancy | 0.06771825 |
| Quadratic-alignment discrepancy | 1.00000000 |
| Exact-alignment discrepancy | 0.05813384 |
| Temporal agreement | 0.04607000 |

The individual Taylor remainder bounds were substantially
less conservative than the combined normalized bounds.

The quadratic-alignment inequality was numerically
saturated to the displayed precision.

The maximum curvature remainder ratio was approximately
1.000142. This small excess passed the prespecified
absolute-plus-relative numerical tolerance.

Consequently, zero detected violations refers to the
numerical tolerance criterion, not exact floating-point
inequalities without tolerance.

## Configuration-level correlations

All correlations use 1,000 configuration-level averages.

| Comparison | Pearson | Spearman |
|---|---:|---:|
| Exact entropy versus quadratic entropy | 0.999289 | 0.998339 |
| Exact entropy versus twice Fisher alignment | 0.995392 | 0.992063 |
| Exact entropy versus relative Fisher norm | 0.602417 | 0.552306 |
| Quadratic entropy versus twice Fisher alignment | 0.998106 | 0.997468 |
| Quadratic entropy versus relative Fisher norm | 0.605756 | 0.553890 |
| Twice Fisher alignment versus relative Fisher norm | 0.611016 | 0.557935 |

These correlations are descriptive and do not establish
independent predictive validation or causal relationships.

## Geometry by warning class

Entries are medians of configuration-level averages.

| Quantity | Divergence negative | Divergence positive | Trajectory negative | Trajectory positive |
|---|---:|---:|---:|---:|
| Relative Fisher norm | 0.010680681 | 0.021338445 | 0.009638782 | 0.024990964 |
| Signed alignment cosine | 0.36590976 | 0.61461340 | 0.56557503 | 0.46222235 |
| Absolute alignment cosine | 0.37233474 | 0.61504707 | 0.56557503 | 0.46352041 |
| Exact entropy correction | 0.006804499 | 0.020370178 | 0.008662505 | 0.019553444 |
| Twice Fisher alignment | 0.006702528 | 0.020567336 | 0.008556699 | 0.019734084 |

Both relative correction magnitude and absolute alignment
were higher in the divergence-warning class.

Relative correction magnitude was higher, but absolute
alignment was lower, in the trajectory-warning class.

These descriptive patterns are consistent with a
target-dependent geometric interpretation of the
previously established diagnostic ranking reversal.

They do not establish a causal mechanism.

## Independent consistency checks

A separate calculation using the generated CSV files,
without importing the assessment implementation, verified:

- Expected record counts and unique identifiers.
- Forty sampled positive times per configuration.
- Applicability at every recorded point.
- Quadratic entropy, alignment and Fisher norm identities.
- All six pointwise inequalities within numerical tolerance.
- Agreement with the archived warning labels.
- Configuration-level diagnostic aggregation.
- The temporal inequality within numerical tolerance.

All checks passed.

This was an independent consistency check of the
generated evidence, not an independent implementation
of the ODE solver or a complete recomputation of every
summary statistic.

## Archived outputs

- geometry_validation_pointwise.csv: 40,000 data rows.
- geometry_validation_cases.csv: 1,000 data rows.
- geometry_validation_summary.csv: 268 data rows.

## Interpretation limits

The results concern sampled positive times in the existing
synthetic compartment ensemble.

They do not establish continuous-time applicability between
sampled points, universal bound tightness, generalization
to other systems, or certified nonlinear trajectory errors.

The pointwise records within each configuration are
dependent. No independence-based pointwise significance
tests or confidence intervals were calculated.
