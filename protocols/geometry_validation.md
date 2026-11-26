# Retrospective empirical assessment of local information geometry

## Objective

Assess the empirical relevance of the local information-geometric
bounds established in the manuscript and quantify the relationship
between exact entropy correction, quadratic entropy correction,
Fisher alignment and relative Fisher correction magnitude.

## Analysis status

This is a retrospective explanatory analysis of an existing
prospective confirmation ensemble.

It is not an independent confirmation experiment, a new prospective
test, or a formal preregistration.

The predictive performance of the existing diagnostics is already known.
The new pointwise geometry summaries have not yet been computed.

## Input evidence

code/compartment/target_confirmation_cases.csv

Use all 1,000 archived prospective confirmation configurations.

Retain the archived model parameters, initial conditions, perturbation
amplitudes and warning targets.

No new random configurations will be generated.

No classifiers will be refitted.

## Numerical reconstruction

Reuse the existing compartment model and integration implementation.

Reconstruct the reference state X, first variation u and second
variation w using the original augmented differential equations.

Use the original 41-point time grid on [0,5], DOP853 integration,
relative tolerance 1e-10 and absolute tolerance 1e-12.

Define

    p1 = X + epsilon*u
    h  = epsilon**2*w
    p2 = p1 + h
    q  = p1 - X

Use the original initial-state component weights.

Analyze the 40 strictly positive sampled times per configuration,
giving 40,000 expected pointwise records.

## Reconstruction checks

Require unique case identifiers and complete input records.

Check positivity and finite numerical values.

Reconstruct the original five diagnostic averages and compare them
against the archived case-level values.

Use absolute tolerance 1e-10 and relative tolerance 1e-7 for
reconstruction comparisons.

Reconstruct the nonlinear perturbed trajectory and verify the
archived divergence and trajectory errors and warning labels.

Report reconstruction discrepancies.

Stop on a failed reconstruction check rather than replacing cases.

## Pointwise geometry

For each sampled positive time, compute

    G = diag(1 / (omega*p1))
    Q = q.T @ G @ q
    H = h.T @ G @ h
    L = q.T @ G @ h

    r = sqrt(H/Q)
    c = L/sqrt(Q*H)
    A_F = abs(L)/Q

    E_quad = abs(2*L + H)/Q

Compute the exact normalized entropy correction using the original
Bregman decomposition.

Record both the implemented regularized diagnostics and the
unregularized mathematical quantities where defined.

When H is zero, signed alignment is undefined but A_F is zero.
Do not replace undefined signed alignment with an artificial value.

When Q is zero, normalized mathematical quantities are undefined.

## Sufficient conditions

Compute

    a = max(abs(q)/p1)
    b = max(abs(h)/p1)

Record separately whether

    a < 1
    b < 1
    Q > 0
    d1 >= 1e-12
    Q >= 1e-12

Compute the manuscript remainder bounds B_S, B_C and B_D only
where their assumptions hold.

Record whether

    2*B_D < Q

The normalized entropy-alignment proposition is applicable only
when all its stated sufficient conditions hold.

Failure of a sufficient condition does not imply failure of the
underlying approximation.

## Bound evaluation

Where applicable, compute

    B_quad =
        (2*(B_S+B_C) + 2*E_quad*B_D)/(Q-2*B_D)

    B_alignment = r**2 + B_quad

Evaluate the observed discrepancies

    error_quad = abs(E_loc-E_quad)
    error_alignment_quad = abs(E_quad-2*A_F)
    error_alignment_exact = abs(E_loc-2*A_F)

Check the corresponding mathematical inequalities.

Use numerical comparison tolerance

    1e-10 + 1e-8*abs(bound)

Report the number of violations and maximum numerical excess.

Do not silently remove violations.

For strictly positive bounds, report observed-error-to-bound
ratios as descriptive measures of bound tightness.

## Configuration-level aggregation

Compute the fraction of sampled positive times satisfying each
condition and all sufficient conditions jointly.

Report the number of configurations satisfying the complete
proposition at every sampled positive time.

For those configurations, evaluate the temporal corollary.

Compute time-averaged exact entropy correction, quadratic entropy
correction, Fisher alignment and relative Fisher norm.

Report the absolute discrepancy between the averaged entropy
correction and twice the averaged Fisher alignment.

## Descriptive analysis

Report counts, medians, 5th and 95th percentiles, and maxima
for applicable approximation discrepancies and bound ratios.

Report Pearson and Spearman correlations across configurations
between exact entropy correction, quadratic entropy correction,
twice Fisher alignment and relative Fisher norm.

Summarize relative correction magnitude and signed alignment
by the two existing warning targets separately.

Use the original strict warning thresholds:

    divergence error > 0.01
    trajectory error > 0.001

No threshold optimization will be performed.

Pointwise summaries are descriptive. Sampled times within a
configuration are not independent observations.

No pointwise significance tests or independence-based confidence
intervals will be calculated.

## Outputs

code/compartment/geometry_validation_pointwise.csv
code/compartment/geometry_validation_cases.csv
code/compartment/geometry_validation_summary.csv

Preserve the original archived evidence without modification.

Verify the existing archive checksums before and after execution.

## Interpretation limits

The analysis assesses the local geometry of the existing synthetic
benchmark. It does not establish universal bound tightness,
causal explanations of predictive performance, generalization to
other dynamical systems, or independently validated forecasting.

The mathematical bounds concern diagnostic approximation, not
certified bounds on the exact nonlinear trajectory error.
