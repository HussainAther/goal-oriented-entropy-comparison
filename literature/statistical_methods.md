# Statistical evaluation and dependence literature

## Owen (2007)

Art B. Owen.
The Pigeonhole Bootstrap.
Annals of Applied Statistics, 1(2), 386--411.
DOI: 10.1214/07-AOAS122

https://arxiv.org/abs/0712.1111

### Relevance

Separate resampling of crossed factors motivates
our independent vessel and posterior-draw
resampling within ecological treatment groups.

The published theoretical results concern
specific crossed random-effects settings.

They do not automatically establish nominal
coverage for our nonlinear AUC statistic.

## Owen and Eckles (2012)

Art B. Owen and Dean Eckles.
Bootstrapping Data Arrays of Arbitrary Order.
Annals of Applied Statistics, 6(3), 895--927.
DOI: 10.1214/12-AOAS547

### Relevance

Independent resampling of factor levels produces
multiplicative observation weights.

Our crossed bootstrap similarly retains the
Cartesian product and multiplicities of sampled
vessels and posterior draws.

### Limitations

The ecological configurations originate from
a shared fitted posterior, with only six
confirmation draws per treatment group.

The bootstrap is conditional on the archived
predictions and fitted calibrators.

It does not establish population-level biological
uncertainty or guaranteed interval coverage.

## DeLong et al. (1988)

Elizabeth R. DeLong, David M. DeLong,
and Daniel L. Clarke-Pearson.

Comparing the Areas Under Two or More Correlated
Receiver Operating Characteristic Curves:
A Nonparametric Approach.

Biometrics, 44(3), 837--845.
DOI: 10.2307/2531595

https://pubmed.ncbi.nlm.nih.gov/3203132/

### Relevance

Diagnostic AUCs computed from identical cases
are statistically dependent.

Our comparisons preserve case pairing across
diagnostic methods in every bootstrap replicate.

We do not implement the DeLong variance estimator.

## Gneiting and Raftery (2007)

Tilmann Gneiting and Adrian E. Raftery.

Strictly Proper Scoring Rules,
Prediction, and Estimation.

Journal of the American Statistical Association,
102(477), 359--378.
DOI: 10.1198/016214506000001437

### Relevance

The binary Brier loss is a strictly proper
probabilistic scoring rule.

We report the mean squared difference between
warning probability and binary outcome.

Lower Brier loss indicates better probabilistic
prediction under this scoring rule.

Brier score is not a pure measure of calibration.
It complements AUC, which measures discrimination.

## Implications for our experiments

- Compare diagnostic methods on identical cases.
- Preserve paired predictions during resampling.
- Preserve dependence across ecological factors.
- Distinguish discrimination from probability quality.
- Report conditional uncertainty explicitly.
- Do not interpret intervals containing zero
  as proof of diagnostic equivalence.
- Do not interpret bootstrap intervals as
  guarantees of biological external validity.

## Ecological robustness findings

Entropy-minus-Fisher-norm AUC difference:

    Estimate: +0.000775

Unstratified vessel bootstrap:

    95% interval: [-0.014764, +0.020252]

Stratified crossed bootstrap:

    95% interval: [-0.032656, +0.041507]

Interval width ratio: approximately 2.12.

Corresponding Brier interval width ratio:
approximately 2.04.

All four resampling schemes include zero for
the principal AUC and Brier contrasts.

These results support caution concerning
superiority claims, not formal equivalence.

## Intended manuscript citations

- Methods: discrimination and proper scoring rules.
- Methods: paired and crossed resampling design.
- Results: ecological uncertainty robustness.
- Discussion: conditional statistical inference.
