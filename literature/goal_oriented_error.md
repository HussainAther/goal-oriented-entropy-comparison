# Goal-oriented and relative-entropy error estimation

## Review status

Initial targeted review of publication records
and abstracts. Detailed theorem-level comparison
remains pending.

## Cao and Petzold (2004)

A Posteriori Error Estimation and Global Error
Control for Ordinary Differential Equations
by the Adjoint Method.

DOI: 10.1137/S1064827503420969

Full text:
https://people.cs.vt.edu/~ycao/publication/siam_42096.pdf

### Established contribution

Adjoint-based estimation and control of global
numerical approximation errors in ODE initial
value problems.

### Relevance and distinction

Our benchmarks also involve ODE initial value
problems. However, our principal approximation
error arises from truncating the initial-condition
Taylor expansion, rather than from integration
discretization.

Our diagnostics do not replace adjoint-based
global numerical error estimation.

Existing bibliography key: cao2004.

## Endtmayer, Langer and Wick (2020)

Two-Side a Posteriori Error Estimates for
the Dual-Weighted Residual Method.

DOI: 10.1137/18M1227275

### Established contribution

Two-sided error estimation for nonlinear PDEs
and nonlinear goal functionals, including
analysis of remainder terms.

### Relevance and distinction

Provides context for nonlinear goal-oriented
error estimation and the importance of
controlling neglected terms.

Our local Taylor remainder inequalities bound
the approximation of the entropy diagnostic.
They are not certified estimates of the complete
model trajectory error.

## Giesselmann, Makridakis and Pryer (2015)

A Posteriori Analysis of Discontinuous Galerkin
Schemes for Systems of Hyperbolic Conservation Laws.

DOI: 10.1137/140970999

Open-access preprint:
https://arxiv.org/abs/1405.7616

### Established contribution

Reliable a posteriori error estimation using
solution reconstruction and relative entropy
stability for nonlinear conservation laws.

### Relevance and distinction

Establishes important precedent for using
relative entropy in numerical error estimation.

Our work instead compares normalized warning
features derived from shared first- and
second-order perturbation information.

Our weighted I-divergence is not assumed to
satisfy the entropy stability structure used
in conservation-law analysis.

## Implications for manuscript positioning

Do not claim:

- The first use of relative entropy for error
  estimation.
- The invention of goal-oriented error analysis.
- A certified a posteriori estimator.
- Universal superiority over norm-based methods.
- A replacement for adjoint-based methods.

Emphasize:

- Matched-information diagnostic comparison.
- Controlled local geometric interpretation.
- Exact counterexamples to diagnostic ordering.
- Prospective confirmation of target-dependent
  predictive ranking reversal.
- Dependence-aware uncertainty analysis.
- Executable and reproducible benchmarks.

## Remaining uncertainty

These references establish important precedents
but do not constitute an exhaustive novelty search.

More directly related comparisons of diagnostic
ranking under alternative approximation-error
targets must still be investigated.
