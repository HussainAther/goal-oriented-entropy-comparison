# ENT-04F — Path A submission gate (2026-11-08)

## Decision: HOLD, manuscript and analysis substantially consolidated

**Scientific contribution:** A single nonlinear ecological test bed shows that apparent entropy-warning superiority depends on the baseline. For +20 µmol/L, the posterior-draw confirmation AUC is 0.977709 (entropy), 0.896492 (component-scaled state norm), and 0.976934 (Fisher-weighted correction). Equal-information comparison therefore removes nearly all the original apparent advantage. The exact Bregman goal-gradient-plus-curvature expression is algebraically identical to the entropy correction; it is not a new method. This is a useful comparator-design caution, but significance to the statistical/machine-learning community remains an editorial question.

## Reproducibility audit

| Item | Result | Evidence |
|---|---|---|
| Perturbed ecological ODE, variations, 576 cases | PASS: direct rerun | `code/perturbation/run_ent03p.py`, `reports/full_perturbation.log` |
| Discovery reference-only features | PASS: direct rerun | `code/features/run_ent03q.py`, `reports/features.log` |
| 144 additional posterior draw confirmations | PASS: direct rerun; numeric errors identical | `code/confirmation/run_ent03v.py`, `reports/confirmation_constrained.log` |
| Matched-information features, calibration, uncertainty | PASS: direct rerun | `code/comparator/run_ent04b.py`, `reports/comparator_after_v.log` |
| Figure generation, numerical assertions | PASS | `build_release.py`, `reports/verification.json` |
| Single uninterrupted chained command | INCOMPLETE: environment timeout, not a scientific mismatch | `reports/full_pipeline_output.log` and stage logs |
| Original Bayesian fitting and archived source quantile reconstruction | PARTIAL: original reconstruction data checked; complete Stan refit not performed in this pass | `data/reproduction_full.csv`, `code/reproduce_published_curves.py` |
| Clean conda environment creation | NOT VERIFIED | `environment.yml` |
| Official journal Quarto HTML/PDF render | NOT VERIFIED: Quarto not available here | `main.qmd`; editable Word conversion rendered |
| Public GitHub, official CI, data rights and permanent archive | NOT DONE | see README |

## Venue decision: Computo candidate, not yet submission-compliant

The journal explicitly welcomes neutral comparisons showing how, when and why statistical/ML methods behave differently. Its submissions require a Quarto notebook, package-environment freezing, GitHub-based CI with accessible HTML/PDF, and code and data. We have a manuscript *source* and locally checked code results, but not an official Computo template repository/extension, clean-environment install, rendered Quarto output, or public CI. Its remit is statistics and machine learning; whether a single nonlinear-ODE diagnostic comparison provides sufficient generalizable insight for readers is uncertain. Author requirements: https://computo-journal.org/site/guidelines-authors.html ; aims: https://computo-journal.org/site/about.html .

## Content edits completed

1. One central hypothesis and contribution boundary: effect of information matching on the perceived superiority of a target-specific warning.
2. First/second-order amplitude dependence and fair comparison as the main results.
3. Biomass forecast failures confined to scope/limitations, with 23-vessel matched-case qualification; latent-state certificate experiments mentioned only as unresolved future work.
4. Exact Bregman identity explicitly flagged as existing mathematics; threshold described as exploratory; draw holdout distinguished from independent ecology.

## Non-negotiable steps before actual submission

1. Install the official Computo Python template and extension; produce local Quarto HTML and PDF from an unmodified official template layout and pass repository CI.
2. Create and verify a clean environment lock with identical end-to-end outputs on an independent machine; break the pipeline into reproducible stages as documented.
3. Resolve and cite source posterior-data redistribution permissions, deposit sharable files with permanent identifiers, and reproduce published curves from their upstream source R/Stan archive or explicitly provide reproducible archived-stage provenance.
4. Resolve authorship, conflicts/funding, source metadata, and final reference/style checks.
5. Get an editorial scope go/no-go: the computational comparison is real, but novelty and breadth are still uncertain; no new experiment is required to make this decision.

**Final status:** NOT ready to click Submit. Further work is packaging/compliance and scientific significance review, not an open-ended new simulation branch.
