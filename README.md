# Goal-oriented relative-entropy warning comparison

Source repository for the article **When Goal-Oriented Error Warnings Agree: A Reproducible Comparison of Relative-Entropy and Information-Weighted Diagnostics**. Based on the official [Computo Python template](https://github.com/computorg/template-computo-python).

## Scientific result

The entropy-correction indicator discriminated internal first-order relative-entropy approximation errors with ROC AUC **0.977709** in 144 posterior-draw confirmation configurations. The matched Fisher-weighted indicator achieved **0.976934**, while a less-informed component-scaled state norm achieved **0.896492**. This is a neutral comparison, *not* a claim that a new entropy estimator or mathematical theorem was discovered. The confirmation draws come from a posterior fitted using all observed chemostats; this is **not** independent biological validation.

## Reproduction

Create the environment using `conda env create -f environment.yml`, activate `entropy-comparator`, and run these stages sequentially from the repository root:

```bash
./run_stages.sh perturbation
./run_stages.sh features
./run_stages.sh confirmation
./run_stages.sh comparator
./run_stages.sh figures
```

These stages regenerate approximation experiments (576 configurations), reference-only diagnostic features, posterior-draw confirmation (144 configurations), matched warning comparisons, figures, and recorded assertions. The original R/Stan posterior fitting and published-curve reproduction require the upstream study archive and are **not rerun by this Python pipeline**. Existing CSV outputs are supplied to ease inspection, not as substitutes for source execution.

## Quarto / Computo

The repository keeps Computo's `_quarto.yml`, `main.qmd`, `README.qmd` and upstream reusable publishing workflow. Install Quarto and run `quarto add computorg/computo-quarto-extension`, then `quarto render` to produce official layouts. **The official build and clean CI still require verification after publishing this repository.** Review PDFs in `manuscript/` are not Computo-produced PDFs.

## Provenance, licensing and limitations

The original data are available at [Dryad DOI 10.5061/dryad.51c59zwj5](https://doi.org/10.5061/dryad.51c59zwj5) and distributed under CC0. Source data, posterior-export provenance and reuse limits are documented in [DATA_PROVENANCE_AND_LICENSE.md](DATA_PROVENANCE_AND_LICENSE.md). Article text is intended for CC BY 4.0; do not assume the research code or third-party scripts have identical licensing without checking. See [RELEASE_GATE.md](RELEASE_GATE.md) for the outstanding build and submission requirements.

## AI assistance disclosure

AI was used in analysis-code development, paper drafting and checking. The author is responsible for all scientific results and verification. The manuscript contains an explicit disclosure.
