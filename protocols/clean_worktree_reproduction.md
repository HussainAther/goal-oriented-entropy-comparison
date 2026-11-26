# Clean-worktree computational reproduction

## Scope

This report documents same-environment reproduction of
the prospective synthetic target-reversal confirmation
and ecological bootstrap robustness analyses.

The reproduction was conducted in a detached Git worktree
at commit 39c19aa, separate from the primary repository.

No archived evidence in the primary repository was
overwritten.

## Computational environment

| Component | Version |
|---|---|
| Operating system | macOS 27.2.1 |
| Architecture | arm64 |
| Python | 3.11.18 |
| NumPy | 2.1.3 |
| SciPy | 1.15.1 |
| pandas | 2.2.3 |
| scikit-learn | 1.6.1 |
| Matplotlib | 3.10.0 |
| pytest | 9.2.0 |
| Jupyter | 1.1.1 |
| ipykernel | 7.4.0 |
| tabulate | 0.10.0 |
| Quarto | 1.11.9 |
| Git | 2.54.0 (Apple Git-157) |

## Procedure

1. Create a detached worktree at commit 39c19aa.
2. Verify the nine archived evidence files.
3. Execute the complete regression test suite.
4. Regenerate the prospective synthetic confirmation.
5. Verify archive integrity.
6. Regenerate the ecological bootstrap robustness analysis.
7. Verify archive integrity and inspect Git modifications.

Commands executed within the isolated worktree:

    python verify_archive.py
    python -m pytest tests/ -q
    python -u code/compartment/run_target_confirmation.py
    python verify_archive.py
    python -u code/comparator/run_bootstrap_robustness.py
    python verify_archive.py
    git diff --stat
    git status --short

## Results

| Verification | Outcome |
|---|---|
| Regression tests | 75 passed |
| Initial archive verification | 9/9 passed |
| Prospective synthetic regeneration | Byte-identical |
| Ecological bootstrap regeneration | Byte-identical |
| Final archive verification | 9/9 passed |
| Final worktree modifications | None |

### Prospective synthetic confirmation

- Fresh confirmation configurations: 1,000.
- Paired bootstrap replicates: 5,000.
- Undefined bootstrap replicates: 0.
- All six regenerated output files matched their
  archived SHA-256 checksums.
- Execution time: 13.510 seconds.

The reproduced principal AUC contrasts were:

| Contrast | Estimate | 95% bootstrap interval |
|---|---:|---|
| Divergence error | 0.249664 | [0.218876, 0.281867] |
| Trajectory error | -0.156068 | [-0.183496, -0.129397] |
| Reversal contrast | 0.405732 | [0.354269, 0.459020] |

Both prespecified directional comparisons remained
supported.

### Ecological bootstrap robustness

- Confirmation configurations: 144.
- Resampling schemes: 4.
- Bootstrap replicates per scheme: 5,000.
- Total bootstrap replicates: 20,000.
- Undefined bootstrap replicates: 0.
- Both regenerated output files matched their
  archived SHA-256 checksums.
- Execution time: 42.335 seconds.

The reproduced entropy-minus-Fisher-norm AUC
intervals were:

| Scheme | Lower | Upper |
|---|---:|---:|
| Unstratified vessel | -0.014764 | 0.020252 |
| Stratified vessel | -0.014286 | 0.019079 |
| Stratified draw | -0.022072 | 0.023961 |
| Stratified crossed | -0.032656 | 0.041507 |

The archived ecological predictions remained unchanged.

## Interpretation and limitations

The experiment demonstrates byte-identical regeneration
of eight generated output files from the committed
source code in a separate local worktree.

The ninth verified file is the archived ecological
prediction input to the bootstrap analysis.

This is distinct from the previously documented isolated
ecological trajectory and prediction regeneration,
which used numerical probability tolerances.

The present verification uses the same local software
environment and hardware architecture. It does not
establish cross-platform byte identity, independent
third-party replication, upstream posterior refitting,
or independent physical model validation.

Execution times are descriptive local observations,
not formal computational performance benchmarks.
