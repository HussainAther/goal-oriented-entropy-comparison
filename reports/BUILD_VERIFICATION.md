# ENT-04G build verification (2026-11-08)

- Numerical assertions: PASS. `python build_release.py` returned success; saved `reports/ent04g_release_assertions.log` and `reports/verification.json`.
- Bibliography: PASS. All manuscript citation keys resolve to `references.bib`; citations are processed by Pandoc for review output.
- Manuscript review previews: PASS. Pandoc produced `manuscript/ENT04G_review.html`, `manuscript/ENT04G_review.pdf`, and `manuscript/ENT04G_review.docx`. Word version was converted into six page PNGs and layout reviewed. Pandoc direct render and official Quarto render are not equivalent.
- Full source-to-results stages: reproduced individually in ENT-04F. Not rerun anew for ENT-04G because changes are manuscript/bibliography/provenance only. The prior chained command timed out; uninterrupted clean build is not verified.
- Public data licensing: Dryad's general published-data policy is CC0; original author R scripts may have separate software licenses. See `DATA_PROVENANCE_AND_LICENSE.md`.
- Environment isolation: FAIL / NOT ATTEMPTED. Conda / micromamba unavailable in this session; `environment.yml` is a target manifest, not a tested lock. Current versions recorded in `reports/environment_current.json`.
- Official Computo rendering: NOT VERIFIED. Quarto and official template absent from this environment; no public GitHub Actions/Pages deployment.
- Editorial scope: NOT PREAPPROVED. The study compares ODE reliability diagnostics on one model, and the journal is oriented to stats/ML computational work.
- Personal submission metadata: corresponding address, ORCID, funding/conflicts and explicit code license require author approval.

Decision: HOLD for official Computo submission, but the manuscript and review deliverables are prepared. Avoid more numerical experiments until the reproducibility/template gates are cleared.
