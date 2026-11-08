# ENT-04G submission release gate

Status: HOLD (not yet formally submittable to Computo)

## Passed in this release
- Downstream poster-draw computation stages previously rerun separately, see archived audit reports.
- Metrics and plots regenerate from packaged stage outputs, with exact numerical assertions.
- BibTeX source added, citation keys validated.
- Source Dryad data CC0 reuse policy checked and provenance statement supplied.
- Review HTML/PDF generated with Pandoc + LaTeX; these are previews, NOT Quarto/Computo outputs.

## Required before submitting
1. Copy content into the current official Computo Python template (including its **unmodified** official CI workflow), verify `quarto render` produces HTML and PDF locally, then verify public GitHub Pages and CI. This environment has no Quarto, DNS network access, or micromamba.
2. Recreate the environment from pinned package versions on a clean machine and rerun all four science stages; this environment differs from environment.yml and its chained run timed out in previous testing.
3. Complete upstream dataset-to-posterior-export provenance or transparently state the fixed-source draw-table dependency in the public software/data-availability statement.
4. Author confirms affiliations, code licensing, funding/conflicts, corresponding address, ORCID and release approval.
5. Editorial scope determination for a single-ecological-model neutral comparator study; article novelty is bounded.
