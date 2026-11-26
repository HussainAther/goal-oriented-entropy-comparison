from pathlib import Path
import csv
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]

class ReleaseDataChecks(unittest.TestCase):
    def test_core_files(self):
        for name in ["main.qmd", "references.bib", "run_stages.sh", "code/comparator/run_ent04b.py", "figures/matched_auc.png"]:
            self.assertTrue((ROOT/name).is_file(), name)

    def test_comparator_claims(self):
        p = ROOT / "code/comparator/matched_metrics.csv"
        with p.open(newline="") as f:
            rows = list(csv.DictReader(f))
        values = {r["method"]: float(r["auc"]) for r in rows}
        self.assertAlmostEqual(values["entropy_exact"], 0.9777088583058733, places=10)
        self.assertAlmostEqual(values["fisher_norm"], 0.9769335142469471, places=10)
        self.assertAlmostEqual(values["state_norm"], 0.8964915681333592, places=10)

    def test_reference_integrity(self):
        qmd = (ROOT/"main.qmd").read_text()
        bib = (ROOT/"references.bib").read_text()
        used = set(re.findall(r"@([a-zA-Z][a-zA-Z0-9_-]+)", qmd))
        defined = set(re.findall(r"@(?:article|misc|book|inproceedings|techreport)\s*\{([^,]+)", bib))
        # Quarto cross-references are distinct from
        # bibliographic citations.
        document_labels = set(re.findall(
            r"(?m)^\s*#\|\s*label:\s*([a-zA-Z][a-zA-Z0-9_-]+)\s*$",
            qmd,
        ))

        # Mathematical statements use Quarto fenced divs,
        # rather than executable-cell label declarations.
        fenced_div_labels = set(re.findall(
            r"(?m)^\s*:{3,}\s*\{#([a-zA-Z][a-zA-Z0-9_-]+)(?=[\s}])",
            qmd,
        ))

        document_labels.update(fenced_div_labels)

        unresolved = used - defined - document_labels

        self.assertFalse(
            unresolved,
            f"Unresolved citations or document references: {unresolved}",
        )

if __name__ == "__main__": unittest.main()
