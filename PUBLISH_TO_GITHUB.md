# Publish the project to GitHub

The connected GitHub integration could not create a new repository automatically. The intended location is **HussainAther/goal-oriented-entropy-comparison**.

1. Open https://github.com/new and choose owner **HussainAther**.
2. Enter repository name **goal-oriented-entropy-comparison**. Select **Private** for initial review, unless you are ready to release the manuscript, posterior extracts and code publicly. Do **not** initialize with a README, license or .gitignore (these files are included).
3. Download and unzip the prepared project, open a terminal inside `goal-oriented-entropy-comparison`, then run:

```bash
git init -b main
git add .
git commit -m "Initialize Computo comparative-methods submission"
git remote add origin https://github.com/HussainAther/goal-oriented-entropy-comparison.git
git push -u origin main
```

If GitHub prompts for authentication, use a credential manager or GitHub CLI. Do not put a token in the repository.

4. Check **Actions** for the `Verify archived scientific assertions` workflow. The `Main Workflows` workflow is inherited from the official Computo template and may need repository Actions permissions and Pages configuration. Verify a successful Quarto HTML/PDF build before manuscript submission.
5. When appropriate for submission, make the repository public and enable GitHub Pages. Verify the original data/code redistribution permissions, DOI metadata, author details and AI disclosure first.

This project is not automatically published by downloading this ZIP; GitHub repository creation and pushing are required.
