# Changelog

## [2025–2026] – ANIS Upgrade
- >60% of labs updated or new.
- Added Generative AI, Federated Learning, GNNs, Digital Twins.
- Introduced Docker/Devcontainer, MLflow, GitHub Actions CI.
- Split Lab 6 into 6a (WGCNA) and 6b (Networks + GNN).
- Linked projects to Oncohelp research and dissertations.
- **Governance & Documentation:**
  - Created `RESOURCES.md` (textbooks, articles, online tools).
  - Created `POLICIES.md` (governance, GDPR, GA4GH, sustainability).
  - Created `CONTRIBUTING.md` (pair work, roles A/B, PR workflow).
  - Added `LICENSE` (MIT) and `LICENSES.md` (third-party dependencies).
  - Added `scripts/generate_licenses.py` and `LICENSES-THIRD-PARTY.md` (auto-generated Python license inventory).
  - Makefile target `licenses` / PyCharm External Tool for license regeneration.
  - Restructured `README.md` for clarity and modular references.
  - Added repository hygiene:
    - `.gitignore` for data/outputs and notebook checkpoints
    - `.gitattributes` for readable notebook diffs
    - (Optional) pre-commit with `nbstripout` to clear notebook outputs
    - PR & Issue templates under `.github/`
    - `CITATION.cff` for academic citation

## [2024–2025] – Initial Release
- Base labs (Databases, Alignment, NGS, Phylogenetics, Clustering, GCNs, ML).
