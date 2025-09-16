# BDHB — Big Data in Health & Bioinformatics Labs

[![CI](https://github.com/bozdogalex/BDHB-lab/actions/workflows/ci.yml/badge.svg)](https://github.com/bozdogalex/BDHB-lab/actions/workflows/ci.yml)
[![Open in GitHub Codespaces](https://img.shields.io/badge/Open%20in-Codespaces-000?logo=github)](https://github.com/codespaces/new?hide_repo_select=true&ref=main&repo=bozdogalex%2FBDHB-lab)

> Master’s-level laboratories blending classical bioinformatics with modern ML, networks, and GenAI. The environment is CPU‑only and identical across Codespaces and Docker via the prebuilt image `ghcr.io/bozdogalex/bdhb:base`.


## Labs (index)

- 00 — Smoke: `labs/00_smoke/`
- 01 — Databases & GitHub: `labs/01_databases/`
- 02 — Sequence Alignment: `labs/02_alignment/`
- 03 — NGS: `labs/03_ngs/`
- 04 — Phylogenetics: `labs/04_phylogenetics/`
- 05 — Clustering: `labs/05_clustering/`
- 06a — WGCNA (+ Diseasome): `labs/06a_wgcna/`
- 06b — Network Viz & GNN: `labs/06b_network_viz/`
- 07 — Federated Learning: `labs/07_ml_flower/`
- 08 — Drug Repurposing: `labs/08_repurposing/`
- 09 — Integrative + Digital Twin: `labs/09_integrative/`
- 10 — Multi‑omics (+ optional Quantum): `labs/10_multiomics/`
- 11 — Generative AI (PubMed vs regex; protein embeddings): `labs/11_genai/`

---


> Full onboarding (screenshots, tips): **docs/onboarding.md**

---

## Repo map

- `labs/` — all weekly lab content
- `docs/` — onboarding, ANIS pack (before/after, one‑pagers, screenshots)
- `mlops/` — MLflow helpers
- `.devcontainer/` — Codespaces/Devcontainer (pulls GHCR image)
- `.github/workflows/` — CI + image publish
- `Dockerfile`, `requirements.txt` — env definition
- `dev.ps1`, `Makefile` — local helpers

---

## Contributing / Policies / Citation

- `CONTRIBUTING.md` (root) — contribution rules & PR tips  
- `CODE_OF_CONDUCT.md` (optional) — community standards  
- `CITATION.cff` (root) — how to cite this work  
- `docs/changelog.md` — changelog (linked from releases)

