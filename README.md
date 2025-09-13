# BDHB – Big Data in Health and Bioinformatics Laboratory

[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Open in GitHub Codespaces](https://img.shields.io/badge/Codespaces-Open-blue?logo=github)](https://github.com/codespaces/new?hide_repo_select=true&ref=main&repo=bozdogalex/BDHB-lab)

## 📖 Overview
The **BDHB Laboratory** is a master’s-level course at Politehnica University of Timișoara (Faculty of Automation and Computing).  
It introduces students to cutting-edge methods in **big data analytics, bioinformatics, and health informatics**, combining classical bioinformatics workflows with modern machine learning, network analysis, and generative AI.

- Format: **Hybrid** – Tuesdays, 18:00–20:00  
- Platform: Microsoft Teams (link in Moodle)  
- Repository: contains weekly labs, resources, and assignments.

---

## 📚 Lab Outline
- **Week 1**: Databases & GitHub Onboarding  
- **Week 2**: Sequence Alignment  
- **Week 3**: Next-Generation Sequencing (NGS) Analysis  
- **Week 4**: Phylogenetics  
- **Week 5**: Clustering Techniques  
- **Week 6a**: Gene Co-Expression Networks (WGCNA)  
- **Week 6b**: Network Visualization & Graph Neural Networks  
- **Week 7**: Machine Learning + Federated Learning  
- **Week 8**: Network-Based Approaches in Drug Repurposing  
- **Week 9**: Integrative Genomics + Digital Twin Narrative  
- **Week 10**: Multi-Omics Integration (+ optional Quantum demo)  
- **Week 11**: Advanced ML & Generative AI in Bioinformatics  
- **Week 12**: Project Work & Consultations  
- **Week 13**: Final Project Presentations  

---

## 🛠️ Setup

### Option A: GitHub Codespaces (recommended)
Click “Open in Codespaces” → Jupyter and dependencies install automatically.  

### Option B: Local with Docker
```bash
docker build -t bdhb:base .
docker run -it -p 8888:8888 -v $PWD:/work bdhb:base
```

### Option C: Manual local install
pip install -r requirements.txt

## 📂 Resources
See [RESOURCES.md](RESOURCES.md) for the full list of textbooks, articles, and online tools.  

---

## 🤝 Contributing
- Pair work is highly suggested (roles A/B).  
- Contributions follow the rules in [CONTRIBUTING.md](CONTRIBUTING.md).  

---

## ⚖️ Policies
Governance, licensing, GDPR, and sustainability are described in [POLICIES.md](POLICIES.md).  

---

## 📈 Continuity
- Content updated annually (>60% refreshed labs).  
- Best student projects are merged into this repo.  
- Integration with **Oncohelp clinical collaborations** → potential dissertations and research publications.  

