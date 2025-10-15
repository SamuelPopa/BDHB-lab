# Week 1 — Databases & GitHub

## Objectives:
1) All students run the same environment (Codespaces or Docker + Jupyter).  
2) All students go once through the GitHub workflow (fork → branch → PR) with a short task.

> In weeks 2–12 we will **reuse the data** collected early on (GEO/TCGA/NCBI/Ensembl) for alignment, NGS, phylogeny, co-expression, ML, etc. (see the course calendar).

---

## Part 0 — Environment check

Choose **one** option:

### A) Codespaces
1. Open a Codespace on branch `main`.
2. In the terminal, run:
   ```bash
   python labs/00_smoke/smoke.py
   ``` 
Expected output: ok.

### B) Local Docker (Windows PowerShell)
```docker pull ghcr.io/bozdogalex/BDHB-lab:base
docker run -it --rm -p 8890:8888 -v "${PWD}:/work" -w /work ghcr.io/bozdogalex/BDHB-lab:base `
  bash -lc "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --IdentityProvider.token='' --allow-root"
  ```

- Open `http://localhost:8890/lab`  and run print("ok").
- Note (paths with & on Windows): when running scripts from this folder, put the path in quotes:
`python "labs/01_intro&databases/demo_entrez_brca1.py"`

## Part 1 — PR Task

- Read git-workflow.md
- Add your GitHub handle to labs/01_intro&databases/roster/handles.csv (format: FirstName LastName,github_handle) via a Pull Request.
- Steps: 
  - Fork the repo → create branch feat/roster-<handle>.
  - Edit labs/01_intro&databases/roster/handles.csv → add one line.
  - Run git commit -m "Add <handle> to roster" → git push.
  - Complete the PR checklist (template: week1_roster.md).

## Part 2 — Demo / Exercises
**Run**
`demo01_entrez_brca1.py` — search + download BRCA1 (GenBank) and compute GC summary.
`demo02_seq_ops.py` — basic sequence operations (transcription, translation, GC, motif).
`demo03_dbsnp.py` — quick dbSNP query and summary.
**Complete and run**
- `ex01_multifasta_gc.py` — fill in the TODOs for FASTA download (Entrez) and run GC calculation on the DOWNLOADED file.
- Example (query):
`python labs/01_intro&databases/ex01_multifasta_gc.py --email student@example.com --query "TP53[Gene] AND Homo sapiens[Organism]" --retmax 3 --out data/work/<handle>/lab01/my_tp53.fa`
- Example (accession):
`python labs/01_intro&databases/ex01_multifasta_gc.py --email student@example.com --accession NM_000546 --out data/work/<handle>/lab01/nm000546.fa`

- Note: save your own files in data/work/<handle>/lab01/. This folder is ignored by git — do NOT upload data to the repository.
**Deliverables**
- Your Pull Request must include:
  - A new line in labs/01_intro&databases/roster/handles.csv.
  - A file labs/01_intro&databases/submissions/<github_handle>/notes.md containing:
  - Confirmation that you have run all demo scripts and the exercise.
  - A simple result observed (e.g., “GC fraction = 0.47”).
  - The completed exercise saved in:
`labs/01_intro&databases/submissions/<github_handle>/ex01_multifasta_gc.py`
  - The completed PR checklist (template file).

---

## Next Week
- We will use the downloaded FASTA files to perform sequence alignments (global and local, NW/SW).
- See [Week 2 — Sequence Alignment](../02_alignment/README.md)

---

## Skills:

- Running a reproducible environment (Codespaces/Docker).
- Correctly opening and completing a PR (fork → branch → PR).
- Performing first queries and basic operations on biological sequences.

---

## Resources : 
- [NCBI](https://www.ncbi.nlm.nih.gov/)  
- [Ensembl Genome Browser](https://www.ensembl.org/)  
- [dbSNP (NCBI)](https://www.ncbi.nlm.nih.gov/snp/)  
- [TCGA (The Cancer Genome Atlas) Portal](https://portal.gdc.cancer.gov/)  
- Book: [Pevsner, *
- Bioinformatics and Functional Genomics*, 3rd ed., Wiley Blackwell, 2015](https://genetics.elte.hu/oktatasi_anyag/archivum/bioinfo/Bioinformatika_2018-2019/book.pdf)  
- Book: [Lesk, *Introduction to Bioinformatics*, 5th ed., Oxford University Press, 2019](https://edscl.in/pluginfile.php/3340/mod_folder/content/0/Introduction%20To%20Bioinformatics.pdf?forcedownload=1)  



