# Week 2 — Assignment (Sequence Alignment)

## General Instructions
- Use **ONLY** your own sequences downloaded in Lab 1 from NCBI (stored locally in `data/work/<handle>/lab01/`).
- You may work in Jupyter (notebook) or `.py` files.
- **Submission (Moodle)**: upload a ZIP file named `lab02_alignment_<handle>.zip` containing:
  - your scripts / notebooks;
  - a `README.txt` file (max 10 lines) with **execution steps**, Python version, and dependencies;
  - a `notes.pdf` file (max 1 page) with the requested **answers / interpretations** below.

---

## Task 1 — Pairwise distances (3p)
- Implement and run **Hamming** (only for pairs of **equal length**) **or** **p-distance** (the proportion of differing positions) for all pairs from a subset of **≥3 sequences** from your file.
- Produce a **distance matrix** (the upper triangle is sufficient).
- In `notes.pdf` (2–3 lines): specify **which two sequences are the closest** and **why** (biological / logical argument).

> **Hint:** if the lengths differ, **do not** use Hamming. Choose:  
> (a) p-distance on raw pairwise alignments (`globalxx` in Biopython)  
> or  
> (b) truncate to the minimum length and **justify your choice**.

---

## Task 2 — Pairwise alignments (4p)
- Choose two sequences from the dataset (specify their IDs).
- Run two pairwise alignments using Biopython:
  - global (e.g. `pairwise2.align.globalxx` or a variant with match/mismatch/gap scoring),
  - local (e.g. `pairwise2.align.localxx`).
- In `notes.pdf` (max 6–7 lines):
  - Compare **global vs. local** (aligned regions, number/position of gaps, score).
  - Include a **small fragment** from the alignment where the local alignment finds a match that the global one “forces” with gaps (if you identify such a case).

---

## Task 3 — Online MSA (3p)
- Choose ≥3 sequences from the dataset.
- Run a multiple sequence alignment (MSA) with **Clustal Omega (EBI)** or another equivalent online tool.
- Export the MSA result (text) and include a **relevant excerpt** in `notes.pdf` (or a permanent link, if available).
- In `notes.pdf` (max 6–7 lines):
  - Mark a **conserved region** (motif / identical segment) and explain why you think it is conserved.
  - Comment on when **MSA helps** interpretation compared to pairwise alignments.

---

## Bonus — Semiglobal (+1p)
- Run a **semiglobal alignment** (a simple implementation or a tool setting that does not penalize gaps at the ends).
- In `notes.pdf` (3–4 lines): explain **when** you would prefer semiglobal over global/local.

---

## Grading & Criteria
- **Task 1:** 3p — correctness of calculations + clear matrix  
- **Task 2:** 4p — correct execution + **interpretation** of global vs. local with example fragment  
- **Task 3:** 3p — correct MSA + **identification of conserved motif** + comparison with pairwise  
- **Bonus:** +1p — well-motivated semiglobal scenario  

---

## Academic Integrity
- All work may be done individually or in pairs, according to [docs/policies.md](../../docs/policies.md).  
- If working in pairs, list both names in `notes.pdf` and include both authors in the `.zip` archive.  
- If you use external resources (including AI), briefly note the source in `README.txt`.
