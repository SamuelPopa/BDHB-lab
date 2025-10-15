# Week 2 — Sequence Alignment

## Objectives
-  Understand the types of alignments: global (Needleman–Wunsch), local (Smith–Waterman), and semi-global.
- Learn to use substitution matrices (PAM, BLOSUM).
- Practice with Biopython and external tools (BLAST, Clustal Omega).
- Develop skills for implementing the basic alignment algorithms.
---

## Part 1 — Demo / Exercises
**Run**
- demo01_pairwise_biopython.py — global and local alignment using Biopython (pairwise2).
- demo02_distance_matrix.py — compute sequence distances (p-distance, Hamming) from FASTA files.
**Complete and run**
- ex02_global_nw.py — implementation skeleton for global alignment (TODO).
- ex03_local_sw.py — implementation skeleton for local alignment (TODO).
Note: Use the datasets downloaded in Lab 1 (from `data/work/<handle>/lab01/`).

---

## Deliverables

**Your Pull Request (PR) must include:*
- The file
   `labs/02_alignment/<github_handle>_notes.md` containing:
   - which datasets you used (e.g., TP53 vs. BRCA1),
   - a short reflection: “When is global alignment preferred over local alignment?”
   - The completed exercises, saved under:
- `labs/02_alignment/submissions/<github_handle>/ex02_global_nw.py`
- `labs/02_alignment/submissions/<github_handle>/ex03_local_sw.py`
- The completed checklist from the PR template.

--- 

## Next Week
- We will extend the analysis to NGS reads (FASTQ → mapping → variant calling).
- The alignments obtained this week will be used to validate mapping and NGS analyses.
- See Week 3 — NGS Analysis

---

## Learning Outcomes
- Understand the difference between global and local alignment.
- Use Biopython for simple pairwise alignments.
- Implement the basic NW and SW algorithms.
- Interpret results and compare them with BLAST / Clustal outputs.

---

## Resources 
- [Global Alignment (Needleman–Wunsch)](../../docs/presentations/alignment1.pdf)  
- [Local Alignment (Smith–Waterman)](../../docs/presentations/alignment2.pdf)  
- [Applied Bioinformatics of Nucleic Acids — Cap. 1](../../docs/papers/Applied_Bioinformatics.pdf)  
- [Scoring Matrix Development (BLOSUM62) (pdf în /papers)](../../docs/papers/Scoring_matrix_development_BLOSUM62.pdf)  
- Substitution matrices: [BLOSUM62 (NCBI)](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62)  
- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)  
- [Clustal Omega — Multiple Sequence Alignment](https://www.ebi.ac.uk/Tools/msa/clustalo/)  
- [Biopython pairwise2](https://biopython.org/docs/1.75/api/Bio.pairwise2.html)  