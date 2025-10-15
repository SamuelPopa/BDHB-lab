#!/usr/bin/env python
"""
Computes p-distance / Hamming distance for all pairs in a FASTA file.
- Reuses data/sample/* ; for sequences of different lengths, truncates to min(len).
Run:
  python labs/02_alignment/demo02_distance_matrix.py --fasta data/sample/tp53_dna_multi.fasta
"""
import argparse
from itertools import combinations
from Bio import SeqIO


def hamming_equal(a, b):
    """Compute the Hamming distance for two sequences of equal length."""
    return sum(x != y for x, y in zip(a, b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="path to the FASTA file")
    args = ap.parse_args()

    recs = list(SeqIO.parse(args.fasta, "fasta"))
    ids = [r.id for r in recs]
    seqs = [str(r.seq) for r in recs]

    print("pair,hamming,p_distance,len_used")
    for (i, j) in combinations(range(len(seqs)), 2):
        a, b = seqs[i], seqs[j]
        L = min(len(a), len(b))
        a2, b2 = a[:L], b[:L]  # truncate for simple comparison
        d_h = hamming_equal(a2, b2)
        d_p = d_h / float(L) if L > 0 else 0.0
        print(f"{ids[i]}-{ids[j]},{d_h},{d_p:.4f},{L}")


if __name__ == "__main__":
    main()
