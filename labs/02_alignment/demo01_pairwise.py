#!/usr/bin/env python
"""
Demo: local and global alignment with Biopython (pairwise).
- Reuse data from data/sample/; extract short subsequences for debugging.
Run:
  python labs/02_alignment/demo01_pairwise_biopython.py --fasta data/sample/tp53_dna_multi.fasta
"""
import argparse
from Bio import SeqIO, pairwise2


def take_two_short_subseqs(fasta_path, k=7):
    """
    Extract the first k nucleotides from the first two sequences in the FASTA file.
    Used to create small test cases for debugging.
    """
    recs = [r for r in SeqIO.parse(fasta_path, "fasta")]
    if len(recs) < 2:
        raise ValueError("Need at least 2 sequences in the FASTA file.")
    a = str(recs[0].seq)[:k]
    b = str(recs[1].seq)[:k]
    return a, b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="path to the FASTA file")
    ap.add_argument("--k", type=int, default=7, help="substring length")
    args = ap.parse_args()

    A, B = take_two_short_subseqs(args.fasta, k=args.k)

    # Simple scoring scheme: +1 for match, -1 for mismatch, -1 for gap
    global_alignments = pairwise2.align.globalms(A, B, 1, -1, -1, -1)
    local_alignments = pairwise2.align.localms(A, B, 1, -1, -1, -1)

    print("[INPUT]")
    print("A:", A)
    print("B:", B)

    print("\n[GLOBAL] top alignment:")
    print(pairwise2.format_alignment(*global_alignments[0]))

    print("[LOCAL] top alignment:")
    print(pairwise2.format_alignment(*local_alignments[0]))


if __name__ == "__main__":
    main()
