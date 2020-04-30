#!/usr/bin/env python

import argparse
from collections import OrderedDict

import pysam

from ragoo2_utilities.utilities import reverse_complement


def main():
    parser = argparse.ArgumentParser(description="Build scaffolds from an 'orderings.bed' file")
    parser.add_argument("orderings", metavar="<orderings.bed>", type=str, help="'orderings.bed' file produced by 'ragoo2_scaffold.py'")
    parser.add_argument("query", metavar="<query.fasta>", type=str, help="query fasta file to be scaffolded. must not be gzipped")
    parser.add_argument("out_fasta_file", metavar="<ragoo.fasta>", type=str, help="output fasta file name")

    args = parser.parse_args()
    orderings_file = args.orderings
    query_file = args.query
    out_file = args.out_fasta_file

    # Organize the orderings
    orderings = OrderedDict()
    with open(orderings_file, "r") as f:
        for line in f:
            l = line.rstrip().split()
            ref, start, end, seq_type, header, strand = l[0], int(l[1]), int(l[2]), l[3], l[4], l[5]
            if ref not in orderings:
                orderings[ref] = []
            orderings[ref].append((start, end, seq_type, header, strand))

    # Sort the orderings and write the fasta output
    out_fasta = open(out_file, "w")
    x = pysam.FastaFile(query_file)
    placed_q_seqs = set()
    for i in orderings:
        out_fasta.write(">" + i + "\n")
        orderings[i] = sorted(orderings[i])

        # Iterate through each sequence of this scaffold
        for j in orderings[i]:
            q_header = j[3]
            if j[2] == "s":
                # This is a query sequence
                placed_q_seqs.add(q_header)
                if j[4] == "+":
                    q_seq = x.fetch(q_header)
                    if not len(q_seq) == j[1] - j[0]:
                        raise RuntimeError("Inconsistency between query fasta and placement file.")
                    out_fasta.write(q_seq)
                else:
                    out_fasta.write(reverse_complement(x.fetch(q_header)))
            else:
                # This is a gap
                assert j[2] == "g"
                gap_len = j[1] - j[0]
                out_fasta.write("N" * gap_len)

        out_fasta.write("\n")
    out_fasta.close()


if __name__ == "__main__":
    main()
