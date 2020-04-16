#!/usr/bin/env python

import argparse

import pysam


def main():
    parser = argparse.ArgumentParser(description="Build scaffolds from an 'orderings.bed' file")
    parser.add_argument("breaks", metavar="<breaks.bed>", type=str, help="'breaks.bed' file produced by 'ragoo_correct.py'")
    parser.add_argument("query", metavar="<query.fasta>", type=str, help="query fasta file to be scaffolded. must not be gzipped")
    parser.add_argument("out_fasta_file", metavar="<query.break.fasta>", type=str, help="output fasta file name")

    args = parser.parse_args()
    breaks_file = args.breaks
    query_file = args.query
    out_file = args.out_fasta_file

    x = pysam.FastaFile(query_file)
    fout = open(out_file, "w")

    # Iterate through the bed file
    with open(breaks_file, "r") as f:
        for line in f:
            header, start, end, stype, qheader, strand = line.rstrip().split()
            start, end = int(start), int(end)
            fout.write(">" + qheader + "\n")
            fout.write(x.fetch(header, start, end) + "\n")

    fout.close()


if __name__ == "__main__":
    main()
