#!/usr/bin/env python

import argparse

import pysam


def main():
    parser = argparse.ArgumentParser(description="Build scaffolds from an 'orderings.bed' file")
    parser.add_argument("agp", metavar="<ragoo2.correction.agp>", type=str, help="AGP v2.1 file produced by 'ragoo2.py correct'")
    parser.add_argument("query", metavar="<query.fasta>", type=str, help="query fasta file to be scaffolded. must not be gzipped")
    parser.add_argument("out_fasta_file", metavar="<query.break.fasta>", type=str, help="output fasta file name")

    args = parser.parse_args()
    agp_file = args.agp
    query_file = args.query
    out_file = args.out_fasta_file

    x = pysam.FastaFile(query_file)
    fout = open(out_file, "w")

    # Iterate through the agp file
    with open(agp_file, "r") as f:
        for line in f:
            obj, obj_start, obj_end, pid, ctype, comp, comp_beg, comp_end, strand = line.rstrip().split("\t")
            start, end = int(obj_start)-1, int(obj_end)
            fout.write(">" + comp + "\n")
            fout.write(x.fetch(obj, start, end) + "\n")

    fout.close()


if __name__ == "__main__":
    main()
