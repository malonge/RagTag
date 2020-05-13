#!/usr/bin/env python

import argparse

import pysam

from ragoo2_utilities.AGPFile import AGPFile


def main():
    parser = argparse.ArgumentParser(description="Build scaffolds from an 'orderings.bed' file")
    parser.add_argument("agp", metavar="<ragoo2.correction.agp>", type=str, help="AGP v2.1 file produced by 'ragoo2.py correct'")
    parser.add_argument("query", metavar="<query.fasta>", type=str, help="query fasta file to be scaffolded. must not be gzipped")

    args = parser.parse_args()
    agp_file = args.agp
    query_file = args.query

    x = pysam.FastaFile(query_file)
    agp = AGPFile(agp_file)

    # Iterate through the agp file
    for line in agp.iterate_lines():
        if line.is_gap:
            raise ValueError("The AGP file should have no gaps.")
        if line.orientation == "-":
            raise ValueError("No sequences should have a '-' orientation.")
        start, end = int(line.obj_beg) - 1, int(line.obj_end)
        print(">" + line.comp)
        print(x.fetch(line.obj, start, end))


if __name__ == "__main__":
    main()
