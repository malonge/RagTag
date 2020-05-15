#!/usr/bin/env python

import sys
import argparse

import pysam

from ragoo2_utilities.utilities import reverse_complement
from ragoo2_utilities.AGPFile import AGPFile


def main():
    parser = argparse.ArgumentParser(description="Build sequences in FASTA format from an AGP v2.1 file. ")
    parser.add_argument("agp", metavar="<scaffolds.agp>", type=str, help="AGP v2.1 file")
    parser.add_argument("components", metavar="<components.fasta>", type=str, help="FASTA file with component sequences to be scaffolded. must not be gzipped")

    args = parser.parse_args()
    agp_file = args.agp
    components_file = args.components

    fai = pysam.FastaFile(components_file)
    agp = AGPFile(agp_file)

    # Iterate over the lines of the AGP file
    prev_obj = None
    is_first = True
    for agp_line in agp.iterate_lines():
        if agp_line.obj != prev_obj:
            if is_first:
                print(">" + agp_line.obj)
                is_first = False
            else:
                print("\n>" + agp_line.obj)

            prev_obj = agp_line.obj

        if agp_line.is_gap:
            sys.stdout.write("N"*agp_line.gap_len)
        else:
            if agp_line.orientation == "-":
                sys.stdout.write(reverse_complement(fai.fetch(agp_line.comp)))
            else:
                sys.stdout.write(fai.fetch(agp_line.comp))

    # End the FASTA file with a newline
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
