#!/usr/bin/env python

"""
MIT License

Copyright (c) 2020 Michael Alonge <malonge11@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import sys
import argparse

import pysam

from ragtag_utilities.utilities import reverse_complement
from ragtag_utilities.AGPFile import AGPFile


def main():
    parser = argparse.ArgumentParser(description="Build sequences in FASTA format from an AGP v2.1 file.", usage="ragtag.py agp2fasta <scaffolds.agp> <components.fasta>")
    parser.add_argument("agp", metavar="<scaffolds.agp>", nargs='?', default="", type=str, help="AGP v2.1 file")
    parser.add_argument("components", metavar="<components.fasta>", nargs='?', default="", type=str, help="component FASTA file (can be uncompressed or bgzipped)")

    args = parser.parse_args()
    if not args.agp or not args.components:
        parser.print_help()
        sys.exit()

    agp_file = args.agp
    components_file = args.components

    fai = pysam.FastaFile(components_file)
    agp = AGPFile(agp_file, mode="r")

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
                sys.stdout.write(reverse_complement(fai.fetch(agp_line.comp, agp_line.comp_beg-1, agp_line.comp_end)))
            else:
                sys.stdout.write(fai.fetch(agp_line.comp, agp_line.comp_beg-1, agp_line.comp_end))

    # End the FASTA file with a newline
    sys.stdout.write("\n")
    fai.close()


if __name__ == "__main__":
    main()
