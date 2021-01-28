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

import os
import sys
import argparse
from collections import defaultdict

from intervaltree import IntervalTree

from ragtag_utilities.utilities import log, get_ragtag_version
from ragtag_utilities.AGPFile import AGPFile


def sub_update(gff_file, agp_file):
    # Make a dictionary associating each original sequence with an interval tree of component sequences
    trans = defaultdict(IntervalTree)
    agp = AGPFile(agp_file, mode="r")
    for agp_line in agp.iterate_lines():

        # Check that the agp file looks correct for this task
        if agp_line.orientation == "-":
            raise ValueError("The placement BED file is not formatted correctly. No sequences should be reverse complemented for misassembly correction.")
        if not agp_line.comp_type == "W":
            raise ValueError("The placement BED file is not formatted correctly. All lines should be WGS contig (W).")
        if agp_line.is_gap:
            raise ValueError("There should be no gaps in the correction AGP file.")

        start, end = agp_line.obj_beg - 1, agp_line.obj_end
        trans[agp_line.obj][start:end] = agp_line.comp

    # Iterate through the gff intervals and update them according to trans
    with open(gff_file, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("#"):
                print(line)  # Print this comment line
            else:
                fields = line.split("\t")
                h, s, e = fields[0], int(fields[3]), int(fields[4])
                s -= 1  # Keep everything zero-indexed

                if h not in trans:
                    raise ValueError("Inconsistent input files.")

                ovlps = trans[h][s:e]
                if len(ovlps) > 1:
                    raise ValueError(
                        "%s:%d-%d in the gff file overlaps two sub sequences in the placement file. Make sure to run 'ragtag.py correct' with '--gff'" % (h, s, e)
                    )
                if len(ovlps) < 1:
                    raise ValueError("The placement BED file is not formatted correctly.")

                # Get the data from the overlapping interval and print the new line
                o = list(ovlps)[0]
                new_s = s - o.begin
                new_e = e - o.begin
                fields[0] = o.data
                fields[3] = str(new_s + 1)  # back to one-based indexing for gff format
                fields[4] = str(new_e)
                print("\t".join(fields))


def sup_update(gff_file, agp_file):
    # Make a dictionary associating each original sequence with the destination sequence
    trans = {}
    strands = {}
    seq_lens = {}
    agp = AGPFile(agp_file, mode="r")
    for agp_line in agp.iterate_lines():
        if not agp_line.is_gap:
            start, end = agp_line.obj_beg - 1, agp_line.obj_end
            trans[agp_line.comp] = (start, end, agp_line.obj)
            strands[agp_line.comp] = agp_line.orientation
            seq_lens[agp_line.comp] = end - start

    # Iterate through the gff intervals and update them according to trans
    with open(gff_file, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("#"):
                print(line)  # Print this comment line
            else:
                fields = line.split("\t")
                h, s, e, st = fields[0], int(fields[3]), int(fields[4]), fields[6]
                s -= 1  # Keep everything zero-indexed

                if h not in trans:
                    print()
                    print(line)
                    raise ValueError("Inconsistent input files.")

                # Check if the original sequence has been reverse complemented
                if strands[h] == "-":
                    l = seq_lens[h]
                    s, e = l-e, l-s
                    if st == "+":
                        st = "-"
                    else:
                        st = "+"

                new_s = trans[h][0] + s
                new_e = trans[h][0] + e
                fields[0] = trans[h][2]
                fields[3] = str(new_s + 1)  # back to one-based indexing for gff format
                fields[4] = str(new_e)
                fields[6] = st
                print("\t".join(fields))


def main():
    parser = argparse.ArgumentParser(description="Update gff intervals given a RagTag AGP file", usage="ragtag.py updategff [-c] <genes.gff> <ragtag.agp>")
    parser.add_argument("gff", nargs='?', default="", metavar="<genes.gff>", type=str, help="gff file")
    parser.add_argument("agp", nargs='?', default="", metavar="<ragtag.*.agp>", type=str, help="agp file")
    parser.add_argument("-c", action="store_true", default=False, help="update for misassembly correction (ragtag.correction.agp)")

    args = parser.parse_args()

    if not args.gff or not args.agp:
        parser.print_help()
        sys.exit()

    log("RagTag " + get_ragtag_version())
    log("CMD: ragtag.py updategff " + " ".join(sys.argv[1:]))

    gff_file = os.path.abspath(args.gff)
    agp_file = os.path.abspath(args.agp)
    is_sub = args.c

    if is_sub:
        sub_update(gff_file, agp_file)
    else:
        sup_update(gff_file, agp_file)

    log("Goodbye")


if __name__ == "__main__":
    main()
