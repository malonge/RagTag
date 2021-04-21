#!/usr/bin/env python

"""
MIT License

Copyright (c) 2021 Michael Alonge <malonge11@gmail.com>

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

import re
import sys
import argparse

import pysam

from ragtag_utilities.utilities import get_ragtag_version
from ragtag_utilities.AGPFile import AGPFile


def main():
    parser = argparse.ArgumentParser(description='Split sequencs at gaps', usage="ragtag.py splitasm <asm.fa>")
    parser.add_argument("asm", metavar="<asm.fa>", default="", type=str, help="assembly fasta file (uncompressed or bgzipped)")
    parser.add_argument("-n", metavar="INT", type=int, default=0, help="minimum gap size [0]")
    parser.add_argument("-o", metavar="PATH", type=str, default="ragtag.splitasm.agp", help="output AGP file path [./ragtag.splitasm.agp]")

    # Parse the command line arguments
    args = parser.parse_args()
    if not args.asm:
        parser.print_help()
        print("\n** The assembly FASTA file is required **")
        sys.exit()

    asm_fn = args.asm
    min_gap_size = args.n
    agp_fn = args.o

    # Initialize the AGP file
    agp = AGPFile(agp_fn, mode="w")
    agp.add_pragma()
    agp.add_comment("# AGP created by RagTag {}".format(get_ragtag_version()))

    # Process the FASTA file
    new_header_idx = 0
    fai = pysam.FastaFile(asm_fn)
    for header in sorted(fai.references):
        seq = fai.fetch(header).upper()
        seq_len = fai.get_reference_length(header)
        gap_coords = [(i.start(), i.end()) for i in re.finditer(r'N+', seq) if i.end() - i.start() > min_gap_size]

        if not gap_coords:
            new_header = "seq{0:08}".format(new_header_idx)
            new_header_idx += 1
            agp.add_seq_line(header, "1", seq_len, "1", "W", new_header, "1", seq_len, "+")
        else:
            gap_coords.append((seq_len, seq_len + 1))
            pid = 1
            if gap_coords[0][0]:
                # The sequence doesn't start with a gap
                new_header = "seq{0:08}".format(new_header_idx)
                agp.add_seq_line(header, "1", str(gap_coords[0][0]), str(pid), "W", new_header, "1", str(gap_coords[0][0]), "+")
                new_header_idx += 1
                pid += 1

            for i in range(1, len(gap_coords)):
                # Add the gap line
                gap_start, gap_end = gap_coords[i-1][0], gap_coords[i-1][1]
                gap_len = gap_end - gap_start
                agp.add_gap_line(header, str(gap_start + 1), str(gap_end), str(pid), "N", str(gap_len), "scaffold", "yes", "align_genus")
                pid += 1

                # Add the sequence line
                obj_start, obj_end = gap_coords[i-1][1], gap_coords[i][0]
                comp_len = obj_end - obj_start
                new_header = "seq{0:08}".format(new_header_idx)
                if gap_coords[i-1][1] != seq_len:
                    agp.add_seq_line(header, str(obj_start + 1), obj_end, pid, "W", new_header, "1", str(comp_len), "+")
                    new_header_idx += 1
                    pid += 1

        agp.write()

    # Iterate over the AGP file and print the sequences
    agp = AGPFile(agp_fn, mode="r")
    for line in agp.iterate_lines():
        if not line.is_gap:
            obj, comp, obj_beg, obj_end = line.obj, line.comp, line.obj_beg, line.obj_end
            print(">" + comp)
            print(fai.fetch(obj, obj_beg-1, obj_end))

    fai.close()


if __name__ == "__main__":
    main()
