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

import argparse

import pysam

from ragtag_utilities.AGPFile import AGPFile


def main():
    parser = argparse.ArgumentParser(description="Break corrected query sequences (objects) into components.")
    parser.add_argument("agp", metavar="<ragtag.correction.agp>", type=str, help="AGP v2.1 file produced by 'ragtag.py correct'")
    parser.add_argument("query", metavar="<query.fasta>", type=str, help="query fasta file corresponding to objects in <ragtag.correction.agp> (can be uncompressed or bgzipped")

    args = parser.parse_args()
    agp_file = args.agp
    query_file = args.query

    fai = pysam.FastaFile(query_file)
    agp = AGPFile(agp_file, mode="r")

    # Iterate through the agp file
    for line in agp.iterate_lines():
        if line.is_gap:
            raise ValueError("The AGP file should have no gaps.")
        if line.orientation == "-":
            raise ValueError("No sequences should have a '-' orientation.")
        start, end = int(line.obj_beg) - 1, int(line.obj_end)
        print(">" + line.comp)
        print(fai.fetch(line.obj, start, end))

    fai.close()


if __name__ == "__main__":
    main()
