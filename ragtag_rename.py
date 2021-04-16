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
    parser = argparse.ArgumentParser(description="Rename FASTA records.", usage="ragtag_rename.py <seqs.fa> [-p PREFIX]")
    parser.add_argument("fasta_fn", metavar="<seqs.fa>", default="", type=str, help="FASTA file (uncompressed or bgzipped)")
    parser.add_argument("-p", metavar="STR", type=str, default="", help="prefix")
    parser.add_argument("-o", metavar="PATH", type=str, default="ragtag.rename.agp", help="output AGP file path [./ragtag.rename.agp]")

    args = parser.parse_args()
    fasta_fn = args.fasta_fn
    prefix = args.p
    agp_fn = args.o

    agp = AGPFile(agp_fn, "w")
    record_idx = 0
    fai = pysam.FastaFile(fasta_fn)
    for reference in fai.references:
        agp.add_seq_line(
            prefix + "{0:08}".format(record_idx),
            1,
            fai.get_reference_length(reference),
            "1",
            "W",
            reference,
            1,
            fai.get_reference_length(reference),
            "+"
        )
        print(">" + prefix + "{0:08}".format(record_idx))
        print(fai.fetch(reference))

        record_idx += 1

    agp.write()
    fai.close()


if __name__ == "__main__":
    main()