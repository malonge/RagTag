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
import re
import sys
import argparse

import pysam

from ragtag_utilities.AGPFile import AGPFile


def lens_from_fa(fname):
    fasta_file = pysam.FastaFile(fname)
    l = fasta_file.lengths
    fasta_file.close()
    return l


def lens_from_fai(fname):
    """ For any file that lists headers and lengths as the first two tab delimited columns. """
    with open(fname) as f:
        l = [int(i.split("\t")[1]) for i in f.read().rstrip().split("\n")]
    return l


def lens_from_gfa(fname):
    re_ln = re.compile(r'LN:i:\d+')
    l = []

    # Iterate through the GFA looking for sequence lines
    with open(fname) as f:
        for line in f:
            if line.startswith("S"):
                re_match = re.search(re_ln, line)
                if re_match is not None:
                    l.append(int(line[re_match.start():re_match.end()][5:]))
                else:
                    fields = line.rstrip().split("\t")
                    if fields[2] == "*":
                        raise ValueError("GFA file lacks sequence length information: {}".format(fname))
                    l.append(len(fields[2]))
    return l


def lens_from_agp(fname):
    agp_file = AGPFile(fname, mode="r")
    return [obj.obj_len for obj in agp_file.iterate_objs()]


def print_header(NG=False):
    header = ["n", "sum", "genome_size", "min", "max", "auN"]
    if NG:
        header = header + ["NG50", "NG90", "LG50", "LG90"]
    else:
        header = header + ["N50", "N90", "L50", "L90"]
    header.append("file")
    print("\t".join(header))


def print_stats(fname, seq_lens, genome_size=0):
    seq_lens = sorted(seq_lens, reverse=True)

    # Add the initial stats
    stats = [
        len(seq_lens),  # n
        sum(seq_lens)  # sum
    ]

    if not genome_size:
        stats.append("NA")
    else:
        stats.append(genome_size)

    stats.append(min(seq_lens))
    stats.append(max(seq_lens))

    # Calculate {N/L}X and auN
    csum = 0
    total = sum(seq_lens)
    if genome_size:
        total = genome_size

    aun = 0
    n50 = None
    l50 = None
    n90 = None
    l90 = None
    found_n50 = False
    found_n90 = False

    for i, l in enumerate(seq_lens):
        aun += l * (l/total)
        csum += l

        # N50
        if csum > total * 0.5:
            if not found_n50:
                l50 = i + 1
                n50 = l
                found_n50 = True

        # N90
        if csum > total * 0.9:
            if not found_n90:
                l90 = i + 1
                n90 = l
                found_n90 = True

    stats.append(round(aun))  # auN
    stats.append(n50)  # N50
    stats.append(n90)  # N90
    stats.append(l50)  # L50
    stats.append(l90)  # L90
    stats.append(fname)  # file

    print("\t".join([str(i) for i in stats]))


def main():
    parser = argparse.ArgumentParser(description="Calculate assembly sequence stats.", usage="ragtag.py asmstats [options] <in.fa>|<in.fa.fai>|<in.fa.sizes>|<in.gfa>|<in.agp> [...]")
    parser.add_argument("asm", metavar="<in.fa>|<in.fa.fai>|<in.fa.sizes>|<in.gfa>|<in.agp> [...]", nargs='*', default=[], type=str, help="assembly files (bgzipped or uncompressed for FASTA)")
    parser.add_argument("-g", metavar="INT", default=0, type=int, help="genome size for NGx [null]")

    args = parser.parse_args()
    asm_file_list = [os.path.abspath(i) for i in args.asm]
    if not asm_file_list:
        parser.print_help()
        print("\n** At least one assembly file is required **")
        sys.exit()

    print_header(NG=args.g)

    for asm_file in asm_file_list:
        if not os.path.isfile(asm_file):
            raise FileNotFoundError("File not found: {}".format(asm_file))

        if any([
            asm_file.endswith(".fasta"),
            asm_file.endswith(".fasta.gz"),
            asm_file.endswith(".fa"),
            asm_file.endswith(".fa.gz"),
            asm_file.endswith(".fna"),
            asm_file.endswith(".fna.gz")
        ]):
            seq_lens = lens_from_fa(asm_file)

        elif any([
            asm_file.endswith(".fai"),
            asm_file.endswith(".sizes"),
            asm_file.endswith(".lens")
        ]):
            seq_lens = lens_from_fai(asm_file)

        elif asm_file.endswith(".gfa"):
            seq_lens = lens_from_gfa(asm_file)
        elif asm_file.endswith(".agp"):
            seq_lens = lens_from_agp(asm_file)
        else:
            raise ValueError("Unrecognized file format: {}".format(asm_file))

        if not seq_lens:
            raise ValueError("Error computing sequence lengths: {}".format(asm_file))

        print_stats(asm_file, seq_lens, genome_size=args.g)


if __name__ == "__main__":
    main()
