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
import gzip
import argparse
from collections import defaultdict


class Alignment:
    def __init__(self, in_r_start, in_r_end, in_q_start, in_q_end, in_cigar, in_strand, in_num_matches, in_aln_len, in_num_mismatches):
        self.r_start = int(in_r_start) + 1
        self.r_end = int(in_r_end)
        self.q_start = int(in_q_start) + 1
        self.q_end = int(in_q_end)
        self.cigar = in_cigar
        self.strand = in_strand
        self.num_matches = int(in_num_matches)
        self.aln_len = int(in_aln_len)
        self.num_mismatches = in_num_mismatches

        self.parsed_cigar = []
        self._parse_cigar()

        if self.strand == "-":
            self.q_start, self.q_end = self.q_end, self.q_start

    def _parse_cigar(self):
        """ Given a CIGAR string, create a list of CIGAR operations. """
        re_cg = re.compile(r'(\d+)([MIDNSHP=X])')
        self.parsed_cigar = [(int(i[0]), i[1]) for i in list(re.findall(re_cg, self.cigar))]


def write_delta(in_alns, in_r_lens, in_q_lens, ref_file, query_file):
    if not ref_file:
        ref_file = "file1"
    if not query_file:
        query_file = "file2"
    print(ref_file + " " + query_file)
    print('NUCMER')

    # Iterate over each reference-query header pair
    for r_header in in_alns.keys():
        for q_header in in_alns[r_header].keys():
            print('>%s %s %d %d' % (r_header, q_header, in_r_lens[r_header], in_q_lens[q_header]))
            for z in in_alns[r_header][q_header]:
                print('%d %d %d %d %d %d %d' % (
                    z.r_start,
                    z.r_end,
                    z.q_start,
                    z.q_end,
                    z.num_mismatches,
                    z.num_mismatches,
                    0
                ))
                # Continue with the cigar string
                offsets = []
                cigar = z.parsed_cigar
                if cigar[0][1] == 'S' or cigar[0][1] == 'H':
                    cigar = cigar[1:-1]
                else:
                    cigar = cigar[:-1]

                counter = 1
                for op in cigar:
                    if op[1] == "M":
                        counter += op[0]
                    elif op[1] == "D":
                        offsets.append(counter)
                        num_I = op[0]
                        for i in range(1, num_I):
                            offsets.append(1)
                        counter = 1
                    elif op[1] == 'I':
                        offsets.append(-1 * counter)
                        num_I = op[0]
                        for i in range(1, num_I):
                            offsets.append(-1)
                        counter = 1
                    else:
                        raise ValueError('Unexpected CIGAR code')
                offsets.append(0)
                offsets = [str(a) for a in offsets]
                print('\n'.join(offsets))


def paf2delta():
    parser = argparse.ArgumentParser(description="Convert a PAF file to a Nucmer delta file.\nPAF file lines must have CIGAR strings ('-c' when using Minimap2)", usage="ragtag.py paf2delta <with-cg.paf>")
    parser.add_argument("paf_file", metavar="<with-cg.paf>", type=str, help="PAF file to convert (gzip allowed).")
    parser.add_argument("-r", default="", metavar="PATH", type=str, help="PATH to there reference fasta file")
    parser.add_argument("-q", default="", metavar="PATH", type=str, help="PATH to there query fasta file")

    args = parser.parse_args()
    paf_file = args.paf_file
    ref_file = args.r
    query_file = args.q
    alns = dict()

    # Dictionaries to store reference and query sequence lengths
    r_chr_lens = dict()
    q_chr_lens = dict()

    if paf_file[-3:] == ".gz":
        f = gzip.open(paf_file)
    else:
        f = open(paf_file, 'r')

    for line in f:
        if not isinstance(line, str):
            line = line.decode("utf-8")

        fields = line.split('\t')

        # Get the reference/query sequence lengths
        r_header = fields[5]
        q_header = fields[0]
        if r_header not in r_chr_lens:
            r_chr_lens[r_header] = int(fields[6])
        if q_header not in q_chr_lens:
            q_chr_lens[q_header] = int(fields[1])

        # Get the rest of the info and instantiate the Alignment object
        cs = None
        nm = None
        for i in fields[12:]:
            if i.startswith("cg:Z:"):
                cs = i[5:]

            if i.startswith("NM:i:"):
                nm = int(i[5:])

        if cs is None:
            raise ValueError("PAF file must contain a CIGAR string. Use 'minimap2 -c'")

        if nm is None:
            raise ValueError('PAF file must include NM tag.')

        x = Alignment(
            fields[7],
            fields[8],
            fields[2],
            fields[3],
            cs,
            fields[4],
            fields[9],
            fields[10],
            nm
        )

        # Add the alignments to the nested dictionary (first key=reference header, second key=query header)
        if r_header not in alns:
            alns[r_header] = defaultdict(list)
        alns[r_header][q_header].append(x)
    f.close()

    write_delta(alns, r_chr_lens, q_chr_lens, ref_file, query_file)


if __name__ == "__main__":
    paf2delta()