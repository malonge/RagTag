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

import pysam

from ragtag_utilities.RestrictionFragmentMap import RestrictionFragmentMap, RestrictionEnzymes
from ragtag_utilities.utilities import log

"""
This code is mostly inspired/influenced by SALSA2.

https://doi.org/10.1371/journal.pcbi.1007273
"""


def count_links(bam_fname, l_cutoffs, r_cutoffs):
    """
    Count the raw links between contig halves.
    :param bam_fname: Hi-C BAM file name
    :param l_cutoffs: sequence -> cutoff for the first half of the sequence
    :param r_cutoffs: sequence -> cutoff for the second half of the sequence
    :return: sequence half -> sequence half -> number of Hi-C links
    """
    links = dict()
    bam = pysam.AlignmentFile(bam_fname, "rb")

    # Check the SAM metadata
    if "HD" not in bam.header:
        raise RuntimeError("Alignment file must be sorted be read name (samtools sort -n)")

    if "SO" not in bam.header["HD"]:
        raise RuntimeError("Alignment file must be sorted be read name (samtools sort -n)")

    if bam.header["HD"]["SO"] != "queryname":
        raise RuntimeError("Alignment file must be sorted be read name (samtools sort -n) but it is currently sorted by: {}".format(bam.header["HD"]["SO"]))

    # Iterate through each alignment
    prev_r_name = ""
    prev_q_name = ""
    prev_mate = True
    prev_pos = 0

    for aln in bam.fetch(until_eof=True):
        if not aln.is_unmapped and not aln.is_secondary:
            r_name = bam.get_reference_name(aln.reference_id)
            q_name = aln.query_name
            mate = aln.is_read1
            pos = ((aln.reference_end - aln.reference_start) // 2) + aln.reference_start

            first_seq, second_seq = "", ""

            # Check if this is the second mate of a pair and it maps to a distinct reference sequence
            if q_name == prev_q_name and mate != prev_mate and r_name != prev_r_name:
                # Assign the nodes that we are incrementing
                if prev_pos <= l_cutoffs[prev_r_name]:
                    first_seq = prev_r_name + "_b"

                if prev_pos > r_cutoffs[prev_r_name]:
                    first_seq = prev_r_name + "_e"

                if pos <= l_cutoffs[r_name]:
                    second_seq = r_name + "_b"

                if pos > r_cutoffs[r_name]:
                    second_seq = r_name + "_e"

                assert first_seq != second_seq
                if first_seq and second_seq:
                    if second_seq < first_seq:
                        first_seq, second_seq = second_seq, first_seq

                    # Increment the link between these two sequences
                    if first_seq not in links:
                        links[first_seq] = dict()
                    if second_seq not in links[first_seq]:
                        links[first_seq][second_seq] = 0

                    links[first_seq][second_seq] += 1

            prev_r_name = r_name
            prev_q_name = q_name
            prev_mate = mate
            prev_pos = pos

    return links


def normalize_links(raw_links, l_norm_factors, r_norm_factors):
    """

    :param raw_links:
    :param l_norm_factors:
    :param r_norm_factors:
    :return:
    """
    norm_links = dict()

    # Iterate through reach link
    for i in raw_links:
        i_base = i[:-2]
        norm_links[i] = dict()

        for j in raw_links[i]:
            j_base = j[:-2]
            # Get the normalization factors for both contig ends
            first_norm = l_norm_factors[i_base] if i.endswith("_b") else r_norm_factors[i_base]
            second_norm = l_norm_factors[j_base] if j.endswith("_b") else r_norm_factors[j_base]

            norm_links[i][j] = raw_links[i][j] / (first_norm + second_norm)

    return norm_links


def write_links(raw_links, norm_links, out_file=sys.stdout):
    for i in raw_links:
        for j in raw_links[i]:
            out_file.write("\t".join([
                i,
                j,
                str(norm_links[i][j]),
                str(raw_links[i][j])
            ]) + "\n")


def main():
    description = """  
    """
    parser = argparse.ArgumentParser(description="Quantify links from a Hi-C BAM file.", usage="ragtag_create_links.py -c components.fasta -b <hic.bam> -r <RE_site>")
    parser.add_argument("-a", metavar="FILE", default="", type=str, help="assembly fasta file [null]")
    parser.add_argument("-b", metavar="FILE", default="", type=str, help="Hi-C alignments in BAM format, sorted by read name [null]")
    parser.add_argument("-r", metavar="STR", default="GATC", type=str, help="CSV list of restriction sites or 'DNase' [GATC]")
    parser.add_argument("-p", metavar="FLOAT", default=1.0, type=float, help="portion of the sequence termini to consider for links [1.0]")
    parser.add_argument("--list-enzymes", action='store_true', default=False, help="list all available restriction enzymes/sites")

    args = parser.parse_args()

    # Print a restriction enzyme help message if requested
    if args.list_enzymes:
        RestrictionEnzymes.get_info()
        sys.exit()

    # Continue with normal functionality if no restriction enzyme help message is requested
    if not args.a or not args.b or not args.r:
        parser.print_help()
        sys.exit()

    # Set the terminus portion
    portion = args.p
    if not 1 >= portion > 0:
        raise ValueError("portion must be between 0 (exclusive) and 1 (inclusive)")

    asm_file = os.path.abspath(args.a)
    bam_file = os.path.abspath(args.b)

    dnase_mode = False
    re_string = args.r.upper()
    if "DNASE" in re_string:
        dnase_mode = True
        log("Running in DNase mode.")

    re_set = set()
    if not dnase_mode:
        re_set = set(filter(None, args.r.split(",")))
        if not re_set:
            raise ValueError("At least one restriction enzyme/site is needed (-r) if not using 'DNase'.")

    # Store the sequence lengths
    asm_lens = dict()
    fai = pysam.FastaFile(asm_file)
    for ref in fai.references:
        asm_lens[ref] = fai.get_reference_length(ref)
    fai.close()

    # Get the left and right cutoff positions for each sequence
    l_cutoffs = dict()
    r_cutoffs = dict()
    for ref in asm_lens:
        l = asm_lens[ref] // 2
        r = asm_lens[ref] - l

        l_cutoffs[ref] = round(l * portion)
        r_cutoffs[ref] = asm_lens[ref] - round(r * portion)

    # Get the raw Hi-C links
    log("Computing raw Hi-C links from: {}".format(bam_file))
    raw_links = count_links(bam_file, l_cutoffs, r_cutoffs)

    # Normalize the Hi-C links
    l_norm_factors = l_cutoffs
    r_norm_factors = r_cutoffs

    # Normalize by the number of restriction sites if not using DNase
    if not dnase_mode:
        l_norm_factors = dict()
        r_norm_factors = dict()

        # Set the restriction enzymes
        RE = RestrictionEnzymes(re_set)
        log("Using the following restriction sites:\n{}".format(str(RE)))

        log("Counting restriction sites")
        rfm = RestrictionFragmentMap(asm_file, RE)

        # Get the number of sites for each contig terminus (l/b and r/e)
        for ref in l_cutoffs:
            l_norm_factors[ref] = rfm.count_sites_lte(ref, l_cutoffs[ref])
            r_norm_factors[ref] = rfm.count_sites_gt(ref, r_cutoffs[ref])

    log("Normalizing raw Hi-C links")
    norm_links = normalize_links(raw_links, l_norm_factors, r_norm_factors)
    write_links(raw_links, norm_links)


if __name__ == "__main__":
    main()
