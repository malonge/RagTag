#!/usr/bin/env python

import os
import argparse

from ragoo_utilities.utilities import log, run
from ragoo_utilities.Aligner import Minimap2Aligner
from ragoo_utilities.Aligner import NucmerAligner


def main():
    parser = argparse.ArgumentParser(description='Correct contigs according to alignments to a reference (v2.0.0)')
    parser.add_argument("reference", metavar="<reference.fasta>", type=str, help="reference fasta file. must not be gzipped.")
    parser.add_argument("query", metavar="<query.fasta>", type=str, help="query fasta file to be scaffolded. must not be gzipped.")
    parser.add_argument("-o", metavar="STR", type=str, default="ragoo_correct_output", help="output directory name [ragoo_correct_output]")
    parser.add_argument("--aligner", metavar="PATH", type=str, default="minimap2", help="Aligner ('nucmer' or 'minimap2') to use for correction. PATHs allowed [minimap2]")
    parser.add_argument("--mm2-params", metavar="STR", type=str, default="-k19 -w19 -t3", help="Space delimted parameters to pass directly to minimap2 ['-k19 -w19 -t3']")
    parser.add_argument("--nucmer-params", metavar="STR", type=str, default="-l 100 -c 500", help="Space delimted parameters to pass directly to nucmer ['-l 100 -c 500']")
    parser.add_argument("-f", metavar="INT", type=int, default=10000, help="minimum unique alignment length to use for scaffolding [10000]")
    parser.add_argument("-q", metavar="INT", type=int, default=-1, help="minimum mapping quality value for alignments. only pertains to minimap2 alignments [-1]")
    parser.add_argument("-d", metavar="INT", type=int, default=100000, help="merge contig alignments within this distance [100000]")
    parser.add_argument("-e", metavar="<exclude.txt>", type=str, default="", help="single column text file of reference headers to ignore")
    parser.add_argument("-j", metavar="<skip.txt>", type=str, default="", help="List of contigs to automatically leave uncorrected")
    parser.add_argument("--inter", action="store_true", default=False, help="Only break misassemblies between reference sequences")
    parser.add_argument("--intra", action="store_true", default=False, help="Only break misassemblies within reference sequences")
    parser.add_argument("--gff", metavar="<features.gff>", type=str, default="", help="Avoid breaking query sequences within any specified gff interval")
    parser.add_argument("-w", action='store_true', default=False,help="overwrite pre-existing intermediate files. ragoo.fasta will always be overwritten")

    # If one uses the same output directory name for both correct and scaffold, could that be ok without having to worry
    # about re-writing files? e.g.
    # ragoo_correct.py -o test_out
    # ragoo_scaffold.py -o test_out
    # Would that be a problem?

    args = parser.parse_args()
    reference_file = os.path.abspath(args.reference)
    query_file = os.path.abspath(args.query)
    output_path = args.o.replace("/", "").replace(".", "")
    min_ulen = args.f
    min_q = args.q
    overwrite_files = args.w

    skip_file = args.j
    if skip_file:
        skip_file = os.path.abspath(args.j)

    exclude_file = args.e
    if exclude_file:
        exclude_file = os.path.abspath(args.e)

    # Get aligner arguments
    aligner_path = args.aligner
    aligner = aligner_path.split("/")[-1]
    if aligner.split("/")[-1] not in {'minimap2', 'nucmer'}:
        raise ValueError("Must specify either 'minimap2' or 'nucmer' (PATHs allowed) with '--aligner'.")
    mm2_params = args.mm2_params
    nucmer_params = args.nucmer_params

    # Make sure no quality filtering for nucmer alignments
    if aligner == "nucmer":
        min_q = -1

    # Get the skip and exclude sets
    query_blacklist = set()
    if skip_file:
        with open(skip_file, "r") as f:
            for line in f:
                query_blacklist.add(line.rstrip())

    ref_blacklist = set()
    if exclude_file:
        with open(exclude_file, "r") as f:
            for line in f:
                ref_blacklist.add(line.rstrip())

    # Get the current working directory and output path
    cwd = os.getcwd()
    output_path = cwd + "/" + output_path + "/"
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # Align the query to the reference
    log("Aligning the query to the reference")
    if aligner == "minimap2":
        al = Minimap2Aligner(reference_file, query_file, aligner_path, mm2_params, output_path + "c_query_against_ref",
                             in_overwrite=overwrite_files)
    else:
        al = NucmerAligner(reference_file, query_file, aligner_path, nucmer_params, output_path + "c_query_against_ref",
                           in_overwrite=overwrite_files)
    al.run_aligner()


if __name__ == "__main__":
    main()