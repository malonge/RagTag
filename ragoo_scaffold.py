#!/usr/bin/env python

import os
import argparse
from collections import defaultdict

from ragoo2_utilities.utilities import log, run
from ragoo2_utilities.Aligner import Minimap2Aligner
from ragoo2_utilities.Aligner import NucmerAligner
from ragoo2_utilities.AlignmentReader import AlignmentReader
from ragoo2_utilities.ContigAlignment import ContigAlignment


def write_orderings(ordering_dict, ctg_dict, gap_dict, overwrite, out_path):
    out_file = out_path + "orderings.bed"

    # Check if the output file already exists
    if os.path.isfile(out_file):
        if not overwrite:
            log("retaining pre-existing file: " + out_file)
            return

    # Proceed with writing the intermediate output
    gap_id = 0
    all_out_lines = []

    # Go through the reference sequences in sorted order
    sorted_ref_headers = sorted(list(ordering_dict.keys()))
    for ref_header in sorted_ref_headers:
        pos = 0
        new_ref_header = ref_header + "_RaGOO"
        q_seqs = ordering_dict[ref_header]
        gap_seqs = gap_dict[ref_header]

        # Iterate through the query sequences for this reference header
        for i in range(len(q_seqs)):
            out_line = []
            q = q_seqs[i][2]
            qlen = ctg_dict[q].query_len
            strand = ctg_dict[q].orientation
            gc, lc, oc = ctg_dict[q].grouping_confidence, ctg_dict[q].location_confidence, ctg_dict[q].orientation_confidence
            out_line.append(new_ref_header)
            out_line.append(str(pos))
            pos += qlen
            out_line.append(str(pos) + "\ts")
            out_line.append(q)
            out_line.append(strand)
            out_line.append(str(gc))
            out_line.append(str(lc))
            out_line.append(str(oc))
            all_out_lines.append("\t".join(out_line))

            if i < len(gap_seqs):
                out_line = []
                # Print the gap line
                out_line.append(new_ref_header)
                out_line.append(str(pos))
                pos += gap_seqs[i]
                out_line.append(str(pos) + "\tg\t" + str(gap_id))
                gap_id += 1
                out_line.append("NA\tNA\tNA\tNA")
                all_out_lines.append("\t".join(out_line))

    log("Writing: " + out_file)
    with open(out_file, "w") as f:
        f.write("\n".join(all_out_lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description='Scaffold contigs according to alignments to a reference (v2.0.0)')
    parser.add_argument("reference", metavar="<reference.fasta>", type=str, help="reference fasta file. must not be gzipped.")
    parser.add_argument("query", metavar="<query.fasta>", type=str, help="query fasta file to be scaffolded. must not be gzipped.")
    parser.add_argument("-o", metavar="STR", type=str, default="ragoo_output", help="output directory name [ragoo_output]")
    parser.add_argument("--aligner", metavar="PATH", type=str, default="minimap2", help="Aligner ('nucmer' or 'minimap2') to use for scaffolding. PATHs allowed [minimap2]")
    parser.add_argument("--mm2-params", metavar="STR", type=str, default="-k19 -w19 -t3", help="Space delimted parameters to pass directly to minimap2 ['-k19 -w19 -t3']")
    parser.add_argument("--nucmer-params", metavar="STR", type=str, default="-l 100 -c 500", help="Space delimted parameters to pass directly to nucmer ['-l 100 -c 500']")
    parser.add_argument("-e", metavar="<exclude.txt>", type=str, default="", help="single column text file of reference headers to ignore")
    parser.add_argument("-j", metavar="<skip.txt>", type=str, default="", help="List of contigs to automatically leave unplaced")
    parser.add_argument("-g", metavar="INT", type=int, default=100, help="gap size for padding in pseudomolecules [100]")
    parser.add_argument("-l", metavar="INT", type=int, default=1000, help="minimum alignment length to use for scaffolding [1000]")
    parser.add_argument("-f", metavar="INT", type=int, default=0, help="minimum unique alignment length to use for scaffolding [0]")
    parser.add_argument("-q", metavar="INT", type=int, default=-1, help="minimum mapping quality value for alignments. only pertains to minimap2 alignments [-1]")
    parser.add_argument("-i", metavar="FLOAT", type=float, default=0.2, help="minimum grouping confidence score needed to be localized [0.2]")
    parser.add_argument("-a", metavar="FLOAT", type=float, default=0.0, help="minimum location confidence score needed to be localized [0.0]")
    parser.add_argument("-d", metavar="FLOAT", type=float, default=0.0, help="minimum orientation confidence score needed to be localized [0.0]")
    parser.add_argument("-C", action='store_true', default=False, help="write unplaced contigs individually instead of making a chr0")
    parser.add_argument("-r", action='store_true', default=False, help=argparse.SUPPRESS) # Infer gaps from reference - not ready
    parser.add_argument("-w", action='store_true', default=False, help="overwrite pre-existing intermediate files. ragoo.fasta will always be overwritten")

    # Get the command line arguments and ensure all paths are absolute.
    args = parser.parse_args()
    reference_file = os.path.abspath(args.reference)
    query_file = os.path.abspath(args.query)
    output_path = args.o.replace("/", "").replace(".", "")
    min_len = args.l
    min_ulen = args.f
    min_q = args.q
    gap_size = args.g
    group_score_thresh = args.i
    loc_score_thresh = args.a
    orient_score_thresh = args.d
    make_chr0 = not args.C
    infer_gaps = args.r
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
        al = Minimap2Aligner(reference_file, query_file, aligner_path, mm2_params, output_path + "query_against_ref", in_overwrite=overwrite_files)
    else:
        al = NucmerAligner(reference_file, query_file, aligner_path, nucmer_params, output_path + "query_against_ref", in_overwrite=overwrite_files)
    al.run_aligner()

    # If alignments are from Nucmer, need to convert from delta to paf
    if aligner == "nucmer":
        # TODO make direct call to executable
        cmd = ["python3", "delta2paf.py", output_path + "query_against_ref.delta", ">", output_path + "query_against_ref.paf"]
        run(" ".join(cmd))

    # Read and organize the alignments
    log('Reading alignments')
    # Read all of the alignments into a temporary structure to save time on ContigAlignment instantiation
    tmp_ctg_alns = dict()
    aln_reader = AlignmentReader(output_path + "query_against_ref.paf")
    for aln_line in aln_reader.parse_alignments():
        # Check that the contig and reference in this alignment are allowed.
        if aln_line.query_header not in query_blacklist and aln_line.ref_header not in ref_blacklist:
            if aln_line.query_header not in tmp_ctg_alns:
                tmp_ctg_alns[aln_line.query_header] = [aln_line.query_header, aln_line.query_len,
                                                                  [aln_line.ref_header], [aln_line.ref_len],
                                                                  [aln_line.ref_start], [aln_line.ref_end],
                                                                  [aln_line.query_start], [aln_line.query_end],
                                                                  [aln_line.strand], [aln_line.aln_len],
                                                                  [aln_line.mapq]]
            else:
                tmp_ctg_alns[aln_line.query_header][2].append(aln_line.ref_header)
                tmp_ctg_alns[aln_line.query_header][3].append(aln_line.ref_len)
                tmp_ctg_alns[aln_line.query_header][4].append(aln_line.ref_start)
                tmp_ctg_alns[aln_line.query_header][5].append(aln_line.ref_end)
                tmp_ctg_alns[aln_line.query_header][6].append(aln_line.query_start)
                tmp_ctg_alns[aln_line.query_header][7].append(aln_line.query_end)
                tmp_ctg_alns[aln_line.query_header][8].append(aln_line.strand)
                tmp_ctg_alns[aln_line.query_header][9].append(aln_line.aln_len)
                tmp_ctg_alns[aln_line.query_header][10].append(aln_line.mapq)

    ctg_alns = dict()
    for i in tmp_ctg_alns:
        ctg_alns[i] = ContigAlignment(
            tmp_ctg_alns[i][0],
            tmp_ctg_alns[i][1],
            tmp_ctg_alns[i][2],
            tmp_ctg_alns[i][3],
            tmp_ctg_alns[i][4],
            tmp_ctg_alns[i][5],
            tmp_ctg_alns[i][6],
            tmp_ctg_alns[i][7],
            tmp_ctg_alns[i][8],
            tmp_ctg_alns[i][9],
            tmp_ctg_alns[i][10]
        )

    # Filter the alignments
    log("Filtering alignments")
    for i in ctg_alns:
        ctg_alns[i] = ctg_alns[i].filter_mapq(min_q)
        if ctg_alns[i] is not None:
            ctg_alns[i] = ctg_alns[i].filter_lengths(min_len)
            if ctg_alns[i] is not None:
                ctg_alns[i] = ctg_alns[i].unique_anchor_filter(min_ulen)

    # Remove query sequences which have no more qualifying alignments
    fltrd_ctg_alns = dict()
    for i in ctg_alns:
        if ctg_alns[i] is not None:
            if all([
                ctg_alns[i].grouping_confidence > group_score_thresh,
                ctg_alns[i].location_confidence > loc_score_thresh,
                ctg_alns[i].orientation_confidence > orient_score_thresh
            ]):
                fltrd_ctg_alns[i] = ctg_alns[i]

    # For each reference sequence which has at least one assigned query sequence, get the list of
    # all query sequences assigned to that reference sequence.
    log("Ordering and orienting query sequences")
    mapped_ref_seqs = defaultdict(list)
    for i in fltrd_ctg_alns:
        best_ref = fltrd_ctg_alns[i].best_ref_header
        ref_start, ref_end = fltrd_ctg_alns[i].get_best_ref_pos()
        mapped_ref_seqs[best_ref].append((ref_start, ref_end, i))

    # Sort the query sequences for each reference sequence and define the padding sizes between adjacent query seqs
    pads_sizes = dict()
    for i in mapped_ref_seqs:
        mapped_ref_seqs[i] = sorted(mapped_ref_seqs[i])
        if infer_gaps:
            # Infer the gap sizes between adjacent query seqs
            pads_sizes[i] = []
            for j in range(1, len(mapped_ref_seqs[i])):
                left_ctg = mapped_ref_seqs[i][j-1][2]
                left_min, left_max = fltrd_ctg_alns[left_ctg].get_best_ref_flanks()
                right_ctg = mapped_ref_seqs[i][j][2]
                right_min, right_max = fltrd_ctg_alns[right_ctg].get_best_ref_flanks()

                # If the contigs overlap, revert to the fixed pre-defined gap size
                if right_min - left_max >= 0:
                    pads_sizes[i].append(right_min - left_max)
                else:
                    pads_sizes[i].append(gap_size)
        else:
            pads_sizes[i] = [gap_size for i in range(len(mapped_ref_seqs[i])-1)]

    # Write the intermediate output file
    write_orderings(mapped_ref_seqs, fltrd_ctg_alns, pads_sizes, overwrite_files, output_path)

    # Write the scaffolds
    log("Writing scaffolds")

    # TODO make direct call to the executable
    cmd = [
        "build_scaffolds.py",
        output_path + "orderings.bed",
        query_file,
        output_path + "ragoo.fasta",
        output_path + "unplaced.txt",
        str(gap_size)
    ]
    if not make_chr0:
        cmd.append("-C")
    run(" ".join(cmd))

    # Calculate the stats
    cmd = [
        "ragoo_stats.py",
        output_path + "orderings.bed",
        output_path + "unplaced.txt",
        output_path + "localization_stats.txt"
    ]
    run(" ".join(cmd))

    # Make the AGP file
    cmd = [
        "bed2agp.py",
        output_path + "orderings.bed",
        output_path + "unplaced.txt",
        output_path + "ragoo.agp",
        str(gap_size)
    ]
    if not make_chr0:
        cmd.append("-C")
    run(" ".join(cmd))


if __name__ == "__main__":
    main()

