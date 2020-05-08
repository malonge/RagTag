#!/usr/bin/env python

import os
import sys
import argparse
from collections import defaultdict

import pysam

from ragoo2_utilities.utilities import log, run, run_o
from ragoo2_utilities.Aligner import Minimap2Aligner
from ragoo2_utilities.Aligner import NucmerAligner
from ragoo2_utilities.AlignmentReader import AlignmentReader
from ragoo2_utilities.ContigAlignment import ContigAlignment


def remove_contained(a):
    """
    remove contained intervals
    :param a: list of tuples (start, end, header)
    :return: intervals with contained intervals removed
    """
    o = []
    a = sorted(a, key=lambda x: (x[0], -x[1]))

    max_end = -1
    for i in a:
        if i[1] > max_end:
            max_end = i[1]
            o.append(i)
    return o


def write_orderings(out_file, query_file, ordering_dict, ctg_dict, gap_dict, gap_type_dict, make_chr0, overwrite):
    # Check if the output file already exists
    if os.path.isfile(out_file):
        if not overwrite:
            log("Retaining pre-existing file: " + out_file)
            return

    # Proceed with writing the intermediate output
    placed_seqs = set()
    gap_id = 0
    all_out_lines = []

    # Go through the reference sequences in sorted order
    sorted_ref_headers = sorted(list(ordering_dict.keys()))
    for ref_header in sorted_ref_headers:
        pos = 0
        new_ref_header = ref_header + "_RaGOO2"
        q_seqs = ordering_dict[ref_header]
        gap_seqs = gap_dict[ref_header]
        gap_types = gap_type_dict[ref_header]

        # Iterate through the query sequences for this reference header
        for i in range(len(q_seqs)):
            out_line = []
            q = q_seqs[i][2]
            placed_seqs.add(q)
            qlen = ctg_dict[q].query_len
            strand = ctg_dict[q].orientation
            gc, lc, oc = ctg_dict[q].grouping_confidence, ctg_dict[q].location_confidence, ctg_dict[q].orientation_confidence
            out_line.append(new_ref_header)
            out_line.append(str(pos))
            pos += qlen
            out_line.append(str(pos) + "\tS")
            out_line.append(q)
            out_line.append(strand)
            out_line.append(str(gc))
            out_line.append(str(lc))
            out_line.append(str(oc))
            all_out_lines.append("\t".join(out_line))

            if i < len(gap_seqs):
                # Print the gap line
                out_line = []
                out_line.append(new_ref_header)
                out_line.append(str(pos))
                pos += gap_seqs[i]
                gap_type = gap_types[i]
                out_line.append(str(pos) + "\t" + gap_type + "\t" + str(gap_id))
                gap_id += 1
                out_line.append("NA\tNA\tNA\tNA")
                all_out_lines.append("\t".join(out_line))

    # Write unplaced sequences
    fai = pysam.FastaFile(query_file)
    all_seqs = set(fai.references)
    unplaced_seqs = all_seqs - placed_seqs
    if unplaced_seqs:
        if make_chr0:
            pos = 0
            new_ref_header = "Chr0_RaGOO2"
            for q in unplaced_seqs:
                out_line = []
                qlen = fai.get_reference_length(q)
                out_line.append(new_ref_header)
                out_line.append(str(pos))
                pos += qlen
                out_line.append(str(pos) + "\tS")
                out_line.append(q)
                out_line.append("+")
                out_line.append("NA")
                out_line.append("NA")
                out_line.append("NA")
                all_out_lines.append("\t".join(out_line))

                # Now for the gap, since we are making a chr0
                out_line = []
                out_line.append(new_ref_header)
                out_line.append(str(pos))
                pos += 100
                out_line.append(str(pos) + "\tU\t" + str(gap_id))
                gap_id += 1
                out_line.append("NA\tNA\tNA\tNA")
                all_out_lines.append("\t".join(out_line))

            # Remove the final unecessary gap
            all_out_lines = all_out_lines[:-1]
        else:
            # List the unplaced contigs individually
            for q in unplaced_seqs:
                out_line = []
                qlen = fai.get_reference_length(q)
                gc, lc, oc = "NA", "NA", "NA"
                out_line.append(q + "_RaGOO2")
                out_line.append("0")
                out_line.append(str(qlen) + "\tS")
                out_line.append(q)
                out_line.append("+")
                out_line.append("NA")
                out_line.append("NA")
                out_line.append("NA")
                all_out_lines.append("\t".join(out_line))

    with open(out_file, "w") as f:
        f.write("\n".join(all_out_lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description='Reference-guided scaffolding', usage="ragoo2.py scaffold <reference.fa> <query.fa>")
    scaf_options = parser.add_argument_group("scaffolding options")
    scaf_options.add_argument("reference", metavar="<reference.fa>", nargs='?', default="", type=str, help="reference fasta file. must not be gzipped.")
    scaf_options.add_argument("query", metavar="<query.fa>", nargs='?', default="", type=str, help="query fasta file. must not be gzipped.")
    scaf_options.add_argument("-e", metavar="<exclude.txt>", type=str, default="", help="list of reference headers to ignore")
    scaf_options.add_argument("-j", metavar="<skip.txt>", type=str, default="", help="list of contigs to leave unplaced")
    scaf_options.add_argument("-f", metavar="INT", type=int, default=1000, help="minimum unique alignment length [1000]")
    scaf_options.add_argument("-d", metavar="INT", type=int, default=100000, help="alignment merge distance [100000]")
    scaf_options.add_argument("-i", metavar="FLOAT", type=float, default=0.2, help="minimum grouping confidence score [0.2]")
    scaf_options.add_argument("-a", metavar="FLOAT", type=float, default=0.0, help="minimum location confidence score [0.0]")
    scaf_options.add_argument("-s", metavar="FLOAT", type=float, default=0.0, help="minimum orientation confidence score [0.0]")
    scaf_options.add_argument("-C", action='store_true', default=False, help="concatenate unplaced contigs and make 'chr0'")
    scaf_options.add_argument("-r", action='store_true', default=False, help="infer gap sizes. if not, all gaps are 100 bp")
    scaf_options.add_argument("-g", metavar="INT", type=int, default=100, help="minimum inferred gap size [100]")
    scaf_options.add_argument("-m", metavar="INT", type=int, default=100000, help="maximum inferred gap size [100000]")

    io_options = parser.add_argument_group("input/output options")
    io_options.add_argument("-o", metavar="STR", type=str, default="ragoo2_output", help="output directory [ragoo2_output]")
    io_options.add_argument("-w", action='store_true', default=False, help="overwrite intermediate files")

    aln_options = parser.add_argument_group("mapping options")
    aln_options.add_argument("-t", metavar="INT", type=int, default=1, help="number of minimap2 threads [1]")
    aln_options.add_argument("--aligner", metavar="PATH", type=str, default="minimap2", help="aligner executable('nucmer' or 'minimap2') [minimap2]")
    aln_options.add_argument("--mm2-params", metavar="STR", type=str, default="-k19 -w19", help="space delimted minimap2 parameters ['-k19 -w19 -t1']")
    aln_options.add_argument("--nucmer-params", metavar="STR", type=str, default="-l 100 -c 500", help="space delimted nucmer parameters ['-l 100 -c 500']")

    # Get the command line arguments and ensure all paths are absolute.
    args = parser.parse_args()
    if not args.reference or not args.query:
        parser.print_help()
        sys.exit()

    log("CMD:" + " ".join(sys.argv))
    reference_file = os.path.abspath(args.reference)
    query_file = os.path.abspath(args.query)

    output_path = args.o.replace("/", "").replace(".", "")
    min_ulen = args.f
    merge_dist = args.d
    group_score_thresh = args.i
    loc_score_thresh = args.a
    orient_score_thresh = args.s
    make_chr0 = args.C
    infer_gaps = args.r
    overwrite_files = args.w
    num_threads = args.t

    min_gap_size = args.g
    max_gap_size = args.m
    if min_gap_size < 1:
        raise ValueError("the minimum gap size must be positive")

    if max_gap_size < 1:
        raise ValueError("the maximum gap size must be positive")

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

    # Add the number of mm2 threads if the mm2 params haven't been overridden.
    if mm2_params == "-k19 -w19":
        mm2_params += " -t" + str(num_threads)

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
    log("Mapping the query genome to the reference genome")
    if aligner == "minimap2":
        al = Minimap2Aligner(reference_file, query_file, aligner_path, mm2_params, output_path + "query_against_ref", in_overwrite=overwrite_files)
    else:
        al = NucmerAligner(reference_file, query_file, aligner_path, nucmer_params, output_path + "query_against_ref", in_overwrite=overwrite_files)
    al.run_aligner()

    # If alignments are from Nucmer, need to convert from delta to paf
    if aligner == "nucmer":
        cmd = ["ragoo2_delta2paf.py", output_path + "query_against_ref.delta", ">", output_path + "query_against_ref.paf"]
        run(cmd)

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
    log("Filtering and merging alignments")
    for i in ctg_alns:
        ctg_alns[i] = ctg_alns[i].unique_anchor_filter(min_ulen)
        if ctg_alns[i] is not None:
            ctg_alns[i] = ctg_alns[i].merge_alns(merge_dist=merge_dist)

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
    g_inferred = 0
    g_small = 0
    g_large = 0
    pad_sizes = dict()
    gap_types = dict()
    for i in mapped_ref_seqs:
        # Remove contained contigs and sort the rest
        non_contained = remove_contained(mapped_ref_seqs[i])
        mapped_ref_seqs[i] = sorted(non_contained)
        if infer_gaps:
            # Infer the gap sizes between adjacent query seqs
            # Use the primary alignments to infer gap sizes
            pad_sizes[i] = []
            gap_types[i] = []
            for j in range(1, len(mapped_ref_seqs[i])):
                # Get info for the upstream alignment
                left_ctg = mapped_ref_seqs[i][j - 1][2]
                left_ref_start, left_ref_end = fltrd_ctg_alns[left_ctg].get_best_ref_pos()
                left_qdist_start, left_qdist_end = fltrd_ctg_alns[left_ctg].get_best_q_dist()

                # Get info for the downstream alignment
                right_ctg = mapped_ref_seqs[i][j][2]
                right_ref_start, right_ref_end = fltrd_ctg_alns[right_ctg].get_best_ref_pos()
                right_qdist_start, right_qdist_end = fltrd_ctg_alns[right_ctg].get_best_q_dist()

                # Get the inferred gap size
                i_gap_size = (right_ref_start - right_qdist_start) - (left_ref_end + left_qdist_end)

                # Check if the gap size is too small or too large
                if i_gap_size <= min_gap_size:
                    pad_sizes[i].append(100)
                    gap_types[i].append("U")
                    g_small += 1
                elif i_gap_size > max_gap_size:
                    pad_sizes[i].append(100)
                    gap_types[i].append("U")
                    g_large += 1
                else:
                    pad_sizes[i].append(i_gap_size)
                    gap_types[i].append("N")
                    g_inferred += 1
        else:
            pad_sizes[i] = [100 for i in range(len(mapped_ref_seqs[i])-1)]
            gap_types[i] = ["U" for i in range(len(mapped_ref_seqs[i])-1)]

    if infer_gaps:
        log("%d inferred gap" % g_inferred)
        log("%d adjacent contig within min distance (%d) of each other" % (g_small, min_gap_size))
        log("%d inferred gaps exceed length threshold (%d)" % (g_large, max_gap_size))

    log("Writing: " + output_path + "scaffolding.placement.bed")
    # Write the intermediate output file
    write_orderings(output_path + "scaffolding.placement.bed", query_file, mapped_ref_seqs, fltrd_ctg_alns, pad_sizes, gap_types, make_chr0, overwrite_files)

    # Write the scaffolds
    log("Writing scaffolds")

    # Make the AGP file
    cmd = [
        "ragoo2_bed2agp.py",
        output_path + "scaffolding.placement.bed",
        output_path + "ragoo2.agp"
    ]
    run(cmd)

    # Build a FASTA from the AGP
    cmd = [
        "ragoo2_agp2fasta.py",
        output_path + "ragoo2.agp",
        query_file
    ]
    run_o(cmd, output_path + "ragoo2.fasta")

    # Calculate the stats
    cmd = [
        "ragoo2_stats.py",
        output_path + "scaffolding.placement.bed",
        output_path + "localization_stats.txt"
    ]
    run(cmd)


if __name__ == "__main__":
    main()
