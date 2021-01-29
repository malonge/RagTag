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

import pysam

from ragtag_utilities.utilities import log, run_oae, get_ragtag_version
from ragtag_utilities.AlignmentReader import PAFReader
from ragtag_utilities.ContigAlignment import ContigAlignment
from ragtag_utilities.AGPFile import AGPFile
from ragtag_utilities.Aligner import Minimap2Aligner
from ragtag_utilities.Aligner import NucmerAligner


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


def write_orderings(out_agp_file, out_confidence_file, query_file, ordering_dict, ctg_dict, gap_dict, gap_type_dict, make_chr0, overwrite, add_suffix):
    # Check if the output file already exists
    if os.path.isfile(out_agp_file):
        if not overwrite:
            log("Retaining pre-existing file: " + out_agp_file)
            return

        else:
            log("Overwriting pre-existing file: " + out_agp_file)

    # Proceed with writing the intermediate output
    placed_seqs = set()
    all_out_cs_lines = []  # For confidence scores
    agp = AGPFile(out_agp_file, mode="w")

    agp.add_pragma()
    agp.add_comment("# AGP created by RagTag {}".format(get_ragtag_version()))

    # Go through the reference sequences in sorted order
    sorted_ref_headers = sorted(list(ordering_dict.keys()))
    for ref_header in sorted_ref_headers:
        pid = 1
        pos = 0
        new_ref_header = ref_header + "_RagTag"
        q_seqs = ordering_dict[ref_header]
        gap_seqs = gap_dict[ref_header]
        gap_types = gap_type_dict[ref_header]

        # Iterate through the query sequences for this reference header
        for i in range(len(q_seqs)):
            out_agp_line = []
            out_cs_line = []
            q = q_seqs[i][2]
            placed_seqs.add(q)
            qlen = ctg_dict[q].query_len
            strand = ctg_dict[q].orientation
            gc, lc, oc = ctg_dict[q].grouping_confidence, ctg_dict[q].location_confidence, ctg_dict[q].orientation_confidence
            out_agp_line.append(new_ref_header)
            out_agp_line.append(str(pos+1))
            pos += qlen
            out_agp_line.append(str(pos))
            out_agp_line.append(str(pid))
            out_agp_line.append("W")
            out_agp_line.append(q)
            out_agp_line.append("1")
            out_agp_line.append(str(ctg_dict[q].query_len))
            out_agp_line.append(strand)

            # Save the confidence score info
            out_cs_line.append(q)
            out_cs_line.append(str(gc))
            out_cs_line.append(str(lc))
            out_cs_line.append(str(oc))

            agp.add_seq_line(*out_agp_line)
            all_out_cs_lines.append("\t".join(out_cs_line))
            pid += 1

            if i < len(gap_seqs):
                # Print the gap line
                out_agp_line = []
                out_agp_line.append(new_ref_header)
                out_agp_line.append(str(pos+1))
                pos += gap_seqs[i]
                out_agp_line.append(str(pos))
                out_agp_line.append(str(pid))
                gap_type = gap_types[i]
                out_agp_line.append(gap_type)
                out_agp_line.append(str(gap_seqs[i]))
                out_agp_line.append("scaffold")
                out_agp_line.append("yes")
                out_agp_line.append("align_genus")
                pid += 1
                agp.add_gap_line(*out_agp_line)

    # Write unplaced sequences
    fai = pysam.FastaFile(query_file)
    all_seqs = set(fai.references)
    unplaced_seqs = sorted(list(all_seqs - placed_seqs))
    if unplaced_seqs:
        if make_chr0:
            pos = 0
            pid = 1
            new_ref_header = "Chr0_RagTag"
            for q in unplaced_seqs:
                out_agp_line = []
                qlen = fai.get_reference_length(q)
                out_agp_line.append(new_ref_header)
                out_agp_line.append(str(pos+1))
                pos += qlen
                out_agp_line.append(str(pos))
                out_agp_line.append(str(pid))
                out_agp_line.append("W")
                out_agp_line.append(q)
                out_agp_line.append("1")
                out_agp_line.append(str(qlen))
                out_agp_line.append("+")

                agp.add_seq_line(*out_agp_line)
                pid += 1

                # Now for the gap, since we are making a chr0
                out_agp_line = []
                out_agp_line.append(new_ref_header)
                out_agp_line.append(str(pos+1))
                pos += 100
                out_agp_line.append(str(pos))
                out_agp_line.append(str(pid))
                out_agp_line.append("U")
                out_agp_line.append("100")
                out_agp_line.append("contig")
                out_agp_line.append("no")
                out_agp_line.append("na")

                agp.add_gap_line(*out_agp_line)
                pid += 1

            # Remove the final unecessary gap
            agp.pop_agp_line()
        else:
            # List the unplaced contigs individually
            for q in unplaced_seqs:
                out_agp_line = []
                qlen = fai.get_reference_length(q)
                if add_suffix:
                    out_agp_line.append(q + "_RagTag")
                else:
                    out_agp_line.append(q)
                out_agp_line.append("1")
                out_agp_line.append(str(qlen))
                out_agp_line.append("1")
                out_agp_line.append("W")
                out_agp_line.append(q)
                out_agp_line.append("1")
                out_agp_line.append(str(qlen))
                out_agp_line.append("+")
                agp.add_seq_line(*out_agp_line)

    agp.write()
    fai.close()

    # Write the confidence scores
    with open(out_confidence_file, "w") as f:
        f.write("query\tgrouping_confidence\tlocation_confidence\torientation_confidence\n")
        f.write("\n".join(all_out_cs_lines) + "\n")


def read_genome_alignments(aln_file, query_blacklist, ref_blacklist):
    tmp_ctg_alns = dict()
    aln_reader = PAFReader(aln_file)
    for aln_line in aln_reader.parse_alignments():
        # Check that the contig and reference in this alignment are allowed.
        if aln_line.query_header not in query_blacklist and aln_line.ref_header not in ref_blacklist:
            if aln_line.query_header not in tmp_ctg_alns:
                tmp_ctg_alns[aln_line.query_header] = [aln_line.query_header, aln_line.query_len,
                                                       [aln_line.query_start], [aln_line.query_end], [aln_line.strand],
                                                       [aln_line.ref_header], [aln_line.ref_len],
                                                       [aln_line.ref_start], [aln_line.ref_end],
                                                       [aln_line.num_match], [aln_line.aln_len],
                                                       [aln_line.mapq]]
            else:
                tmp_ctg_alns[aln_line.query_header][2].append(aln_line.query_start)
                tmp_ctg_alns[aln_line.query_header][3].append(aln_line.query_end)
                tmp_ctg_alns[aln_line.query_header][4].append(aln_line.strand)
                tmp_ctg_alns[aln_line.query_header][5].append(aln_line.ref_header)
                tmp_ctg_alns[aln_line.query_header][6].append(aln_line.ref_len)
                tmp_ctg_alns[aln_line.query_header][7].append(aln_line.ref_start)
                tmp_ctg_alns[aln_line.query_header][8].append(aln_line.ref_end)
                tmp_ctg_alns[aln_line.query_header][9].append(aln_line.num_match)
                tmp_ctg_alns[aln_line.query_header][10].append(aln_line.aln_len)
                tmp_ctg_alns[aln_line.query_header][11].append(aln_line.mapq)

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
            tmp_ctg_alns[i][10],
            tmp_ctg_alns[i][11]
        )
    return ctg_alns


def main():
    parser = argparse.ArgumentParser(description='Reference-guided scaffolding', usage="ragtag.py scaffold <reference.fa> <query.fa>")

    parser.add_argument("reference", metavar="<reference.fa>", nargs='?', default="", type=str, help="reference fasta file (uncompressed or bgzipped)")
    parser.add_argument("query", metavar="<query.fa>", nargs='?', default="", type=str, help="query fasta file (uncompressed or bgzipped)")

    scaf_options = parser.add_argument_group("scaffolding options")
    scaf_options.add_argument("-e", metavar="<exclude.txt>", type=str, default="", help="list of reference headers to ignore [null]")
    scaf_options.add_argument("-j", metavar="<skip.txt>", type=str, default="", help="list of query headers to leave unplaced [null]")
    scaf_options.add_argument("-f", metavar="INT", type=int, default=1000, help="minimum unique alignment length [1000]")
    scaf_options.add_argument("--remove-small", action="store_true", default=False, help="remove unique alignments shorter than -f")
    scaf_options.add_argument("-q", metavar="INT", type=int, default=10, help="minimum mapq (NA for Nucmer alignments) [10]")
    scaf_options.add_argument("-d", metavar="INT", type=int, default=100000, help="alignment merge distance [100000]")
    scaf_options.add_argument("-i", metavar="FLOAT", type=float, default=0.2, help="minimum grouping confidence score [0.2]")
    scaf_options.add_argument("-a", metavar="FLOAT", type=float, default=0.0, help="minimum location confidence score [0.0]")
    scaf_options.add_argument("-s", metavar="FLOAT", type=float, default=0.0, help="minimum orientation confidence score [0.0]")
    scaf_options.add_argument("-C", action='store_true', default=False, help="concatenate unplaced contigs and make 'chr0'")
    scaf_options.add_argument("-r", action='store_true', default=False, help="infer gap sizes. if not, all gaps are 100 bp")
    scaf_options.add_argument("-g", metavar="INT", type=int, default=100, help="minimum inferred gap size [100]")
    scaf_options.add_argument("-m", metavar="INT", type=int, default=100000, help="maximum inferred gap size [100000]")

    io_options = parser.add_argument_group("input/output options")
    io_options.add_argument("-o", metavar="PATH", type=str, default="ragtag_output", help="output directory [./ragtag_output]")
    io_options.add_argument("-w", action='store_true', default=False, help="overwrite intermediate files")
    io_options.add_argument("-u", action='store_true', default=False, help="add suffix to unplaced sequence headers")
    io_options.add_argument("--debug", action='store_true', default=False, help=argparse.SUPPRESS)

    aln_options = parser.add_argument_group("mapping options")
    aln_options.add_argument("-t", metavar="INT", type=int, default=1, help="number of minimap2 threads [1]")
    aln_options.add_argument("--aligner", metavar="PATH", type=str, default="minimap2", help="aligner executable ('nucmer' or 'minimap2') [minimap2]")
    mm2_default = "-x asm5"
    aln_options.add_argument("--mm2-params", metavar="STR", type=str, default=mm2_default, help="space delimited minimap2 parameters ['%s']" % mm2_default)
    aln_options.add_argument("--nucmer-params", metavar="STR", type=str, default="-l 100 -c 500", help="space delimted nucmer parameters ['-l 100 -c 500']")

    args = parser.parse_args()
    if not args.reference or not args.query:
        parser.print_help()
        print("\n** The reference and query FASTA files are required **")
        sys.exit()

    log("RagTag " + get_ragtag_version())
    log("CMD: ragtag.py scaffold " + " ".join(sys.argv[1:]))

    reference_file = os.path.abspath(args.reference)
    query_file = os.path.abspath(args.query)

    # Check that the reference/query file exists
    if not os.path.isfile(reference_file):
        raise ValueError("Could not find file: %s" % reference_file)

    if not os.path.isfile(query_file):
        raise ValueError("Could not find file: %s" % query_file)


    min_ulen = args.f
    keep_small_uniques = not args.remove_small
    merge_dist = args.d
    group_score_thresh = args.i
    loc_score_thresh = args.a
    orient_score_thresh = args.s
    make_chr0 = args.C
    infer_gaps = args.r
    num_threads = args.t

    # I/O options
    output_path = args.o
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = os.path.abspath(output_path) + "/"

    # Setup a log file for external RagTag scripts
    ragtag_log = output_path + "ragtag.scaffold.err"
    open(ragtag_log, "w").close()  # Wipe the log file

    overwrite_files = args.w
    remove_suffix = not args.u
    if remove_suffix:
        log("WARNING: Without '-u' invoked, some component/object AGP pairs might share the same ID. Some external programs/databases don't like this. To ensure valid AGP format, use '-u'.")

    # Gap options
    min_gap_size = args.g
    max_gap_size = args.m
    if min_gap_size < 1:
        raise ValueError("the minimum gap size must be positive")

    if max_gap_size < 1:
        raise ValueError("the maximum gap size must be positive")

    # Skip/exclude options
    query_blacklist = set()
    skip_file = args.j
    if skip_file:
        skip_file = os.path.abspath(args.j)
        with open(skip_file, "r") as f:
            for line in f:
                query_blacklist.add(line.rstrip())

    ref_blacklist = set()
    exclude_file = args.e
    if exclude_file:
        exclude_file = os.path.abspath(args.e)
        with open(exclude_file, "r") as f:
            for line in f:
                ref_blacklist.add(line.rstrip())

    # Get aligner arguments
    aligner_path = args.aligner
    aligner = aligner_path.split("/")[-1]
    if aligner.split("/")[-1] not in {'minimap2', 'nucmer'}:
        raise ValueError("Must specify either 'minimap2' or 'nucmer' (PATHs allowed) with '--aligner'.")

    mm2_params = args.mm2_params
    nucmer_params = args.nucmer_params

    # Mapq filtering params
    min_mapq = args.q
    if aligner == "nucmer":
        min_mapq = 0

    # Add the number of mm2 threads if the mm2 params haven't been overridden.
    if mm2_params == mm2_default:
        mm2_params += " -t " + str(num_threads)

    # Debugging options
    debug_mode = args.debug
    debug_non_fltrd_file = output_path + "ragtag.scaffolds.debug.unfiltered.paf"
    debug_fltrd_file = output_path + "ragtag.scaffolds.debug.filtered.paf"
    debug_merged_file = output_path + "ragtag.scaffolds.debug.merged.paf"
    debug_query_info_file = output_path + "ragtag.scaffolds.debug.query.info.txt"

    # Align the query to the reference
    log("Mapping the query genome to the reference genome")
    if aligner == "minimap2":
        al = Minimap2Aligner(reference_file, [query_file], aligner_path, mm2_params, output_path + "query_against_ref", in_overwrite=overwrite_files)
    else:
        al = NucmerAligner(reference_file, [query_file], aligner_path, nucmer_params, output_path + "query_against_ref", in_overwrite=overwrite_files)
    al.run_aligner()

    # If alignments are from Nucmer, need to convert from delta to paf
    if aligner == "nucmer":
        cmd = ["ragtag_delta2paf.py", output_path + "query_against_ref.delta"]
        run_oae(cmd, output_path + "query_against_ref.paf", ragtag_log)

    # Read and organize the alignments
    log('Reading whole genome alignments')
    # ctg_alns = dict :: key=query header, value=ContigAlignment object
    ctg_alns = read_genome_alignments(output_path + "query_against_ref.paf", query_blacklist, ref_blacklist)

    # Filter the alignments
    if debug_mode:
        # create new empty copies of debugging output files
        open(debug_non_fltrd_file, "w").close()
        open(debug_fltrd_file, "w").close()
        open(debug_merged_file, "w").close()
        open(debug_query_info_file, "w").close()

    log("Filtering and merging alignments")
    for i in ctg_alns:

        # Write unfiltered alignments
        if debug_mode:
            with open(debug_non_fltrd_file, "a") as f:
                f.write(str(ctg_alns[i]))

        ctg_alns[i] = ctg_alns[i].unique_anchor_filter(min_ulen, keep_small=keep_small_uniques)
        if ctg_alns[i] is not None:
            ctg_alns[i] = ctg_alns[i].filter_mapq(min_mapq)
            if ctg_alns[i] is not None:

                # Write filtered alignments
                if debug_mode:
                    with open(debug_fltrd_file, "a") as f:
                        f.write(str(ctg_alns[i]))

                ctg_alns[i] = ctg_alns[i].merge_alns(merge_dist=merge_dist)

    # Remove query sequences which have no more qualifying alignments
    fltrd_ctg_alns = dict()
    for i in ctg_alns:
        if ctg_alns[i] is not None:

            # Write merged alignments and confidence scores
            if debug_mode:
                with open(debug_merged_file, "a") as f:
                    f.write(str(ctg_alns[i]))

                with open(debug_query_info_file, "a") as f:
                    f.write("\t".join([
                        i,
                        ctg_alns[i].best_ref_header,
                        str(ctg_alns[i].grouping_confidence),
                        str(ctg_alns[i].location_confidence),
                        str(ctg_alns[i].orientation_confidence),
                    ]) + "\n")

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

    # Write the scaffolds
    log("Writing scaffolds")

    # Write the intermediate output file in AGP v2.1 format
    log("Writing: " + output_path + "ragtag.scaffolds.agp")
    write_orderings(output_path + "ragtag.scaffolds.agp", output_path + "ragtag.confidence.txt", query_file, mapped_ref_seqs, fltrd_ctg_alns, pad_sizes, gap_types, make_chr0, True, not remove_suffix)

    # Build a FASTA from the AGP
    cmd = [
        "ragtag_agp2fasta.py",
        output_path + "ragtag.scaffolds.agp",
        query_file
    ]
    run_oae(cmd, output_path + "ragtag.scaffolds.fasta", ragtag_log)

    # Calculate the stats
    cmd = [
        "ragtag_stats.py",
        output_path + "ragtag.scaffolds.agp",
        output_path + "ragtag.confidence.txt"
    ]
    run_oae(cmd, output_path + "ragtag.scaffolds.stats", ragtag_log)

    log("Goodbye")


if __name__ == "__main__":
    main()
