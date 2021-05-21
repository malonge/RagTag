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

import os
import sys
import argparse

import pysam

from ragtag_utilities.utilities import log, run_oae, get_ragtag_version
from ragtag_utilities.AlignmentReader import PAFReader
from ragtag_utilities.ContigAlignment import ContigAlignment
from ragtag_utilities.Aligner import Minimap2Aligner
from ragtag_utilities.Aligner import UnimapAligner
from ragtag_utilities.Aligner import NucmerAligner
from ragtag_utilities.ScaffoldGraph import AGPMultiScaffoldGraph
from ragtag_utilities.ScaffoldGraph import PatchScaffoldGraph
from ragtag_utilities.ScaffoldGraph import Alignment


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


def build_aln_scaffold_graph(ctg_alns, components_fn, max_term_dist):
    """
    Build a directed scaffold graph from filtered alignments
    :param ctg_alns: query sequence -> ContigAlignment object
    :param components_fn: name of FASTA file with all relevant sequences
    :param max_term_dist: maximum alignment distance from a sequence terminus
    :return: PatchScaffoldGraph
    """
    sg = PatchScaffoldGraph(components_fn)

    for query_seq in ctg_alns:
        als = ctg_alns[query_seq]
        als.sort_by_query()
        last = None
        last_reversed = False

        # Iterate over each alignment for this query sequence
        for i in range(als.num_alns):
            cur_reversed = False

            # For each reference/query alignment terminus, determine if it is close to the sequence terminus
            ref_left_end, ref_right_end = als.ref_start_end(i, max_term_dist)
            query_left_end, query_right_end = als.query_start_end(i, max_term_dist)

            # Determine if we are reversing the reference sequence of this alignment
            if ref_left_end and ref_right_end:
                # Entire reference aligns - have to look at strand in alignment
                if als.strands[i] == '-':
                    cur_reversed = True
            elif ref_left_end:
                # Beginning of reference - we would expect suffix of read if same strand
                if query_left_end and query_right_end:
                    if als.strands[i] == '-':
                        cur_reversed = True
                elif query_left_end:
                    cur_reversed = True
            else:
                # End of contig - we would expect prefix of read if same strand
                if query_left_end and query_right_end:
                    if als.strands[i] == '-':
                        cur_reversed = True
                elif query_right_end:
                    cur_reversed = True

            if last is not None:
                if als.ref_headers[last] != als.ref_headers[i]:
                    my_query_end_offset = als.ref_lens[last] - als.ref_ends[last]
                    if last_reversed:
                        my_query_end_offset = als.ref_starts[last]

                    their_query_start_offset = als.ref_starts[i]
                    if cur_reversed:
                        their_query_start_offset = als.ref_lens[i] - als.ref_ends[i]

                    my_query_end = als.query_ends[last] + my_query_end_offset
                    their_query_start = als.query_starts[i] - their_query_start_offset
                    overlap = my_query_end - their_query_start

                    if overlap <= als.ref_lens[last] and overlap <= als.ref_lens[i]:

                        # Determine the scaffold graph nodes
                        u = als.ref_headers[last] + "_e"
                        v = als.ref_headers[i] + "_b"
                        if last_reversed:
                            u = als.ref_headers[last] + "_b"
                        if cur_reversed:
                            v = als.ref_headers[i] + "_e"

                        alignment = Alignment(
                            u,
                            v,
                            query_seq,
                            als.query_len,
                            my_query_end,
                            their_query_start,
                            0,  # Always on the query's forward strand
                            is_gap=False
                        )
                        sg.add_edge(u, v, alignment)

            last = i
            last_reversed = cur_reversed

    sg.remove_heavier_than(1)
    return sg


def main():
    description = "Homology-based assembly patching: Make continuous joins and fill gaps " \
                  "in 'target.fa' using sequences from 'query.fa'"

    parser = argparse.ArgumentParser(description=description, usage="ragtag.py patch <target.fa> <query.fa>")

    parser.add_argument("reference", metavar="<target.fa>", nargs='?', default="", type=str, help="target fasta file (uncompressed or bgzipped)")
    parser.add_argument("query", metavar="<query.fa>", nargs='?', default="", type=str, help="query fasta file (uncompressed or bgzipped)")

    patch_options = parser.add_argument_group("patching")
    patch_options.add_argument("-e", metavar="<exclude.txt>", type=str, default="", help="list of target sequences to ignore [null]")
    patch_options.add_argument("-j", metavar="<skip.txt>", type=str, default="", help="list of query sequences to ignore [null]")
    patch_options.add_argument("-f", metavar="INT", type=int, default=1000, help="minimum unique alignment length [1000]")
    patch_options.add_argument("--remove-small", action="store_true", default=False, help="remove unique alignments shorter than '-f'")
    patch_options.add_argument("-q", metavar="INT", type=int, default=10, help="minimum mapq (NA for Nucmer alignments) [10]")
    patch_options.add_argument("-d", metavar="INT", type=int, default=100000, help="maximum alignment merge distance [100000]")
    patch_options.add_argument("-s", metavar="INT", type=int, default=50000, help="minimum merged alignment length [50000]")
    patch_options.add_argument("-i", metavar="FLOAT", type=float, default=0.05, help="maximum merged alignment distance from sequence terminus. fraction of the sequence length if < 1 [0.05]")
    patch_options.add_argument("--fill-only", action="store_true", default=False, help="only fill existing target gaps. do not join target sequences")
    patch_options.add_argument("--join-only", action="store_true", default=False, help="only join and patch target sequences. do not fill existing gaps")

    io_options = parser.add_argument_group("input/output options")
    io_options.add_argument("-o", metavar="PATH", type=str, default="ragtag_output", help="output directory [./ragtag_output]")
    io_options.add_argument("-w", action='store_true', default=False, help="overwrite intermediate files")
    io_options.add_argument("-u", action='store_true', default=False, help="add suffix to unplaced sequence headers")
    io_options.add_argument("--debug", action='store_true', default=False, help=argparse.SUPPRESS)

    aln_options = parser.add_argument_group("mapping options")
    aln_options.add_argument("-t", metavar="INT", type=int, default=1, help="number of minimap2/unimap threads [1]")
    aln_options.add_argument("--aligner", metavar="PATH", type=str, default="nucmer", help="aligner executable ('nucmer' (recommended), 'unimap' or 'minimap2') [nucmer]")
    mm2_default = "-x asm5"
    aln_options.add_argument("--mm2-params", metavar="STR", type=str, default=mm2_default, help="space delimited minimap2 parameters ['%s']" % mm2_default)
    aln_options.add_argument("--unimap-params", metavar="STR", type=str, default=mm2_default, help="space delimited unimap parameters ['%s']" % mm2_default)
    aln_options.add_argument("--nucmer-params", metavar="STR", type=str, default="--maxmatch -l 100 -c 500", help="space delimted nucmer parameters ['--maxmatch -l 100 -c 500']")

    args = parser.parse_args()
    if not args.reference or not args.query:
        parser.print_help()
        sys.exit("\n** The target and query FASTA files are required **")

    log("VERSION", "RagTag " + get_ragtag_version())
    log("WARNING", "This is a beta version of `ragtag patch`")
    log("CMD", "ragtag.py patch " + " ".join(sys.argv[1:]))

    reference_fn = os.path.abspath(args.reference)
    query_fn = os.path.abspath(args.query)

    # Check that the reference/query file exists
    if not os.path.isfile(reference_fn):
        raise FileNotFoundError("Could not find file: %s" % reference_fn)

    if not os.path.isfile(query_fn):
        raise FileNotFoundError("Could not find file: %s" % query_fn)

    # Alignment processing parameters
    min_ulen = args.f
    keep_small_uniques = not args.remove_small
    merge_dist = args.d
    num_threads = args.t

    aligner_path = args.aligner
    aligner = aligner_path.split("/")[-1]
    if aligner.split("/")[-1] not in {'minimap2', 'unimap', 'nucmer'}:
        raise ValueError("Must specify either 'minimap2', 'unimap', or 'nucmer' (PATHs allowed) with '--aligner'.")

    mm2_params = args.mm2_params
    unimap_params = args.unimap_params
    nucmer_params = args.nucmer_params

    # Mapq filtering parameters
    min_mapq = args.q
    if aligner == "nucmer":
        min_mapq = 0

    # Add the number of mm2/unimap threads if the mm2 params haven't been overridden.
    if mm2_params == mm2_default:
        mm2_params += " -t " + str(num_threads)
    if unimap_params == mm2_default:
        unimap_params += " -t " + str(num_threads)

    # Set reference/query sequences to ignore
    ref_blacklist = set()
    exclude_file = args.e
    if exclude_file:
        exclude_file = os.path.abspath(args.e)
        with open(exclude_file, "r") as f:
            for line in f:
                ref_blacklist.add(line.rstrip())

    query_blacklist = set()
    skip_file = args.j
    if skip_file:
        skip_file = os.path.abspath(skip_file)
        with open(skip_file, "r") as f:
            for line in f:
                query_blacklist.add(line.rstrip())

    # Supporting alignment parameters
    min_sup_aln_len = args.s
    max_term_dist = args.i
    if max_term_dist <= 0:
        raise ValueError("-i must be a positive nonzero number.")

    # Task options
    fill_only = args.fill_only
    join_only = args.join_only
    if fill_only and join_only:
        raise ValueError("'--fill-only' and '--join-only' cannot be used together")

    # I/O parameters
    add_suffix = args.u
    if not add_suffix:
        log("WARNING", "Without '-u' invoked, some component/object AGP pairs might share the same ID. Some external programs/databases don't like this. To ensure valid AGP format, use '-u'.")

    overwrite_files = args.w
    output_path = args.o
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = os.path.abspath(output_path) + "/"
    file_prefix = "ragtag.patch"

    # Setup a log file for external RagTag scripts
    ragtag_log = output_path + file_prefix + ".err"
    open(ragtag_log, "w").close()  # Wipe the log file

    # Debugging options
    debug_mode = args.debug

    # Break the reference assembly at gaps
    cmd = [
        "ragtag_splitasm.py",
        "-o",
        output_path + file_prefix + ".ctg.agp",
        reference_fn
    ]
    reference_ctg_fn = output_path + file_prefix + ".ctg.fasta"
    if os.path.isfile(reference_ctg_fn):
        if overwrite_files:
            log("INFO", "Overwriting pre-existing file: " + reference_ctg_fn)
            run_oae(cmd, reference_ctg_fn, ragtag_log)
        else:
            log("INFO", "Retaining pre-existing file: " + reference_ctg_fn)
    else:
        run_oae(cmd, reference_ctg_fn, ragtag_log)

    # Rename the query sequences
    cmd = [
        "ragtag_rename.py",
        query_fn,
        "-p",
        "qseq",
        "-o",
        output_path + file_prefix + ".rename.agp",
    ]
    query_rename_fn = output_path + file_prefix + ".rename.fasta"
    if os.path.isfile(query_rename_fn):
        if overwrite_files:
            log("INFO", "Overwriting pre-existing file: " + query_rename_fn)
            run_oae(cmd, query_rename_fn, ragtag_log)
        else:
            log("INFO", "Retaining pre-existing file: " + query_rename_fn)
    else:
        run_oae(cmd, query_rename_fn, ragtag_log)

    # Combine the reference contigs and query sequences to make a components fasta file
    components_fn = output_path + file_prefix + ".comps.fasta"
    if os.path.isfile(components_fn):
        if overwrite_files:
            log("INFO", "Overwriting pre-existing file: " + components_fn)
            write_comps = True
        else:
            log("INFO", "Retaining pre-existing file: " + components_fn)
            write_comps = False
    else:
        write_comps = True

    if write_comps:
        log("INFO", "Writing: " + components_fn)
        ref_fai = pysam.FastaFile(reference_ctg_fn)
        query_fai = pysam.FastaFile(query_rename_fn)
        with open(components_fn, "w") as f:
            for ref in ref_fai.references:
                f.write(">" + ref + "\n")
                f.write(ref_fai.fetch(ref) + "\n")

            for query in query_fai.references:
                f.write(">" + query + "\n")
                f.write(query_fai.fetch(query) + "\n")

    # Map the query assembly to the reference contigs
    log("INFO", "Mapping the query genome to the target genome")
    if aligner == "minimap2":
        al = Minimap2Aligner(reference_ctg_fn, [query_rename_fn], aligner_path, mm2_params, output_path + file_prefix + ".asm", in_overwrite=overwrite_files)
    elif aligner == "unimap":
        al = UnimapAligner(reference_ctg_fn, [query_rename_fn], aligner_path, unimap_params, output_path + file_prefix + ".asm", in_overwrite=overwrite_files)
    else:
        al = NucmerAligner(reference_ctg_fn, [query_rename_fn], aligner_path, nucmer_params, output_path + file_prefix + ".asm", in_overwrite=overwrite_files)
    al.run_aligner()

    # If alignments are from Nucmer, need to convert from delta to paf
    if aligner == "nucmer":
        cmd = ["ragtag_delta2paf.py", output_path + file_prefix + ".asm.delta"]
        run_oae(cmd, output_path + file_prefix + ".asm.paf", ragtag_log)

    # Read and organize the alignments
    log("INFO", "Reading whole genome alignments")
    # ctg_alns: query header -> ContigAlignment object
    ctg_alns = read_genome_alignments(output_path + file_prefix + ".asm.paf", query_blacklist, ref_blacklist)

    # Check if any alignments are left
    if not ctg_alns:
        raise RuntimeError("There are no alignments. Check '{}'.".format(output_path + file_prefix + ".asm.paf"))

    # Filter the alignments
    unfiltered_strings, filtered_strings, merged_strings, useful_strings = [], [], [], []
    log("INFO", "Filtering and merging alignments")
    fltrd_ctg_alns = dict()
    for i in ctg_alns:
        # Unique anchor filtering
        unfiltered_strings.append(str(ctg_alns[i]))
        ctg_alns[i] = ctg_alns[i].unique_anchor_filter(min_ulen, keep_small=keep_small_uniques)

        # mapq filtering
        if ctg_alns[i] is not None:
            ctg_alns[i] = ctg_alns[i].filter_mapq(min_mapq)
            if ctg_alns[i] is not None:
                filtered_strings.append(str(ctg_alns[i]))

                # alignment merging
                ctg_alns[i] = ctg_alns[i].merge_alns(merge_dist=merge_dist, careful_merge=True)
                if ctg_alns[i] is not None:
                    merged_strings.append(str(ctg_alns[i]))

                    # Length filtering
                    ctg_alns[i] = ctg_alns[i].filter_lengths(min_sup_aln_len)
                    if ctg_alns[i] is not None:
                        # terminal filtering
                        ctg_alns[i] = ctg_alns[i].keep_terminals(max_term_dist)

                        # Save the remaining useful alignments
                        if ctg_alns[i] is not None and ctg_alns[i].num_refs > 1 and not ctg_alns[i].has_internal_ref_cuttings(max_term_dist):
                            useful_strings.append(str(ctg_alns[i]))
                            fltrd_ctg_alns[i] = ctg_alns[i]

    # Write debugging files
    debug_non_fltrd_file = output_path + file_prefix + ".debug.unfiltered.paf"
    debug_fltrd_file = output_path + file_prefix + ".debug.filtered.paf"
    debug_merged_file = output_path + file_prefix + ".debug.merged.paf"
    debug_useful_file = output_path + file_prefix + ".debug.useful.paf"
    if debug_mode:
        with open(debug_non_fltrd_file, "w") as f:
            f.write("".join(unfiltered_strings))

        with open(debug_fltrd_file, "w") as f:
            f.write("".join(filtered_strings))

        with open(debug_merged_file, "w") as f:
            f.write("".join(merged_strings))

        with open(debug_useful_file, "w") as f:
            f.write("".join(useful_strings))

    # Make a Scaffold Graph encoding known reference contigs adjacencies
    log("INFO", "Building a scaffold graph from the contig AGP file")
    agp_multi_sg = AGPMultiScaffoldGraph(reference_ctg_fn)
    agp_multi_sg.add_agps([output_path + file_prefix + ".ctg.agp"])
    agp_sg = agp_multi_sg.merge()

    # As a hack, go through the AGP sg and make the required directed scaffold graph
    agp_psg = PatchScaffoldGraph(components_fn)
    for u, v in agp_sg.edges:
        aln = Alignment(
            u,
            v,
            "",
            agp_sg[u][v]["gap_size"][0],
            0,
            agp_sg[u][v]["gap_size"][0],
            0,
            is_gap=True
        )
        agp_psg.add_edge(u, v, aln)

    # Make a second directed scaffold graph from the alignments
    log("INFO", "Building a scaffold graph from the target/query mappings")
    aln_psg = build_aln_scaffold_graph(fltrd_ctg_alns, components_fn, max_term_dist)

    # Add edges for unfilled gaps
    for u, v in agp_psg.edges:
        if not aln_psg.has_edge(u, v):
            aln_psg.add_edge(u, v, agp_psg[u][v]["alignment"])

    # Remove known false edges
    for u, v in agp_psg.edges:
        for neighbor in aln_psg.neighbors(u):
            if neighbor != v:
                aln_psg.remove_edge(u, neighbor)
                aln_psg.remove_edge(neighbor, u)

        for neighbor in aln_psg.neighbors(v):
            if neighbor != u:
                aln_psg.remove_edge(neighbor, v)
                aln_psg.remove_edge(v, neighbor)

    # Adjust the graph depending on if only fills or joins are requested
    if fill_only:
        psg = PatchScaffoldGraph(components_fn)
        for u, v in agp_psg.edges:
            psg.add_edge(u, v, aln_psg[u][v]["alignment"])
            psg.add_edge(v, u, aln_psg[v][u]["alignment"])
        aln_psg = psg

    if join_only:
        for u, v in agp_psg.edges:
            aln_psg[u][v]["alignment"] = agp_psg[u][v]["alignment"]
            aln_psg[v][u]["alignment"] = agp_psg[v][u]["alignment"]

    if debug_mode:
        aln_psg.write_gml(output_path + file_prefix + ".debug.sg.gml")

    # Compute a matching solution for the graph
    log("INFO", "Computing a matching solution to the scaffold graph")
    match_psg = aln_psg.max_weight_matching()

    if debug_mode:
        match_psg.write_gml(output_path + file_prefix + ".debug.matching.gml")

    # Write the output in AGP format
    log("INFO", "Writing output files")
    match_psg.write_agp(output_path + file_prefix + ".agp", output_path + file_prefix + ".ctg.fasta", add_suffix_to_unplaced=add_suffix)

    # Write the output in fasta format
    cmd = [
        "ragtag_agp2fa.py",
        output_path + file_prefix + ".agp",
        components_fn
    ]
    run_oae(cmd, output_path + file_prefix + ".fasta", ragtag_log)

    log("INFO", "Goodbye")


if __name__ == "__main__":
    main()
