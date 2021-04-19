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
import numpy as np
import networkx as nx

from ragtag_utilities.utilities import log, run_oae, get_ragtag_version
from ragtag_utilities.AlignmentReader import PAFReader
from ragtag_utilities.ContigAlignment import ContigAlignment
from ragtag_utilities.AGPFile import AGPFile
from ragtag_utilities.Aligner import Minimap2Aligner
from ragtag_utilities.Aligner import UnimapAligner
from ragtag_utilities.Aligner import NucmerAligner
from ragtag_utilities.ScaffoldGraph import ScaffoldGraphBase
from ragtag_utilities.ScaffoldGraph import AGPMultiScaffoldGraph
from ragtag_utilities.ScaffoldGraph import MultiScaffoldGraph


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


def build_aln_scaffold_graph(components_fasta_fname, ctg_alns, max_term_dist):
    """
    Build a scaffold graph from alignments in PAF format,
    :param components_fasta_fname: reference FASTA file name
    :param ctg_alns: query header -> ContigAlignment object. Assumed to be filtered to only contain informative alignments
    :param max_term_dist: max distance of alignment to sequence terminus
    :return: A Scaffold Graph
    """
    G = MultiScaffoldGraph(components_fasta_fname)
    for query in ctg_alns:
        # Sort the alignments by query position
        ctg_alns[query].sort_by_query()

        # Iterate over each adjacent pair of alignments
        for i in range(1, ctg_alns[query].num_alns):
            # Only consider adjacent alignments to distinct reference sequences
            if ctg_alns[query].ref_headers[i-1] != ctg_alns[query].ref_headers[i]:
                left_strand, right_strand = ctg_alns[query].strands[i-1], ctg_alns[query].strands[i]

                # Get the potential nodes
                left_node = ctg_alns[query].ref_headers[i-1] + "_e"
                if left_strand == "-":
                    left_node = ctg_alns[query].ref_headers[i-1] + "_b"

                right_node = ctg_alns[query].ref_headers[i] + "_b"
                if right_strand == "-":
                    right_node = ctg_alns[query].ref_headers[i] + "_e"

                # Check that the alignments are close to the correct terminus
                left_pass = True
                if left_strand == "+":
                    if (ctg_alns[query].ref_lens[i-1] - ctg_alns[query].ref_ends[i-1]) > max_term_dist:
                        left_pass = False
                else:
                    if ctg_alns[query].ref_starts[i-1] > max_term_dist:
                        left_pass = False

                right_pass = True
                if right_strand == "-":
                    if (ctg_alns[query].ref_lens[i] - ctg_alns[query].ref_ends[i]) > max_term_dist:
                        right_pass = False
                else:
                    if ctg_alns[query].ref_starts[i] > max_term_dist:
                        right_pass = False

                # If everything is valid, collect the metadata and add nodes and an edge to the graph
                if left_pass and right_pass:
                    left_start_pos, right_end_pos = ctg_alns[query].ref_starts[i-1], ctg_alns[query].ref_ends[i]

                    left_seq, right_seq = ctg_alns[query].ref_headers[i-1], ctg_alns[query].ref_headers[i]
                    left_pos, right_pos = ctg_alns[query].ref_ends[i-1], ctg_alns[query].ref_starts[i]

                    query_header = ctg_alns[query].query_header
                    query_start = ctg_alns[query].query_ends[i-1]
                    query_end = ctg_alns[query].query_starts[i]

                    G.add_edge(
                        left_node,
                        right_node,
                        weight=1,
                        source_fname="alns.paf",
                        is_known_gap_size=True,
                        gap_size=query_end - query_start,
                        gap_type="scaffold",
                        linkage=True,
                        linkage_evidence="align_genus",
                        seqs=[left_seq, right_seq],
                        pos=[left_pos, right_pos],
                        term_pos=[left_start_pos, right_end_pos],
                        query=query_header,
                        query_start=query_start,
                        query_end=query_end
                    )

    SG = G.merge()

    # Remove edges with more than one supporting alignment
    SG.filter_non_one()
    return SG


def get_maximal_matching(sg):
    """ Find a solution to a Scaffold Graph. """
    matching = sg.get_max_weight_matching()

    cover_graph = nx.Graph()
    cover_graph.add_nodes_from(sg.nodes(data=True))

    # Add edges connecting contig ends
    node_base_set = set([i[:-2] for i in list(cover_graph.nodes)])
    for node in node_base_set:
        cover_graph.add_edge(node + "_b", node + "_e", weight=np.inf)

    # Add the Scaffold Graph edges that form the matching
    for u, v in matching:
        cover_graph.add_edge(u, v, **sg[u][v])

    # Remove any potential cycles
    edges_to_delete = []
    for cc in nx.connected_components(G=cover_graph):
        cc = cover_graph.subgraph(cc).copy()
        if cc.number_of_nodes() == cc.number_of_edges():
            assembly_edges = cc.edges(data=True)
            edges_to_delete.append(min(assembly_edges, key=lambda entry: entry[2]["weight"]))
    for u, v, data in edges_to_delete:
        cover_graph.remove_edge(u, v)

    # Remove intra-sequence edges
    for node in node_base_set:
        cover_graph.remove_edge(node + "_b", node + "_e")

    return cover_graph


def write_agp_solution(cover_graph, agp_scaffold_graph, scaffold_graph, agp_fname, add_suffix_to_unplaced=False, join_only=False):
    if not isinstance(agp_scaffold_graph, ScaffoldGraphBase):
        raise TypeError("agp_scaffold_graph must be an instance of ScaffoldGraph")

    if not isinstance(scaffold_graph, ScaffoldGraphBase):
        raise TypeError("scaffold_graph must be an instance of ScaffoldGraph")

    placed_components = set()

    # Initialize the AGP file
    agp = AGPFile(agp_fname, mode="w")
    agp.add_pragma()
    agp.add_comment("# AGP created by RagTag {}".format(get_ragtag_version()))

    # Enter the master loop to find joins until we run out
    obj_header_idx = 0
    while True:
        left_edge, start_comp_node = None, None
        obj_id, obj_pos = 1, 0
        obj_header = "scf" + "{0:08}".format(obj_header_idx) + "_RagTag"

        # Iterate over edges until we find a starting point that we haven't seen yet
        for u, v in cover_graph.edges:
            u_base, v_base = u[:-2], v[:-2]
            if u_base in placed_components:
                continue

            u_base_degree = cover_graph.degree[u_base + "_b"] + cover_graph.degree[u_base + "_e"]
            v_base_degree = cover_graph.degree[v_base + "_b"] + cover_graph.degree[v_base + "_e"]
            assert u_base_degree in {1, 2} and v_base_degree in {1, 2}

            # Check if one of these is a terminal component
            if u_base_degree == 1:
                start_comp_node = u
                left_edge = (u, v)
                break

            if v_base_degree == 1:
                start_comp_node = v
                left_edge = (u, v)
                break

        # If the above 'for' loop doesn't find any starting nodes, we are done.
        if left_edge is None:
            break

        start_comp = start_comp_node[:-2]
        start_comp_len = scaffold_graph.get_component_len(start_comp)

        comp_strand = "+"
        if start_comp_node.endswith("_b"):
            comp_strand = "-"

        # Get the position for the first component
        start_comp_pos_idx = cover_graph[left_edge[0]][left_edge[1]]["seqs"][0].index(start_comp)
        start_comp_pos = cover_graph[left_edge[0]][left_edge[1]]["pos"][0][start_comp_pos_idx]
        start_comp_term_pos = cover_graph[left_edge[0]][left_edge[1]]["term_pos"][0][start_comp_pos_idx]

        comp_start, comp_end = 0, start_comp_pos
        if comp_strand == "-":
            comp_start, comp_end = start_comp_term_pos, start_comp_len

        comp_len = comp_end - comp_start
        agp.add_seq_line(obj_header, obj_pos+1, obj_pos+comp_len, obj_id, "W", start_comp, comp_start+1, comp_end, comp_strand)
        obj_pos += comp_len
        obj_id += 1
        placed_components.add(start_comp)

        # iterate over each pair of connected adjacency edges
        prev_comp_node = start_comp_node
        while True:
            # Look for a right edge
            curr_comp_node = (set(left_edge) - {prev_comp_node}).pop()
            curr_comp = curr_comp_node[:-2]
            next_node_degree = cover_graph.degree[curr_comp + "_b"] + cover_graph.degree[curr_comp + "_e"]
            # Break if we don't find a non-terminal right edge
            if next_node_degree == 1:
                break

            # If we found a right edge, add the resulting patch and ref agp lines
            curr_node_oppo = curr_comp + "_b" if curr_comp_node.endswith("_e") else curr_comp + "_e"
            right_edge = (curr_node_oppo, list(cover_graph[curr_node_oppo].keys())[0])

            # Write the patch provided by the left edge
            fill = True
            if "is_gap" in cover_graph[left_edge[0]][left_edge[1]]:
                if cover_graph[left_edge[0]][left_edge[1]]["is_gap"][0]:
                    if cover_graph[left_edge[0]][left_edge[1]]["is_filled"][0]:
                        if join_only:
                            fill = False
                    else:
                        fill = False

            if fill:
                query = cover_graph[left_edge[0]][left_edge[1]]["query"][0]
                query_start, query_end = cover_graph[left_edge[0]][left_edge[1]]["query_start"][0], cover_graph[left_edge[0]][left_edge[1]]["query_end"][0]
                query_len = query_end - query_start

                # Check if the ref contigs overlap or bookend
                if query_len > 0:
                    agp.add_seq_line(obj_header, obj_pos+1, obj_pos+query_len, obj_id, "W", query, query_start+1, query_end, "+")
                    obj_pos += query_len
                    obj_id += 1
            else:
                # Add a gap line
                query_len = cover_graph[left_edge[0]][left_edge[1]]["gap_size"][0]
                agp.add_gap_line(obj_header,
                                 obj_pos+1,
                                 obj_pos+query_len,
                                 obj_id,
                                 "N" if cover_graph[left_edge[0]][left_edge[1]]["is_known_gap_size"][0] else "U",
                                 query_len,
                                 cover_graph[left_edge[0]][left_edge[1]]["gap_type"][0],
                                 "yes" if cover_graph[left_edge[0]][left_edge[1]]["linkage"][0] else "no",
                                 cover_graph[left_edge[0]][left_edge[1]]["linkage_evidence"][0]
                            )
                obj_pos += query_len
                obj_id += 1

            # Write the reference sequence between the edges
            start_pos_idx = cover_graph[left_edge[0]][left_edge[1]]["seqs"][0].index(curr_comp)
            start_pos = cover_graph[left_edge[0]][left_edge[1]]["pos"][0][start_pos_idx]

            if query_len < 0:
                start_pos = start_pos + abs(query_len)

            end_pos_idx = cover_graph[right_edge[0]][right_edge[1]]["seqs"][0].index(curr_comp)
            end_pos = cover_graph[right_edge[0]][right_edge[1]]["pos"][0][end_pos_idx]

            curr_comp_strand = "+"
            if curr_node_oppo.endswith("_b"):
                curr_comp_strand = "-"

            curr_comp_len = end_pos - start_pos
            agp.add_seq_line(obj_header, obj_pos+1, obj_pos+curr_comp_len, obj_id, "W", curr_comp, start_pos+1, end_pos, curr_comp_strand)
            obj_pos += curr_comp_len
            obj_id += 1
            placed_components.add(curr_comp)

            left_edge = right_edge
            prev_comp_node = curr_node_oppo

        # Finish the object by adding the other terminal component
        fill = True
        if "is_gap" in cover_graph[left_edge[0]][left_edge[1]]:
            if cover_graph[left_edge[0]][left_edge[1]]["is_gap"][0]:
                if cover_graph[left_edge[0]][left_edge[1]]["is_filled"][0]:
                    if join_only:
                        fill = False
                else:
                    fill = False

        if fill:
            query = cover_graph[left_edge[0]][left_edge[1]]["query"][0]
            query_start, query_end = cover_graph[left_edge[0]][left_edge[1]]["query_start"][0], cover_graph[left_edge[0]][left_edge[1]]["query_end"][0]
            query_len = query_end - query_start

            if query_len > 0:
                agp.add_seq_line(obj_header, obj_pos+1, obj_pos + query_len, obj_id, "W", query, query_start+1, query_end, "+")
                obj_pos += query_len
                obj_id += 1
        else:
            # Add a gap line
            query_len = cover_graph[left_edge[0]][left_edge[1]]["gap_size"][0]
            agp.add_gap_line(obj_header,
                             obj_pos + 1,
                             obj_pos + query_len,
                             obj_id,
                             "N" if cover_graph[left_edge[0]][left_edge[1]]["is_known_gap_size"][0] else "U",
                             query_len,
                             cover_graph[left_edge[0]][left_edge[1]]["gap_type"][0],
                             "yes" if cover_graph[left_edge[0]][left_edge[1]]["linkage"][0] else "no",
                             cover_graph[left_edge[0]][left_edge[1]]["linkage_evidence"][0]
                             )
            obj_pos += query_len
            obj_id += 1

        # Write the final terminal component
        comp_len = scaffold_graph.get_component_len(curr_comp)
        comp_pos_idx = cover_graph[left_edge[0]][left_edge[1]]["seqs"][0].index(curr_comp)
        comp_pos = cover_graph[left_edge[0]][left_edge[1]]["pos"][0][comp_pos_idx]
        comp_term_pos = cover_graph[left_edge[0]][left_edge[1]]["term_pos"][0][comp_pos_idx]

        curr_comp_strand = "+"
        start_pos = comp_pos
        end_pos = comp_len
        if curr_comp_node.endswith("_e"):
            curr_comp_strand = "-"
            start_pos = 0
            end_pos = comp_term_pos

        comp_len = end_pos - start_pos
        agp.add_seq_line(obj_header, obj_pos+1, obj_pos+comp_len, obj_id, "W", curr_comp, start_pos+1, end_pos, curr_comp_strand)
        placed_components.add(curr_comp)

        obj_header_idx += 1

    # Write all unplaced contigs
    remaining_components = scaffold_graph.components - placed_components
    for c in remaining_components:
        agp.add_seq_line(
            c + "_RagTag" * add_suffix_to_unplaced,
            "1",
            str(scaffold_graph.get_component_len(c)),
            "1",
            "W",
            c,
            "1",
            str(scaffold_graph.get_component_len(c)),
            "+"
        )

    agp.write()


def main():
    parser = argparse.ArgumentParser(description='Use a query assembly to patch a reference assembly', usage="ragtag.py patch <reference.fa> <query.fa>")

    parser.add_argument("reference", metavar="<reference.fa>", nargs='?', default="", type=str, help="reference fasta file (uncompressed or bgzipped)")
    parser.add_argument("query", metavar="<query.fa>", nargs='?', default="", type=str, help="query fasta file (uncompressed or bgzipped)")

    patch_options = parser.add_argument_group("patching")
    patch_options.add_argument("-e", metavar="<exclude.txt>", type=str, default="", help="list of reference sequences to ignore [null]")
    patch_options.add_argument("-j", metavar="<skip.txt>", type=str, default="", help="list of query sequences to ignore [null]")
    patch_options.add_argument("-f", metavar="INT", type=int, default=1000, help="minimum unique alignment length [1000]")
    patch_options.add_argument("--remove-small", action="store_true", default=False, help="remove unique alignments shorter than '-f'")
    patch_options.add_argument("-q", metavar="INT", type=int, default=10, help="minimum mapq (NA for Nucmer alignments) [10]")
    patch_options.add_argument("-d", metavar="INT", type=int, default=100000, help="maximum reference alignment merge distance [100000]")
    patch_options.add_argument("--careful-merge", action="store_true", default=False, help="apply '-d' to the query and reference coordinates")  # TODO  make this default
    patch_options.add_argument("-s", metavar="INT", type=int, default=50000, help="minimum merged alignment length [50000]")
    patch_options.add_argument("-i", metavar="INT", type=int, default=1000, help="maximum merged alignment distance from sequence terminus [1000]")
    patch_options.add_argument("--fill-only", action="store_true", default=False, help="only fill existing reference gaps. do not join reference sequnces")
    patch_options.add_argument("--join-only", action="store_true", default=False, help="only join and patch reference sequences. do not fill existing gaps")

    io_options = parser.add_argument_group("input/output options")
    io_options.add_argument("-o", metavar="PATH", type=str, default="ragtag_output", help="output directory [./ragtag_output]")
    io_options.add_argument("-w", action='store_true', default=False, help="overwrite intermediate files")
    io_options.add_argument("-u", action='store_true', default=False, help="add suffix to unplaced sequence headers")
    io_options.add_argument("--debug", action='store_true', default=False, help=argparse.SUPPRESS)

    aln_options = parser.add_argument_group("mapping options")
    aln_options.add_argument("-t", metavar="INT", type=int, default=1, help="number of minimap2 threads [1]")
    aln_options.add_argument("--aligner", metavar="PATH", type=str, default="minimap2", help="aligner executable ('nucmer', 'unimap' or 'minimap2') [minimap2]")
    mm2_default = "-x asm5"
    aln_options.add_argument("--mm2-params", metavar="STR", type=str, default=mm2_default, help="space delimited minimap2 parameters ['%s']" % mm2_default)
    aln_options.add_argument("--unimap-params", metavar="STR", type=str, default=mm2_default, help="space delimited unimap parameters ['%s']" % mm2_default)
    aln_options.add_argument("--nucmer-params", metavar="STR", type=str, default="-l 100 -c 500", help="space delimted nucmer parameters ['-l 100 -c 500']")

    args = parser.parse_args()
    if not args.reference or not args.query:
        parser.print_help()
        print("\n** The reference and query FASTA files are required **")
        sys.exit()

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
    careful_merge = args.careful_merge
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

    # Combine the reference contigs and query sequencs to make a components fasta file
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
    log("INFO", "Mapping the query genome to the reference genome")
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
                        if ctg_alns[i] is not None and ctg_alns[i].num_refs > 1:
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

    # Make a scaffold graph from the alignments to the query
    log("INFO", "Building a scaffold graph from the assembly alignments")
    sg = build_aln_scaffold_graph(reference_ctg_fn, fltrd_ctg_alns, max_term_dist)

    # Remove known false edges from the scaffold graph
    for u, v in agp_sg.edges:
        if u in sg.nodes:
            for neighbor in sg.neighbors(u):
                if neighbor != v:
                    sg.remove_edge(u, neighbor)

        if v in sg.nodes:
            for neighbor in sg.neighbors(v):
                if neighbor != u:
                    sg.remove_edge(v, neighbor)

    if debug_mode:
        agp_sg.connect_and_write_gml(output_path + file_prefix + ".debug.agp_sg.gml")
        sg.connect_and_write_gml(output_path + file_prefix + ".debug.sg.gml")

    log("INFO", "Computing solution to the scaffold graph")
    cover_graph = get_maximal_matching(sg)

    # Add all of the agp edges to the cover graph
    for u, v in agp_sg.edges:
        if cover_graph.has_edge(u, v):
            cover_graph[u][v]["is_gap"] = [True]
            cover_graph[u][v]["is_filled"] = [True]
            cover_graph[u][v]["agp_is_known_gap_size"] = agp_sg[u][v]["is_known_gap_size"]
            cover_graph[u][v]["agp_gap_size"] = agp_sg[u][v]["gap_size"]
            cover_graph[u][v]["agp_gap_type"] = agp_sg[u][v]["gap_type"]
            cover_graph[u][v]["agp_linkage"] = agp_sg[u][v]["linkage"]
            cover_graph[u][v]["agp_linkage_evidence"] = agp_sg[u][v]["linkage_evidence"]
        else:
            new_data = {
                "is_gap": [True],
                "is_filled": [False],
                "agp_is_known_gap_size": agp_sg[u][v]["is_known_gap_size"],
                "agp_gap_size": agp_sg[u][v]["gap_size"],
                "agp_gap_type": agp_sg[u][v]["gap_type"],
                "agp_linkage": agp_sg[u][v]["linkage"],
                "agp_linkage_evidence": agp_sg[u][v]["linkage_evidence"]
            }
            cover_graph.add_edge(u, v, **new_data)

    # If only filling gaps, remove all edges that make joins
    if fill_only:
        new_cover_graph = nx.Graph()
        new_cover_graph.add_nodes_from(cover_graph.nodes(data=True))
        for u, v in cover_graph.edges:
            if "is_gap" in cover_graph[u][v]:
                if cover_graph[u][v]["is_gap"][0]:
                    data = cover_graph[u][v]
                    new_cover_graph.add_edge(u, v, **cover_graph[u][v])
        cover_graph = new_cover_graph

    if debug_mode:
        gml_G = cover_graph.copy()
        # networkx doesn't like writing non-string attributes to GML
        for u, v in gml_G.edges:
            for key in gml_G[u][v]:
                gml_G[u][v][key] = str(gml_G[u][v][key])
        nx.readwrite.gml.write_gml(gml_G, output_path + file_prefix + ".debug.covergraph.gml")

    # Write the scaffolding output to an AGP file
    log("INFO", "Writing scaffolding results")
    write_agp_solution(cover_graph, agp_sg, sg, output_path + file_prefix + ".agp", add_suffix_to_unplaced=add_suffix, join_only=join_only)

    # Build a FASTA from the AGP
    cmd = [
        "ragtag_agp2fasta.py",
        output_path + file_prefix + ".agp",
        components_fn
    ]
    run_oae(cmd, output_path + file_prefix + ".fasta", ragtag_log)

    log("INFO", "Goodbye")


if __name__ == "__main__":
    main()
