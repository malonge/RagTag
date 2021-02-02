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

import networkx as nx
import numpy as np

from ragtag_utilities.AGPFile import AGPFile
from ragtag_utilities.ScaffoldGraph import ScaffoldGraphBase
from ragtag_utilities.ScaffoldGraph import MultiScaffoldGraph
from ragtag_utilities.ScaffoldGraph import AGPMultiScaffoldGraph
from ragtag_utilities.RestrictionFragmentMap import RestrictionEnzymes
from ragtag_utilities.utilities import log, get_ragtag_version, run_oae


def build_hic_graph(hic_file, components_fasta_fname):
    G = MultiScaffoldGraph(components_fasta_fname)

    # Iterate over the hic links and build the salsa graph
    with open(hic_file, "r") as f:
        for line in f:
            fields = line.rstrip().split("\t")
            n1, n2, w = fields[0], fields[1], float(fields[2])

            G.add_edge(
                n1,
                n2,
                weight=w,
                source_fname=hic_file,
                is_known_gap_size=False,
                gap_size=100,
                gap_type="scaffold",
                linkage=True,
                linkage_evidence="proximity_ligation"
            )

    SG = G.merge()
    return SG.best_buddy_scale()


def get_maximal_matching(sg):
    """ Find a scaffolding solution given a Scaffold Graph"""
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

    return cover_graph


def get_gap_size(size_list, func):
    if func not in {"MIN", "MAX", "MEAN"}:
        raise ValueError("{} not a valid function for determining gap size".format(func))

    if func == "MIN":
        return min(size_list)

    if func == "MAX":
        return max(size_list)

    return int(np.mean(size_list))


def write_agp_solution(cover_graph, scaffold_graph, agp_fname, gap_func="MIN", add_suffix_to_unplaced=False):
    """
    Here, we work with two graphs: A cover_graph and a scaffold_graph. A covergrpah defines a solution to the scaffold
    graph, and nodes from the same component are connected for convenience.

    We use the scaffold_graph for any original scaffold_graph info/functionality
    """
    if not isinstance(scaffold_graph, ScaffoldGraphBase):
        raise TypeError("scaffold_graph must be an instance of ScaffoldGraph")

    placed_components = set()

    # Iterate over each connected component
    agp = AGPFile(agp_fname, mode="w")
    agp.add_pragma()
    agp.add_comment("# AGP created by RagTag {}".format(get_ragtag_version()))

    # Iterate through the connected components
    for i, cc in enumerate(nx.connected_components(G=cover_graph)):
        # Sort the list of nodes for deterministic output
        cc = sorted(list(cc))
        obj_header = "scf" + "{0:05}".format(i+1) + "_RagTag"
        current_node = None

        # Iterate over each node in the connected component until we find a node with degree=1
        for node in cc:
            if cover_graph.degree[node] == 1:
                current_node = node
                break

        assert current_node is not None

        # Starting with the degree=1 node, build the AGP object from nodes in the path.
        visited_nodes = {current_node}
        degree = 0
        obj_id = 1
        obj_pos = 0

        # Traverse the component until we find the other end node
        while degree != 1:
            conn_nodes = set(cover_graph.neighbors(current_node))
            next_node = (conn_nodes - visited_nodes).pop()
            degree = cover_graph.degree[next_node]
            comp_len = scaffold_graph.get_component_len(next_node[:-2])

            # Check if this is an intra or inter sequence edge
            orientation = "+"
            if next_node[:-2] == current_node[:-2]:
                if next_node.endswith("_b"):
                    orientation = "-"
                    assert current_node.endswith("_e")

                agp.add_seq_line(
                    obj_header,
                    str(obj_pos + 1),
                    str(obj_pos + comp_len),
                    str(obj_id),
                    "W",
                    next_node[:-2],
                    "1",
                    str(comp_len),
                    orientation
                )
                obj_pos += comp_len
                placed_components.add(next_node[:-2])
            else:
                # Organize the gap info
                adjacency_data = scaffold_graph[current_node][next_node]

                # AGP Column 5
                all_is_known_gap_size = adjacency_data["is_known_gap_size"]
                comp_type = "N" if any(all_is_known_gap_size) else "U"

                # AGP column 6b
                gap_size = 100
                all_gap_sizes = adjacency_data["gap_size"]
                fltrd_gap_sizes = [all_gap_sizes[i] for i in range(len(all_gap_sizes)) if all_is_known_gap_size[i]]
                if fltrd_gap_sizes:
                    if len(fltrd_gap_sizes) == 1:
                        gap_size = fltrd_gap_sizes[0]
                    else:
                        gap_size = get_gap_size(fltrd_gap_sizes, gap_func)

                # AGP column 7b
                all_gap_types = set(adjacency_data["gap_type"])
                gap_type = "scaffold"
                if len(all_gap_types) == 1:
                    gap_type = all_gap_types.pop()

                # AGP column 8b
                has_linkage = "yes" if any(adjacency_data["linkage"]) else "no"

                # AGP column 9b
                all_evidences = set(adjacency_data["linkage_evidence"])
                linkage_evidence = "na"
                if has_linkage == "yes":
                    if "na" in all_evidences:
                        all_evidences.remove("na")
                    linkage_evidence = ";".join([str(i) for i in all_evidences])

                agp.add_gap_line(
                    obj_header,
                    str(obj_pos + 1),
                    str(obj_pos + gap_size),
                    str(obj_id),
                    comp_type,
                    str(gap_size),
                    gap_type,
                    has_linkage,
                    linkage_evidence
                )
                obj_pos += gap_size

            obj_id += 1
            visited_nodes.add(next_node)
            current_node = next_node

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
    parser = argparse.ArgumentParser(description="Merging and reconciling scaffolds", usage="ragtag.py merge <asm.fa> <scf1.agp> <scf2.agp> [...]")
    parser.add_argument("components", metavar="<asm.fasta>", nargs='?', default="", type=str, help="assembly fasta file (uncompressed or bgzipped)")
    parser.add_argument("agps", metavar="<scf1.agp> <scf2.agp> [...]", nargs='*', default=[], type=str, help="scaffolding AGP files")

    merge_options = parser.add_argument_group("merging options")
    merge_options.add_argument("-f", metavar="FILE", default="", type=str, help="CSV list of (AGP file,weight) [null]")
    merge_options.add_argument("-j", metavar="<skip.txt>", type=str, default="", help="list of query headers to leave unplaced [null]")
    merge_options.add_argument("-l", metavar="INT", default=100000, type=int, help="minimum assembly sequence length [100000]")
    merge_options.add_argument("-e", metavar="FLOAT", default=0.0, type=float, help="minimum edge weight. NA if using Hi-C [0.0]")
    merge_options.add_argument("--gap-func", metavar="STR", default="min", type=str, help="function for merging gap lengths {'min', 'max', or 'mean'} [min]")

    io_options = parser.add_argument_group("input/output options")
    io_options.add_argument("-o", metavar="PATH", type=str, default="ragtag_output", help="output directory [./ragtag_output]")
    io_options.add_argument("-w", action='store_true', default=False, help="overwrite intermediate files")
    io_options.add_argument("-u", action='store_true', default=False, help="add suffix to unplaced sequence headers")
    io_options.add_argument("--debug", action='store_true', default=False, help=argparse.SUPPRESS)

    hic_options = parser.add_argument_group("Hi-C options")
    hic_options.add_argument("-b", metavar="FILE", default="", type=str, help="Hi-C alignments in BAM format, sorted by read name [null]")
    hic_options.add_argument("-r", metavar="STR", default="GATC", type=str, help="CSV list of restriction enzymes/sites or 'DNase' [GATC]")
    hic_options.add_argument("-p", metavar="FLOAT", default=1.0, type=float, help="portion of the sequence termini to consider for links [1.0]")
    hic_options.add_argument("--list-enzymes", action='store_true', default=False, help="list all available restriction enzymes/sites")

    args = parser.parse_args()

    # Print a restriction enzyme help message if requested
    if args.list_enzymes:
        RestrictionEnzymes.get_info()
        sys.exit()

    if not args.components:
        parser.print_help()
        print("\n** The assembly FASTA file is required **")
        sys.exit()

    if not args.agps and not args.f:
        parser.print_help()
        print("\n** At least two AGP files are required **")
        sys.exit()

    log("RagTag " + get_ragtag_version())
    log("This is a beta version of `ragtag merge`")
    log("CMD: ragtag.py merge " + " ".join(sys.argv[1:]))

    # Check that the components FASTA file exists
    comp_fname = args.components
    if not os.path.isfile(comp_fname):
        raise ValueError("Could not find file: %s" % comp_fname)

    # Optional arguments
    agp_fofn = args.f
    hic_bam_fname = args.b
    re_string = args.r
    portion = args.p

    # Set the minimum component sequence length
    min_comp_len = args.l
    if min_comp_len < 0:
        min_comp_len = 0

    # Set the minimum edge weight
    min_edge_weight = args.e
    if min_edge_weight < 0:
        min_edge_weight = 0

    # Set the gap merging function options
    gap_func = args.gap_func.upper()
    if gap_func not in {"MIN", "MAX", "MEAN"}:
        raise ValueError("Gap merging function must be either 'min', 'max', or 'mean'. Got: {}".format(args.gap_func))

    # Debugging options
    debug_mode = args.debug

    # I/O options
    output_path = args.o
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = os.path.abspath(output_path) + "/"

    overwrite_files = args.w
    add_suffix = args.u
    if not add_suffix:
        log(
            "WARNING: Without '-u' invoked, some component/object AGP pairs might share the same ID. Some external programs/databases don't like this. To ensure valid AGP format, use '-u'.")

    # get the set of contigs to skip
    comp_exclusion_set = set()
    skip_fname = args.j
    if skip_fname:
        skip_fname = os.path.abspath(skip_fname)
        with open(skip_fname, "r") as f:
            for line in f:
                comp_exclusion_set.add(line.rstrip().split()[0])

    # Setup a file for general logging
    merge_log = output_path + "ragtag.merge.err"
    open(merge_log, "w").close()  # Wipe the log file

    # Process the AGP files
    agp_list = [os.path.abspath(i) for i in args.agps]
    weight_list = [1 for _ in range(len(agp_list))]

    # Check for file of AGPs and weights
    if agp_fofn:
        agp_list, weight_list = [], []
        with open(agp_fofn, "r") as f:
            for line in f:
                fields = line.rstrip().split(",")
                agp_list.append(fields[0])
                weight_list.append(float(fields[1]))

    if len(agp_list) < 2:
        raise ValueError("At least two AGP files are required for merging")

    # Build the graph and filter nodes by sequence length
    log("Building the scaffold graph from the AGP files")
    agp_multi_sg = AGPMultiScaffoldGraph(comp_fname)
    agp_multi_sg.add_agps(agp_list, in_weights=weight_list, exclusion_set=comp_exclusion_set)
    if min_comp_len:
        agp_multi_sg.filter_by_seq_len(min_comp_len)
    if debug_mode:
        agp_multi_sg.connect_and_write_gml(output_path + "ragtag.merge.msg.gml")

    # Merge the SAG
    log("Merging the scaffold graph")
    agp_sg = agp_multi_sg.merge()

    # Check if we are using Hi-C links to weight the graph.
    if hic_bam_fname:
        log("Weighting the scaffold graph with Hi-C links")
        if not comp_fname or not re_string:
            raise RuntimeError("Hi-C requires alignments (-b) assembly sequences (-a) and restriction sites (-r)")

        cmd = [
            "ragtag_create_links.py",
            "-a", comp_fname,
            "-b", hic_bam_fname,
            "-r", re_string,
            "-p", str(portion)
        ]

        out_links_fname = output_path + "ragtag.merge.links"
        if os.path.isfile(out_links_fname):
            if not overwrite_files:
                log("Retaining pre-existing file: " + out_links_fname)
            else:
                run_oae(cmd, output_path + "ragtag.merge.links", merge_log)
        else:
            run_oae(cmd, output_path + "ragtag.merge.links", merge_log)

        hic_sg = build_hic_graph(out_links_fname, comp_fname)
        agp_sg = agp_sg.steal_weights_from(hic_sg)

    # Filter by edge weight
    if min_edge_weight and not hic_bam_fname:
        agp_sg.filter_by_weight(min_edge_weight)

    if debug_mode:
        agp_sg.connect_and_write_gml(output_path + "ragtag.merge.sg.gml")

    # Compute a solution to the ScaffoldGraph
    log("Computing a scaffolding solution")
    cover_graph = get_maximal_matching(agp_sg)
    if debug_mode:
        nx.readwrite.gml.write_gml(cover_graph, output_path + "ragtag.merge.covergraph.gml")

    # Write the scaffolding output to an AGP file
    log("Writing results")
    write_agp_solution(cover_graph, agp_sg, output_path + "ragtag.merge.agp", gap_func=gap_func, add_suffix_to_unplaced=add_suffix)

    # Generate a FASTA file corresponding to the AGP
    cmd = [
        "ragtag_agp2fasta.py",
        output_path + "ragtag.merge.agp",
        comp_fname
    ]
    run_oae(cmd, output_path + "ragtag.merge.fasta", merge_log)


if __name__ == "__main__":
    main()