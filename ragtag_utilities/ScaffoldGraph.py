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
import abc
import itertools
import numbers

import pysam
import numpy as np
import networkx as nx

from ragtag_utilities.AGPFile import AGPFile
from ragtag_utilities.utilities import get_ragtag_version


"""
ScaffoldGraph.py defines graph data structures and related objects for assembly scaffolding. These objects are inspired by
concepts presented in CAMSA (https://doi.org/10.1186/s12859-017-1919-y) and
SALSA (https://doi.org/10.1186/s12864-017-3879-z). Specifically, each of these tools uses the "Scaffold Graph" (SG), a
data structure for representing assembly scaffolding adjacencies.

ScaffoldGraphBase is an abstract base class defining the core of the SG. ScaffoldGraphBase is a wrapper around a
networkx graph that enforces certain rules for adding nodes and edges. Nodes must be strings that end in either "_b" or
"_e", representing the beginning or end of a contig, respectively. Additionally, edges must define a weight and other
metadata necessary for AGP output. Inheriting classes can define additional mandatory edge data. ScaffoldGraphBase also
stores information about the sequence components that the graph nodes refer to. Nodes can only be added if they
correspond to sequences found in the components FASTA file. 

Inheriting subclasses can be graphs or mutligraphs, so methods/attributes defined in the base class should be
appropriate for either. Any methods/attributes that are specific to multi or non-multigraphs should be defined in
subclasses.

MultiScaffoldGraph inherits from ScaffoldGraphBase but overrides the graph to be a multigraph. This should be used
scaffolding applications where adjacencies are supported by multiple pieces of evidence, like multiple AGP files
or multiple ultra long reads. The MultiScaffoldGraph defines a merge function which coverts the multi-graph
to a graph, storing multi-edge information appropriately (while still defining a single edge weight).

ScaffoldGraph inherits from ScaffoldGraphBase and adds additional functionality that would not be appropriate for
multi-graphs. This is where utilities for computing scaffolding solutions are defined.

AGPMultiScaffoldGraph inherits from MultiScaffoldGraph to define AGP-specific functionality. It draws upon the
AssemblyPoint object to build the graph, even for unkown scaffolding orientations. 


Common names and abbreviations:

  "SG" ---- Scaffold Graph
  "AP" ---- Assembly Point
  "AGP" --- Refers to AGP files (https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/)
  "GML" --- Graph Modeling Language 
"""


class AssemblyPoint:
    """
    An AssemblyPoint represents a sequence adjacency. This adjacency is between two seqeucnes, each with an orientation.
    Furthermore, the adjacency has some more attributes, like a confidence value. Other adjacency attributes can be
    defined in an AGP file. The object stores the source AGP file yielding the adjacency, as well as the following
    AGP info (assuming a gap separates the sequences):

      1. Is the gap of a known or unknown size?
      2. How long is the gap?
      3. What kind of gap is it (e.g. scaffold, centromere, etc.)?
      4. What evidence was used to justify the sequence adjacency.

    If adjacent sequences don't have a gap in between them, these attributes can store appropriate values as well. For
    example, the gap_size would be 0 and the evidence would be "bookend".

    The AssemblyPoint defines "realizations", or all possible edges that can be formed between sequence nodes of a SAG.
    """

    def __init__(self, seq1, strand1, seq2, strand2, confidence, source_agp_fname, is_known_gap_size, gap_size, gap_type, linkage, linkage_evidence):
        """
        :param seq1: The first sequence of the adjacency.
        :param strand1: The orientation of seq1.
        :param seq2: The second sequence of the adjacency.
        :param strand2: The orientation of seq2.
        :param confidence: A confidence value for this adjacency.
        :param source_agp_fname: The file name of the AGP file that yields this adjacency.
        :param is_known_gap_size: Is the gap size between seq1 and seq2 known?
        :param gap_size: The size of the gap between seq1 and seq2.
        :param gap_type: The type of gap between seq1 and seq2.
        :param linkage: Is there evidence of linkage?
        :param linkage_evidence: The type of evidence used to support this adjacency.
        """
        self.seq1 = seq1
        self.strand1 = strand1
        self.seq2 = seq2
        self.strand2 = strand2
        self.confidence = confidence
        self.source_agp_fname = source_agp_fname
        self.is_known_gap_size = is_known_gap_size  # "N" vs "U" AGP gap component types (AGP column 5)
        self.gap_size = gap_size  # AGP column 6b
        self.gap_type = gap_type  # AGP column 7b
        self.linkage = linkage  # AGP column 8b
        self.linkage_evidence = linkage_evidence  # AGP column 9b

        # Replace all ambiguous AGP orientation values with the '?' value.
        self.strand1 = "?" if self.strand1 in {"na", "0"} else self.strand1
        self.strand2 = "?" if self.strand2 in {"na", "0"} else self.strand2

    def get_nodes(self):
        return [self.seq1 + "_b", self.seq1 + "_e", self.seq2 + "_b", self.seq2 + "_e"]

    def get_realizations(self, candidate_scaling_factor=0.75):
        """
        :return: All possible SAG assembly edges.
        """
        seq1_choices = {"+", "-"} if self.strand1 == "?" else {self.strand1}
        seq2_choices = {"+", "-"} if self.strand2 == "?" else {self.strand2}
        is_candidate = 1 if self.strand1 == "?" or self.strand2 == "?" else 0

        for strand1, strand2 in itertools.product(seq1_choices, seq2_choices):
            start = self.seq1 + ("_e" if strand1 == "+" else "_b")
            end = self.seq2 + ("_b" if strand2 == "+" else "_e")
            start, end = (start, end) if start < end else (end, start)

            yield start, end, (1 - ((1-candidate_scaling_factor) * is_candidate)) * self.confidence


class ScaffoldGraphBase:
    """
    Base class for the Scaffold Graph (SG). SGBase is a wrapper around a generic networkx graph.
    This class defines/enforces two key feautres of the Scaffold Graph:
        1) Nodes must end with either "_b" or "_e" indicating the beginning or end of a contig, respectively.
        2) All edges must have a weight

    Since MultiSG and SG will both inherit from this base class, SGBase should not define any methods that don't
    make sense for both multi and non multigraphs. Such methods can be defined directly in MultiSG and SG.
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, components_fasta_fname):
        self.graph = nx.Graph()
        self._req_edge_data = {
            "weight",
            "source_fname",
            "is_known_gap_size",
            "gap_size",
            "gap_type",
            "linkage",
            "linkage_evidence"
        }

        self.components_fasta_fname = components_fasta_fname
        self.component_lens = dict()  # component name -> length
        self._set_component_lens()

    def __getitem__(self, e):
        return self.graph[e]

    def _set_component_lens(self):
        fai = pysam.FastaFile(self.components_fasta_fname)

        for i in fai.references:
            self.component_lens[i] = fai.get_reference_length(i)

        fai.close()

    def get_component_len(self, c):
        return self.component_lens[c]

    @property
    def components(self):
        return set(self.component_lens.keys())

    @property
    def degree(self):
        return self.graph.degree

    @property
    def nodes(self, **kwargs):
        return self.graph.nodes(**kwargs)

    @property
    def edges(self):
        return self.graph.edges

    def add_node(self, n, **kwargs):
        """ Add a node to the graph. """
        if not isinstance(n, str):
            raise TypeError("Nodes must be strings. Got {}".format(type(n)))

        if not n.endswith("_b") and not n.endswith("_e"):
            raise ValueError("Nodes must end with either '_b' or '_e'.")

        if n[:-2] not in self.component_lens:
            raise ValueError("{} was not found in {}".format(n[:-2], self.components_fasta_fname))

        self.graph.add_node(n, **kwargs)

    def remove_node(self, n):
        self.graph.remove_node(n)

    def neighbors(self, n):
        return self.graph.neighbors(n)

    def add_edge(self, u, v, **kwargs):
        """ Add edges like normal for a networkx graph, but ensure that the edge has necessary data. """
        if not kwargs:
            raise ValueError("Scaffold Graph edges must have the following data: {}".format(",".join(list(self._req_edge_data))))

        for i in self._req_edge_data:
            if i not in kwargs:
                raise ValueError("Scaffold Graph edges must have the following data: {}".format(i))

        # Add the nodes before making the edge.
        # Not strictly necessary for networkx, but it ensures that all nodes match the fasta file
        if u not in self.nodes:
            self.add_node(u)

        if v not in self.nodes:
            self.add_node(v)

        self.graph.add_edge(u, v, **kwargs)

    def add_edges_from(self, ebunch, **kwargs):
        # Check that each edge has the required edge data
        for edge in ebunch:
            if len(edge) < 3:
                raise ValueError("Edges must have edge data.")

            for i in self._req_edge_data:
                if i not in edge[2]:
                    raise ValueError("Edges must have the following property: {}".format(i))

            # Add edges one at a time so that nodes can be checked against the fasta file
            for i in edge[2]:
                kwargs[i] = edge[2][i]
            self.graph.add_edge(edge[0], edge[1], **kwargs)

    def remove_edge(self, u, v):
        self.graph.remove_edge(u, v)

    def has_edge(self, u, v):
        return self.graph.has_edge(u, v)

    def write_gml(self, f):
        """
        Write the graph in GML format.
        :param f: Write the graph to this file handle.
        """
        G = self.graph.copy()

        # networkx doesn't like writing non-string attributes to GML
        for u, v in G.edges:
            for key in G[u][v]:
                G[u][v][key] = str(G[u][v][key])
        nx.readwrite.gml.write_gml(G, f)

    def connect_and_write_gml(self, f):
        """
        Connect scaffold edges, then write the graph in GML format
        :param f: Write the graph to this file handle.
        """
        G = self.graph.copy()
        node_base_set = set([i[:-2] for i in list(G.nodes)])
        for node in node_base_set:
            G.add_edge(node + "_b", node + "_e")

        # networkx doesn't like writing non-string attributes to GML
        for u, v in G.edges:
            for key in G[u][v]:
                G[u][v][key] = str(G[u][v][key])
        nx.readwrite.gml.write_gml(G, f)

    def filter_by_seq_len(self, min_len):
        """
        Remove nodes whose associated component sequence is shorter than min_len.
        :param min_len: minimum component sequence length
        """
        if not isinstance(min_len, numbers.Number):
            raise TypeError("min_len must be a number")

        # Iterate over the nodes and remove any nodes shorter than min_len
        old_nodes = set(self.nodes)
        for n in old_nodes:
            comp_name = n[:-2]
            if self.get_component_len(comp_name) < min_len:
                self.remove_node(n)


class ScaffoldGraph(ScaffoldGraphBase):
    """
    The Scaffold Graph object extends the baseclass to include additional utilities.
    """

    def __init__(self, components_fasta_fname):
        super().__init__(components_fasta_fname)

    def filter_by_weight(self, w):
        """
        Remove edges with weight less than w.
        :param w: The minimum edge weight
        """
        G = nx.Graph()

        for u, v in self.edges:
            if self.graph[u][v]["weight"] >= w:
                # Add the nodes first in case they have data
                G.add_node(u, **self.nodes(data=True)[u])
                G.add_node(v, **self.nodes(data=True)[v])
                G.add_edge(u, v, **self.graph[u][v])

        self.graph = G

    def filter_non_one(self):
        """ Remove edges where the weight is not equal to one. """
        G = nx.Graph()

        for u, v in self.edges:
            if self.graph[u][v]["weight"] == 1:
                # Add the nodes first in case they have data
                G.add_node(u, **self.nodes(data=True)[u])
                G.add_node(v, **self.nodes(data=True)[v])
                G.add_edge(u, v, **self.graph[u][v])

        self.graph = G

    def get_max_weight_matching(self):
        return nx.max_weight_matching(G=self.graph)

    def get_max_incident_edge_weight(self, u, v):
        """
        Given two nodes, find the incident (to either node) edge with the maximum weight.

        Code inspired by SALSA2: https://github.com/marbl/SALSA/blob/master/fast_scaled_scores.py
        :param u:
        :param v:
        :return:
        """
        max_weight = 0
        for node in self.graph.neighbors(u):
            if node != v and self.graph[u][node]["weight"] > max_weight:
                max_weight = self.graph[u][node]["weight"]

        for node in self.graph.neighbors(v):
            if node != u and self.graph[v][node]["weight"] > max_weight:
                max_weight = self.graph[v][node]["weight"]

        return max_weight

    def best_buddy_scale(self):
        """
        Scale weights according to the SALSA2 best buddy scheme:
        https://doi.org/10.1371/journal.pcbi.1007273

        Code inspired by SALSA2: https://github.com/marbl/SALSA/blob/master/fast_scaled_scores.py

        Returns a new instance of the class.
        """
        max_inc_weights = dict()  # node -> maximum incident edge weight

        # Iterate through the edges and get the max weight for each node
        for u, v in self.edges:
            if u not in max_inc_weights:
                max_inc_weights[u] = self.graph[u][v]["weight"]
            else:
                if self.graph[u][v]["weight"] > max_inc_weights[u]:
                    max_inc_weights[u] = self.graph[u][v]["weight"]

            if v not in max_inc_weights:
                max_inc_weights[v] = self.graph[u][v]["weight"]
            else:
                if self.graph[u][v]["weight"] > max_inc_weights[v]:
                    max_inc_weights[v] = self.graph[u][v]["weight"]

        # Iterate through the edges again and scale the weights
        G = ScaffoldGraph(self.components_fasta_fname)
        for u, v in self.edges:
            w = self.graph[u][v]["weight"]
            max_u = max_inc_weights[u]
            max_v = max_inc_weights[v]
            best_alt = max(max_u, max_v)

            if best_alt == w:
                best_alt = self.get_max_incident_edge_weight(u, v)

            if best_alt == 0:
                best_alt = 1

            # Make a copy of the edge data and remove the weight
            new_edge_data = dict(self.graph[u][v])
            new_edge_data.pop("weight")
            G.add_edge(u, v, weight=w/best_alt, **new_edge_data)

        return G

    def steal_weights_from(self, in_sg):
        """ Replace weights with matching weights from a new scaffold graph. All other edge data remains. """
        G = ScaffoldGraph(self.components_fasta_fname)
        for u, v in self.edges:
            if in_sg.has_edge(u, v):
                # Make a copy of the edge data and remove the weight
                new_edge_data = dict(self.graph[u][v])
                new_edge_data.pop("weight")
                G.add_edge(u, v, weight=in_sg[u][v]["weight"], **new_edge_data)

        return G


class MultiScaffoldGraph(ScaffoldGraphBase):
    """
    The Multi Scaffold Graph object is distinct in that it is a multigraph rather than a graph. This is convenient
    when multiple pieces of evidence can be used to support sequence joins. Each piece of evidence can be stored as
    a distinct edge between nodes. Scaffolding solutions are likely to be computed from regular graphs, MultiSG
    defines a merge function that combines data from multi edges and returns a instance of SG.
    """

    def __init__(self, components_fasta_fname):
        super().__init__(components_fasta_fname)
        self.graph = nx.MultiGraph()  # Override graph to be a multigraph

    @staticmethod
    def _merge_edge_dicts(*edge_dicts):
        """
        Given a list of edge data, combine like data into lists.
        Sum the weight list and provide that as the new weight.
        :param edge_dicts: list of edge data dictionaries
        :return: merged edge data dictionary
        """
        data = dict()
        for k in edge_dicts[0]:
            data[k] = [edge_dicts[0][k]]

        for new_dict in edge_dicts[1:]:
            for k in new_dict:
                if k not in data:
                    raise ValueError("Inconsistent edge data.")

                data[k].append(new_dict[k])

        data["original_weights"] = data["weight"]
        data["weight"] = sum(data["weight"])
        return data

    def merge(self):
        """ Combine edges connecting the same nodes. Return a new SG graph. """
        G = ScaffoldGraph(self.components_fasta_fname)

        # TODO implement an add_nodes_from wrapper so direct access is unnecessary
        G.graph.add_nodes_from(self.nodes(data=True))
        for u, v, c in self.graph.edges:
            edge_data_list = [self.graph[u][v][i] for i in self.graph[u][v]]
            G.add_edge(u, v, **self._merge_edge_dicts(*edge_data_list))

        return G


class AGPMultiScaffoldGraph(MultiScaffoldGraph):
    """
    An extension of the MultiSG that is constructed from multiple AGP files.

    https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
    """

    def __init__(self, components_fasta_fname):
        super().__init__(components_fasta_fname)

        # Store agp files and associated weights as they are added
        self.agps = []
        self.weights = []

    def _get_assembly_points(self, agp, weight):
        """
        Find all adjacencies defined in an AGP file
        :param agp: An AGP file defining sequence adjacencies
        :param weight: The weight to assign to each adjacency
        """
        comps = set()
        prev_obj = ""
        seq1 = ""
        strand1 = ""

        # Gap info
        gap_count = 0
        prev_agp_known = None
        prev_gap_size = 0
        prev_gap_type = None
        prev_linkage = ""
        prev_evidence = ""

        # Iterate over the AGP file and yield assembly points
        agp_file = AGPFile(agp)
        for agp_line in agp_file.iterate_lines():
            if not agp_line.is_gap:
                # Add this component to our master list
                if agp_line.comp not in self.component_lens:
                    raise RuntimeError("{} is in {} but not {}.".format(agp_line.comp, agp, self.components_fasta_fname))
                comps.add(agp_line.comp)

                comp_len = agp_line.comp_end
                if comp_len < self.get_component_len(agp_line.comp):
                    raise RuntimeError("only complete components can be added to the graph.")

                if comp_len > self.get_component_len(agp_line.comp):
                    raise RuntimeError("inconsistent component lengths: {} bp in {} and {} bp in {}". format(comp_len, agp, self.get_component_len(agp_line.comp), self.components_fasta_fname))

                if agp_line.obj == prev_obj:
                    # Check if these components are bookended (no gap in between)
                    if not gap_count:
                        prev_evidence = "bookend"

                    # Check if two consecutive gaps preceded this component
                    if gap_count > 1:
                        raise ValueError("Consecutive gaps in the AGP file are not currently supported.")

                    yield AssemblyPoint(
                        seq1,
                        strand1,
                        agp_line.comp,
                        agp_line.orientation,
                        weight,
                        agp,
                        prev_agp_known,
                        prev_gap_size,
                        prev_gap_type,
                        prev_linkage,
                        prev_evidence
                    )

                    # Set this component as the previous component
                    seq1 = agp_line.comp
                    strand1 = agp_line.orientation

                    gap_count = 0
                    prev_agp_known = None
                    prev_gap_size = 0
                    prev_gap_type = None
                    prev_linkage = ""
                    prev_evidence = ""

                else:
                    seq1 = agp_line.comp
                    strand1 = agp_line.orientation

                    prev_obj = agp_line.obj
                    gap_count = 0
                    prev_agp_known = None
                    prev_gap_size = 0
                    prev_gap_type = None
                    prev_linkage = ""
                    prev_evidence = ""
            else:
                if agp_line.obj == prev_obj:
                    gap_count += 1
                    prev_agp_known = True if agp_line.comp_type == "N" else False
                    prev_gap_size = agp_line.gap_len
                    prev_gap_type = agp_line.gap_type
                    prev_linkage = True if agp_line.linkage == "yes" else False
                    prev_evidence = agp_line.linkage_evidence

        if comps != self.components:
            raise ValueError("Input AGPs do not have the same set of components.")

    def add_agps(self, in_agps, in_weights=None, exclusion_set=None):
        """
        Update the graph given AGP files.
        :param in_agps: A list of AGP files
        :param in_weights: A list of weights for the AGP files
        :param exclusion_set: set of components to skip
        """
        if exclusion_set is None:
            exclusion_set = set()

        if in_weights is None:
            in_weights = [1] * len(in_agps)
        else:
            if len(in_agps) != len(in_weights):
                raise ValueError("Must assign one weight per AGP file, or none at all.")

        # Iterate through the AGP files to update the graph
        for agp, weight in zip(in_agps, in_weights):
            agp = os.path.abspath(agp)
            if agp in self.agps:
                raise ValueError("%s appears more than once." % agp)

            # Add this AGP and its weight to the master list
            self.agps.append(agp)
            self.weights.append(weight)

            for ap in self._get_assembly_points(agp, weight):
                for u, v, w in ap.get_realizations():
                    if u[:-2] not in exclusion_set and v[:-2] not in exclusion_set:
                        u_base, v_base = u[:-2], v[:-2]
                        u_len = self.get_component_len(u_base)
                        self.add_edge(
                            u,
                            v,
                            weight=w,
                            source_fname=ap.source_agp_fname,
                            is_known_gap_size=ap.is_known_gap_size,
                            gap_size=ap.gap_size,
                            gap_type=ap.gap_type,
                            linkage=ap.linkage,
                            linkage_evidence=ap.linkage_evidence,
                            seqs=[u_base, v_base],
                            pos=[u_len, 0]
                        )


class Alignment:
    """
    An object representing query a sequence/gap joining two target sequences.

      target_from --------- target "from" scaffold graph node
      target_to ----------- target "to" scaffold graph node
      query --------------- name of the sequence used to join "from" and "to". ignored if gap
      query_len ----------- length of the query sequence
      my_query_end -------- the 0-based position in the query that corresponds to the end of "from"
      their_query_start --- the 0-based position in the query that corresponds to the beginning of "to"
      strand -------------- the strand of the query sequence
      is_gap -------------- boolean indicating a gap or a sequence
    """

    def __init__(self, target_from, target_to, query, query_len, my_query_end, their_query_start, strand, is_gap=False):
        """

        :param target_from:
        :param target_to:
        :param query:
        :param query_len:
        :param my_query_end:
        :param their_query_start:
        :param strand: 0, 1 for forward, reverse
        :param is_gap:
        """
        self.target_from = target_from
        self.target_to = target_to
        self.query = query
        self.query_len = query_len
        self.my_query_end = my_query_end  # 0-based query start coordinate
        self.their_query_start = their_query_start  # 0-based query end coordinate
        self.strand = strand  # 0, 1 for forward, reverse

        self.is_gap = is_gap
        self.query_span = self.their_query_start - self.my_query_end

    def __str__(self):
        return "::".join([
            self.target_from,
            self.target_to,
            self.query
        ])

    def reverse(self):
        return Alignment(
            self.target_to,
            self.target_from,
            self.query,
            self.query_len,
            self.query_len - self.their_query_start,
            self.query_len - self.my_query_end,
            1 - self.strand,
            is_gap=self.is_gap
        )


class PatchScaffoldGraph:
    """
    A directed scaffold graph useful for continuous scaffolding

    Terms:
      Target assembly - a genome assembly that we wish to scaffold
      Target sequence - an individual sequence in a target assembly
      Query sequence - a sequence, independent from the target assembly, used to scaffold the target assembly
      Adjacency - a pair of termini between two distinct target sequences suggested to be joined

    The object keeps track of a directed graph. Nodes in the graph are target sequence termini
    (target sequence + "_b" or "_e" representing the beginning or end of the sequence, respectively). Edges
    connecting these nodes represent an adjacency, and the necessary information must be encoded in an Alignment
    object.

    This class defines a method to produce an AGP file from the scaffold graph, assuming that the edges
    form a matching.
    """

    def __init__(self, components_fn):
        self.graph = nx.DiGraph()

        self.components_fn = components_fn  # name of FASTA file with all relevant target/query sequences
        self.component_lens = dict()  # component name -> length
        self._set_component_lens()

    def __getitem__(self, e):
        return self.graph[e]

    def _set_component_lens(self):
        fai = pysam.FastaFile(self.components_fn)

        for i in fai.references:
            self.component_lens[i] = fai.get_reference_length(i)

    @property
    def nodes(self, **kwargs):
        return self.graph.nodes(**kwargs)

    @property
    def edges(self):
        return self.graph.edges

    def has_edge(self, u, v):
        return self.graph.has_edge(u, v)

    def neighbors(self, n):
        return self.graph.neighbors(n)

    def add_edge(self, u, v, aln):
        """
        :param u: from
        :param v: to
        :param aln: Alignment object
        """
        if not isinstance(aln, Alignment):
            raise TypeError("aln must be an Alignment object.")

        if u != aln.target_from or v != aln.target_to:
            raise RuntimeError("Scaffold graph nodes don't match Alignment nodes.")

        if self.has_edge(u, v):
            self.graph[u][v]["weight"] += 1
            self.graph[u][v]["alignment"] = aln

            self.graph[v][u]["weight"] += 1
            self.graph[v][u]["alignment"] = aln.reverse()
        else:
            self.graph.add_edge(u, v, weight=1, alignment=aln)
            self.graph.add_edge(v, u, weight=1, alignment=aln.reverse())

    def remove_edge(self, u, v):
        self.graph.remove_edge(u, v)

    def write_gml(self, f):
        G = nx.DiGraph()
        for u, v in self.edges:
            G.add_edge(
                u,
                v,
                target_from=self.graph[u][v]["alignment"].target_from,
                target_to=self.graph[u][v]["alignment"].target_to,
                query=self.graph[u][v]["alignment"].query,
                query_len=str(self.graph[u][v]["alignment"].query_len),
                my_query_end=str(self.graph[u][v]["alignment"].my_query_end),
                their_query_start=str(self.graph[u][v]["alignment"].their_query_start),
                strand=str(self.graph[u][v]["alignment"].strand),
                is_gap=str(self.graph[u][v]["alignment"].is_gap)
            )

            G.add_edge(
                v,
                u,
                target_from=self.graph[v][u]["alignment"].target_from,
                target_to=self.graph[v][u]["alignment"].target_to,
                query=self.graph[v][u]["alignment"].query,
                query_len=str(self.graph[v][u]["alignment"].query_len),
                my_query_end=str(self.graph[v][u]["alignment"].my_query_end),
                their_query_start=str(self.graph[v][u]["alignment"].their_query_start),
                strand=str(self.graph[v][u]["alignment"].strand),
                is_gap=str(self.graph[v][u]["alignment"].is_gap)
            )
        nx.readwrite.gml.write_gml(G, f)

    def remove_heavier_than(self, w):
        """
        Remove edges with weight > w
        :param w: maximum edge weight
        """
        G = nx.DiGraph()
        for u, v in self.edges:
            if self.graph[u][v]["weight"] <= w:
                G.add_edge(u, v, weight=self.graph[u][v]["weight"], alignment=self.graph[u][v]["alignment"])
        self.graph = G

    def max_weight_matching(self):
        """
        Return a new PatchScaffoldGraph that is a matching
        :return: Matching PatchScaffoldGraph
        """
        undirected_sg = nx.Graph()
        for u, v in self.edges:
            undirected_sg.add_edge(u, v, weight=self.graph[u][v]["weight"])
        matching = nx.max_weight_matching(G=undirected_sg)

        # Remove any cycles that would result from adding intra-sequence edges
        # Make a new graph
        cover_graph = nx.Graph()
        cover_graph.add_nodes_from(undirected_sg.nodes(data=True))

        # Add edges connecting contig ends
        node_base_set = set([i[:-2] for i in list(cover_graph.nodes)])
        for node in node_base_set:
            cover_graph.add_edge(node + "_b", node + "_e", weight=np.inf)

        # Add the edges that form the matching
        for u, v in matching:
            cover_graph.add_edge(u, v, **undirected_sg[u][v])

        # Remove any potential cycles
        edges_to_delete = []
        for cc in nx.connected_components(G=cover_graph):
            cc = cover_graph.subgraph(cc).copy()
            if cc.number_of_nodes() == cc.number_of_edges():
                assembly_edges = cc.edges(data=True)
                edge_to_delete = min(assembly_edges, key=lambda entry: entry[2]["weight"])
                edges_to_delete.append((edge_to_delete[0], edge_to_delete[1]))

        # Add the edges that form the matching into a new PatchScaffoldGraph
        new_psg = PatchScaffoldGraph(self.components_fn)
        for u, v in matching:
            if (u, v) not in edges_to_delete and (v, u) not in edges_to_delete:
                new_psg.add_edge(u, v, self.graph[u][v]["alignment"])
                new_psg.add_edge(v, u, self.graph[v][u]["alignment"])

        return new_psg

    def write_agp(self, agp_fn, ref_fn, add_suffix_to_unplaced=False):
        """
        Write the AGP file implied by the scaffold graph
        :param agp_fn: AGP file name
        :param ref_fn: reference FASTA file name
        """
        used_components = set()
        used_edges = set()
        obj_header_idx = -1

        agp = AGPFile(agp_fn, "w")
        agp.add_pragma()
        agp.add_comment("# AGP created by RagTag {}".format(get_ragtag_version()))

        while True:
            # Find a starting node
            from_node = None
            to_node = None
            cur_ref = None
            for u, v in sorted(self.edges):
                if (u, v) not in used_edges:
                    u_base = u[:-2]

                    u_degree = 0
                    if u_base + "_b" in self.nodes:
                        u_degree += self.graph.degree[u_base + "_b"]
                    if u_base + "_e" in self.nodes:
                        u_degree += self.graph.degree[u_base + "_e"]

                    assert u_degree in {2, 4}

                    # Check if we have found a starting target sequence
                    if u_degree == 2:
                        cur_ref = u_base
                        from_node = u
                        to_node = v
                        used_edges.add((u, v))
                        used_edges.add((v, u))
                        break

            # If we haven't found a new starting target sequence, we are done
            if from_node is None:
                break

            # Initialize this object
            obj_header_idx += 1
            obj_header = "scf" + "{0:08}".format(obj_header_idx)
            obj_pos = 0
            obj_pid = 1

            # Process the first target sequence
            cur_ref_len = self.component_lens[cur_ref]
            cur_ref_strand = "+"
            if from_node.endswith("_b"):
                cur_ref_strand = "-"
            agp.add_seq_line(obj_header, obj_pos+1, obj_pos+cur_ref_len, obj_pid, "W", cur_ref, 1, cur_ref_len, cur_ref_strand)
            obj_pos += cur_ref_len
            obj_pid += 1
            used_components.add(cur_ref)

            # Process the remaining sequences.
            next_edge_exists = True
            while next_edge_exists:
                # Process the patch
                patch_aln = self.graph[from_node][to_node]["alignment"]
                patch_query = patch_aln.query
                patch_strand = "+"
                if patch_aln.strand:
                    patch_strand = "-"

                patch_len = patch_aln.their_query_start - patch_aln.my_query_end
                if patch_len > 0:
                    if patch_aln.is_gap:
                        agp.add_gap_line(obj_header, obj_pos+1, obj_pos+patch_len, obj_pid, "N", patch_len, "scaffold", "yes", "align_genus")
                    else:
                        agp.add_seq_line(obj_header, obj_pos+1, obj_pos+patch_len, obj_pid, "W", patch_query, patch_aln.my_query_end+1, patch_aln.their_query_start, patch_strand)
                        used_components.add(patch_query)
                    obj_pos += patch_len
                    obj_pid += 1

                # Next, process the reference sequence
                comp_start = min(0, patch_len)
                cur_ref = to_node[:-2]
                cur_ref_len = self.component_lens[cur_ref]
                cur_ref_strand = "+"
                if to_node.endswith("_e"):
                    cur_ref_strand = "-"
                agp.add_seq_line(obj_header, obj_pos+1, obj_pos+(cur_ref_len + comp_start), obj_pid, "W", cur_ref, 1+(-1*comp_start), cur_ref_len, cur_ref_strand)
                obj_pos += cur_ref_len + comp_start
                obj_pid += 1
                used_components.add(cur_ref)

                # Look for the next edge
                from_node = to_node[:-2] + "_b"
                if to_node.endswith("_b"):
                    from_node = to_node[:-2] + "_e"

                if from_node in self.graph.nodes:
                    next_nodes = set(self.graph[from_node])
                    assert len(next_nodes) == 1
                    to_node = next_nodes.pop()
                    used_edges.add((from_node, to_node))
                    used_edges.add((to_node, from_node))
                else:
                    next_edge_exists = False

        # Write unplaced reference sequences
        fai = pysam.FastaFile(ref_fn)
        all_ref_seqs = set(fai.references)
        fai.close()
        remaining_components = all_ref_seqs - used_components
        for c in sorted(remaining_components):
            agp.add_seq_line(
                c + "_RagTag" * add_suffix_to_unplaced,
                "1",
                str(self.component_lens[c]),
                "1",
                "W",
                c,
                "1",
                str(self.component_lens[c]),
                "+"
            )

        agp.write()
