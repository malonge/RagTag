#!/usr/bin/env python


class AlignmentLine:
    """ Object to represent a single Minimap2 or Nucmer alignment. """

    def __init__(self, in_query_header, in_query_len, in_query_start, in_query_end, in_strand, in_ref_header, in_ref_len, in_ref_start, in_ref_end, in_num_match, in_aln_len, in_mapq):
        """
        Start positions should be before end positions for both query and target
        """
        self.query_header = in_query_header
        self.query_len = in_query_len
        self.query_start = in_query_start
        self.query_end = in_query_end
        self.strand = in_strand
        self.ref_header = in_ref_header
        self.ref_len = in_ref_len
        self.ref_start = in_ref_start
        self.ref_end = in_ref_end
        self.num_match =in_num_match
        self.aln_len = in_aln_len
        self.mapq = in_mapq

        assert self.query_start <= self.query_end
        assert self.ref_start <= self.ref_end


class AlignmentReader:

    def __init__(self, aln_file):
        self.aln_file = aln_file

    def parse_alignments(self):
        with open(self.aln_file) as f:
            for line in f:
                line = line.rstrip().split("\t")
                yield AlignmentLine(
                    line[0],
                    int(line[1]),
                    int(line[2]),
                    int(line[3]),
                    line[4],
                    line[5],
                    int(line[6]),
                    int(line[7]),
                    int(line[8]),
                    int(line[9]),
                    int(line[10]),
                    int(line[11])
                )