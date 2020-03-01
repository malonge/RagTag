from ragoo_utilities.utilities import summarize_planesweep


class ContigAlignment:
    """
    Description
    """

    def __init__(self, in_query_header, in_query_len, in_reference_headers, in_ref_lens, in_ref_starts, in_ref_ends, in_query_starts, in_query_ends, in_strands, in_mapqs):
        # Query info
        self.query_header = in_query_header
        self.query_len = in_query_len

        # Attributes describing alignments
        self._ref_headers = in_reference_headers
        self._ref_lens = in_ref_lens
        self._ref_starts = in_ref_starts
        self._ref_ends = in_ref_ends
        self._query_starts = in_query_starts
        self._query_ends = in_query_ends
        self._strands = in_strands
        self._mapqs = in_mapqs

        # Check that the dimensions are valid
        all_lens = self._get_attr_lens()
        if not len(set(all_lens)) == 1:
            raise ValueError("The alignments are incomplete.")

        # Attributes derived from alignments
        self.best_ref_header = None
        self.grouping_confidence = None
        self.orientation = None
        self.orientation_confidence = None

    def __str__(self):
        alns = []
        for i in range(len(self._ref_headers)):
            alns.append("\t".join([
                self.query_header,
                str(self.query_len),
                str(self._query_starts[i]),
                str(self._query_ends[i]),
                self._strands[i],
                self._ref_headers[i],
                str(self._ref_lens[i]),
                str(self._ref_starts[i]),
                str(self._ref_ends[i]),
                str(self._mapqs[i])
            ]))
        return "\n".join(alns)

    def _get_attr_lens(self):
        all_lens = [
            len(self._ref_headers),
            len(self._ref_lens),
            len(self._ref_starts),
            len(self._ref_ends),
            len(self._query_starts),
            len(self._query_ends),
            len(self._strands),
            len(self._mapqs)
        ]
        return all_lens

    def _rearrange_alns(self, hits):
        """ Order the alignments according to 'hits', an ordered list of indices. """
        return ContigAlignment(
            self.query_header,
            self.query_len,
            [self._ref_headers[i] for i in hits],
            [self._ref_lens[i] for i in hits],
            [self._ref_starts[i] for i in hits],
            [self._ref_ends[i] for i in hits],
            [self._query_starts[i] for i in hits],
            [self._query_ends[i] for i in hits],
            [self._strands[i] for i in hits],
            [self._mapqs[i] for i in hits]
        )

    def add_alignment(self, in_reference_header, in_ref_len, in_ref_start, in_ref_end, in_query_start, in_query_end, in_strand, in_mapq):
        """ Add an alignment for this query. """
        return ContigAlignment(
            self.query_header,
            self.query_len,
            self._ref_headers + [in_reference_header],
            self._ref_lens + [in_ref_len],
            self._ref_starts + [in_ref_start],
            self._ref_ends + [in_ref_end],
            self._query_starts + [in_query_start],
            self._query_ends + [in_query_end],
            self._strands + [in_strand],
            self._mapqs + [in_mapq]
        )

    def filter_lengths(self, l):
        if not isinstance(l, int):
            raise ValueError("l must be an integer. ")
        hits = [i for i in range(len(self._ref_headers)) if abs(self._query_ends[i] - self._query_starts[i]) >= l]
        return self._rearrange_alns(hits)

    def filter_mapq(self, q):
        if not isinstance(q, int):
            raise ValueError("q must be an integer. ")
        hits = [i for i in range(len(self._ref_headers)) if self._mapqs[i] >= q]
        return self._rearrange_alns(hits)

    def unique_anchor_filter(self, unique_length):
        """
        The contents of this method are either influenced by or directly copied from "Assemblytics_uniq_anchor.py"
        written by Maria Nattestad. The original script can be found here:

        https://github.com/MariaNattestad/Assemblytics

        And the publication associated with Maria's work is here:

        Nattestad, Maria, and Michael C. Schatz. "Assemblytics: a
        web analytics tool for the detection of variants from an
        assembly." Bioinformatics 32.19 (2016): 3021-3023.
        """
        lines_by_query = []
        for i, j in zip(self._query_starts, self._query_ends):
            lines_by_query.append((i, j))

        hits = summarize_planesweep(lines_by_query, unique_length)
        return self._rearrange_alns(hits)
