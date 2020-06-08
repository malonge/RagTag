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

import operator
from collections import defaultdict

from ragtag_utilities.utilities import summarize_planesweep, p2q, q2p


class ContigAlignment:
    """
    A ContigAlignment object stores and organizes all alignments pertaining to a single query sequence.

    The alignment organization is guided by the PAF file format: https://github.com/lh3/miniasm/blob/master/PAF.md. The
    only difference is that ContigAlignment uses the word "reference" instead of "target" to refer to the sequences
    to which the query is aligned.

    Each instance of a ContigAlignment object will refer to a single query sequence (self.query_header) and associated
    query sequence length (self.query_len). The remaining PAF fields are stored in 10 lists for the 10 remaining fields:

    1. query start,       2. query end,        3. strand
    4. reference header,  5. reference length, 6. reference start,  7. reference end
    8. # residue matches, 9. alignment length, 10. mapq

    The offset in each of these lists corresponds to a single alignment. For example, the first element in each of the
    above-mentioned lists provides the necessary information to describe the first alignment stored in the object.

    The set of alignments organized in a ContigAlignment object is used to calculate metrics for the query sequence. For
    example, self.best_ref_header gives the reference sequence which the query sequence most covers. For this reason,
    the ContigAlignment object is pseudo-immutable. The order of alignments may change (e.g. _sort_by_ref() and
    _sort_by_query()) but any addition/subtraction of alignments is not allowed. Any such addition/removal functionality
    (e.g. add_alignment()) must return a new instance of the class.
    """

    def __init__(self, in_query_header, in_query_len, in_query_starts, in_query_ends, in_strands, in_reference_headers, in_ref_lens, in_ref_starts, in_ref_ends, in_residue_matches, in_aln_lens, in_mapqs):
        self.query_header = in_query_header
        self.query_len = in_query_len
        self._query_starts = in_query_starts
        self._query_ends = in_query_ends
        self._strands = in_strands
        self._ref_headers = in_reference_headers
        self._ref_lens = in_ref_lens
        self._ref_starts = in_ref_starts
        self._ref_ends = in_ref_ends
        self._residue_matches = in_residue_matches
        self._aln_lens = in_aln_lens
        self._mapqs = in_mapqs

        # Check that the dimensions are valid
        all_lens = self._get_attr_lens()
        if not len(set(all_lens)) == 1:
            raise ValueError("The alignments are incomplete.")
        if not all_lens[0]:
            raise ValueError("ContigAlignment must contain at least one alignment.")

        # Get alignment lengths with respect to the reference
        self._ref_aln_lens = None
        self._set_ref_aln_lens()

        # Attributes derived from alignments
        self.best_ref_header = None
        self.grouping_confidence = None
        self._get_best_ref_header()

        self.primary_alignment = None
        self._get_primary_alignment()

        self.orientation = self._strands[self.primary_alignment]
        self.orientation_confidence = None
        self._get_orientation_confidence()
        self.location_confidence = None
        self._get_location_confidence()

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
                str(self._residue_matches[i]),
                str(self._aln_lens[i]),
                str(self._mapqs[i])
            ]))

        return "\n".join(alns) + "\n"

    @staticmethod
    def _average_mapqs(q1, q2):
        """ Calculate the average of two mapping quality values. """
        p1, p2 = q2p(q1), q2p(q2)
        return p2q((p1 + p2) / 2)

    def _get_attr_lens(self):
        """ Return the lengths of all alignment field lists. """
        all_lens = [
            len(self._query_starts),
            len(self._query_ends),
            len(self._strands),
            len(self._ref_headers),
            len(self._ref_lens),
            len(self._ref_starts),
            len(self._ref_ends),
            len(self._residue_matches),
            len(self._aln_lens),
            len(self._mapqs)
        ]
        return all_lens

    def _set_ref_aln_lens(self):
        """ For each alignment, set the alignment length w.r.t the reference coordinates."""
        self._ref_aln_lens = [i-j for i, j in zip(self._ref_ends, self._ref_starts)]

    def _get_best_ref_header(self):
        """ From the alignments, determine the reference sequence that is the most covered by this query sequence. """
        all_ref_headers = set(self._ref_headers)
        if len(all_ref_headers) == 1:
            self.best_ref_header = self._ref_headers[0]
            self.grouping_confidence = 1.0
            return

        # Get all the alignment intervals for each reference sequence
        all_intervals = defaultdict(list)
        for i in range(len(self._ref_headers)):
            this_range = (self._ref_starts[i], self._ref_ends[i])
            this_seq = self._ref_headers[i]
            all_intervals[this_seq].append(this_range)

        # For each reference header, sort the intervals and get the union interval length.
        ranges = defaultdict(int)
        for i in all_intervals.keys():
            sorted_intervals = sorted(all_intervals[i], key=lambda tup: tup[0])
            max_end = -1
            for j in sorted_intervals:
                start_new_terr = max(j[0], max_end)
                ranges[i] += max(0, j[1] - start_new_terr)
                max_end = max(max_end, j[1])

        # Convert to a list and sort the ranges.items() in order to have ties broken in a deterministic way.
        max_seq = max(sorted(list(ranges.items())), key=operator.itemgetter(1))[0]
        self.best_ref_header = max_seq

        # Now get the confidence of this chromosome assignment
        # Equal to the max range over all ranges
        self.grouping_confidence = ranges[max_seq] / sum(ranges.values())

    def _get_ref_alns(self, r):
        """ Provide the offsets for alignments to a specified reference sequence. """
        return [i for i in range(len(self._ref_headers)) if self._ref_headers[i] == r]

    def _get_best_ref_alns(self):
        """ Provide the offsets for alignments to the 'best' reference sequence. """
        return self._get_ref_alns(self.best_ref_header)

    def _get_primary_alignment(self):
        """
        Get the offset corresponding to the 'primary' (longest) alignment.
        Only consider alignments to the 'best' reference sequence.
        """
        max_index = -1
        max_len = -1
        for i in self._get_best_ref_alns():
            this_len = self._ref_aln_lens[i]
            if this_len > max_len:
                max_len = this_len
                max_index = i
        self.primary_alignment = max_index

    def _get_orientation_confidence(self):
        """ Calculate the orientation confidence score. """
        num = 0
        denom = 0
        for i in self._get_best_ref_alns():
            aln_len = self._ref_aln_lens[i]
            if self._strands[i] == self.orientation:
                num += aln_len
            denom += aln_len
        self.orientation_confidence = num/denom

    def _get_location_confidence(self):
        """ Calculate the location confidence score. """
        best_ref_alns = self._get_best_ref_alns()

        # Get all the alignment reference intervals for alignments to the best reference sequence
        aln_intervals = []
        all_positions = []
        for i in best_ref_alns:
            aln_intervals.append((self._ref_starts[i], self._ref_ends[i]))
            all_positions.append(self._ref_starts[i])
            all_positions.append(self._ref_ends[i])

        # The denominator is the max - min alignment positions
        denom = max(all_positions) - min(all_positions)

        # The numerator is the coverage
        num = 0
        sorted_intervals = sorted(aln_intervals, key=lambda tup: tup[0])
        max_end = -1
        for j in sorted_intervals:
            start_new_terr = max(j[0], max_end)
            num += max(0, j[1] - start_new_terr)
            max_end = max(max_end, j[1])

        self.location_confidence = num/denom

    def _update_alns(self, hits):
        """ Order alignments according to 'hits', an ordered list of offsets. Return a new instance of the class. """
        if hits:
            return ContigAlignment(
                self.query_header,
                self.query_len,
                [self._query_starts[i] for i in hits],
                [self._query_ends[i] for i in hits],
                [self._strands[i] for i in hits],
                [self._ref_headers[i] for i in hits],
                [self._ref_lens[i] for i in hits],
                [self._ref_starts[i] for i in hits],
                [self._ref_ends[i] for i in hits],
                [self._residue_matches[i] for i in hits],
                [self._aln_lens[i] for i in hits],
                [self._mapqs[i] for i in hits]
            )
        else:
            return None

    def _rearrange_alns(self, hits):
        """ Order the alignments according to 'hits', an ordered list of indices. """
        if len(hits) != len(self._ref_headers):
            raise ValueError("Can only order alignments. To update, use '_update_alns()'")

        self._query_starts = [self._query_starts[i] for i in hits]
        self._query_ends = [self._query_ends[i] for i in hits]
        self._strands = [self._strands[i] for i in hits]
        self._ref_headers = [self._ref_headers[i] for i in hits]
        self._ref_lens = [self._ref_lens[i] for i in hits]
        self._ref_starts = [self._ref_starts[i] for i in hits]
        self._ref_ends = [self._ref_ends[i] for i in hits]
        self._residue_matches = [self._residue_matches[i] for i in hits]
        self._aln_lens = [self._aln_lens[i] for i in hits]
        self._mapqs = [self._mapqs[i] for i in hits]

        self._set_ref_aln_lens()

    def _sort_by_ref(self):
        """ Sort the alignments by reference header/position. """
        ref_pos = []
        for i in range(len(self._ref_headers)):
            ref_pos.append((self._ref_headers[i], self._ref_starts[i], self._ref_ends[i], i))
        hits = [i[3] for i in sorted(ref_pos)]

        self._rearrange_alns(hits)

    def _sort_by_query(self):
        """ Sort the alignments by query position. """
        q_pos = []
        for i in range(len(self._ref_headers)):
            q_pos.append((self._query_starts[i], self._query_ends[i], i))
        hits = [i[2] for i in sorted(q_pos)]

        self._rearrange_alns(hits)

    def add_alignment(self, in_query_start, in_query_end, in_strand, in_reference_header, in_ref_len, in_ref_start, in_ref_end, in_residue_matches, in_aln_len, in_mapq):
        """ Add an alignment for this query sequence. Return a new instance of the class. """
        return ContigAlignment(
            self.query_header,
            self.query_len,
            self._query_starts + [in_query_start],
            self._query_ends + [in_query_end],
            self._strands + [in_strand],
            self._ref_headers + [in_reference_header],
            self._ref_lens + [in_ref_len],
            self._ref_starts + [in_ref_start],
            self._ref_ends + [in_ref_end],
            self._residue_matches + [in_residue_matches],
            self._aln_lens + [in_aln_len],
            self._mapqs + [in_mapq]
        )

    def filter_lengths(self, l):
        """ Remove alignments shorter than l. """
        hits = [i for i in range(len(self._ref_headers)) if self._ref_aln_lens[i] >= l]
        return self._update_alns(hits)

    def filter_mapq(self, q):
        """ Remove alignments with mapq < q. """
        hits = [i for i in range(len(self._ref_headers)) if self._mapqs[i] >= q]
        return self._update_alns(hits)

    def unique_anchor_filter(self, l, keep_small=True):
        """
        Unique anchor filter the alignments. l is the minimum unique alignment length. small_uniques are retained.

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

        hits = summarize_planesweep(lines_by_query, l, keep_small_uniques=keep_small)
        return self._update_alns(hits)

    def get_best_ref_pos(self):
        """ Return the ref start and ref end for the primary alignment. """
        return self._ref_starts[self.primary_alignment], self._ref_ends[self.primary_alignment]

    def get_best_q_dist(self):
        """ Find the distance from the primary alignment to the beginning and end of the query sequence. """
        ld = self._query_starts[self.primary_alignment]
        ed = self.query_len - self._query_ends[self.primary_alignment]

        if self.orientation == "-":
            ld, ed = ed, ld

        return ld, ed

    def filter_query_contained(self):
        """
        Remove alignments that are contained (w.r.t the query) within other alignments.
        This does not consider alignments contained by chains of alignments.
        Consider merging first (merge_alns()) to account for that.
        """
        cidx = set()
        for i in range(len(self._ref_headers)):
            for j in range(len(self._ref_headers)):
                # Check if j is contained by i
                if i == j:
                    continue
                if self._query_starts[i] <= self._query_starts[j] and self._query_ends[i] >= self._query_ends[j]:
                    if i not in cidx:  # Keeps at least one of many identical coordinates
                        cidx.add(j)

        hits = [i for i in range(len(self._ref_headers)) if i not in cidx]
        return self._update_alns(hits)

    def merge_alns(self, merge_dist=100000):
        """
        Merge adjacent alignments that have the same reference sequence, the same orientation, and are less than
        merge_dist away from each other.
        """
        # Sort the alignments
        self._sort_by_ref()

        # Make a copy of the alignment info
        query_starts = self._query_starts
        query_ends = self._query_ends
        strands = self._strands
        ref_headers = self._ref_headers
        ref_lens = self._ref_lens
        ref_starts = self._ref_starts
        ref_ends = self._ref_ends
        residue_matches = self._residue_matches
        aln_lens = self._aln_lens
        mapqs = self._mapqs

        # Keep track of which alignments we are comparing
        i = 0
        j = 1
        while j < len(ref_headers):
            if all([
                        ref_headers[i] == ref_headers[j],
                        strands[i] == strands[j],
                        ref_starts[j] - ref_ends[i] <= merge_dist
            ]):
                # Merge the alignments in place of the first alignment
                query_starts[i] = min(query_starts[i], query_starts[j])
                query_ends[i] = max(query_ends[i], query_ends[j])
                ref_starts[i] = min(ref_starts[i], ref_starts[j])
                ref_ends[i] = max(ref_ends[i], ref_ends[j])

                perc_r_match_i = residue_matches[i] / aln_lens[i]
                perc_r_match_j = residue_matches[j] / aln_lens[j]
                avg_perc_r_match = (perc_r_match_i + perc_r_match_j) / 2
                aln_lens[i] = ref_ends[i] - ref_starts[i]
                residue_matches[i] = round(avg_perc_r_match * aln_lens[i])

                mapqs[i] = self._average_mapqs(mapqs[i], mapqs[j])

                # Remove the redundant alignment
                query_starts.pop(j)
                query_ends.pop(j)
                strands.pop(j)
                ref_headers.pop(j)
                ref_lens.pop(j)
                ref_starts.pop(j)
                ref_ends.pop(j)
                residue_matches.pop(j)
                aln_lens.pop(j)
                mapqs.pop(j)
            else:
                i += 1
                j += 1

        # Make a new object with the merged data
        x = ContigAlignment(
            self.query_header,
            self.query_len,
            query_starts,
            query_ends,
            strands,
            ref_headers,
            ref_lens,
            ref_starts,
            ref_ends,
            residue_matches,
            aln_lens,
            mapqs
        )

        # remove contained alignments.
        return x.filter_query_contained()

    def get_break_candidates(self, min_dist=5000):
        """
        Return coordinates of the query sequence between consecutive alignments.
        Consider merging alignments first (merge_alns())
        :return: Two lists of coordinates. One where consecutive alignments aligned to the
        same reference (intra) and one where the consecutive alignments aligned to different
        references (inter).
        """
        self._sort_by_query()
        intra_candidates = []
        inter_candidates = []

        # If there are more than two alignments, iterate through all but the first.
        for i in range(1, len(self._ref_headers)):
            if min_dist < self._query_starts[i] < (self.query_len - min_dist):
                if self._ref_headers[i] == self._ref_headers[i-1]:
                    intra_candidates.append(self._query_starts[i])
                else:
                    inter_candidates.append(self._query_starts[i])

        return intra_candidates, inter_candidates

