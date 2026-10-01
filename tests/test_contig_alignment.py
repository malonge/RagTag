import unittest

from ragtag_utilities.ContigAlignment import ContigAlignment


def make_alns(strand, query_starts, query_ends, ref_starts, ref_ends):
    """ Make a ContigAlignment for a 1 Mbp query with all alignments to one reference sequence. """
    n = len(query_starts)
    return ContigAlignment(
        "q",
        1000000,
        list(query_starts),
        list(query_ends),
        [strand] * n,
        ["r"] * n,
        [10000000] * n,
        list(ref_starts),
        list(ref_ends),
        [e - s for s, e in zip(ref_starts, ref_ends)],
        [e - s for s, e in zip(ref_starts, ref_ends)],
        [0] * n
    )


class TestCarefulMerge(unittest.TestCase):

    def test_forward_strand_merge(self):
        """ Three collinear forward-strand alignments with 10 kbp gaps merge into one alignment. """
        x = make_alns(
            "+",
            [0, 210000, 420000],
            [200000, 410000, 620000],
            [1000000, 1210000, 1420000],
            [1200000, 1410000, 1620000]
        )
        m = x.merge_alns(merge_dist=100000, careful_merge=True)
        self.assertEqual(m.num_alns, 1)
        self.assertEqual((m.query_starts[0], m.query_ends[0]), (0, 620000))
        self.assertEqual((m.ref_starts[0], m.ref_ends[0]), (1000000, 1620000))

    def test_reverse_strand_merge(self):
        """ The reverse-strand version of the forward-strand test also merges into one alignment. """
        x = make_alns(
            "-",
            [420000, 210000, 0],
            [620000, 410000, 200000],
            [1000000, 1210000, 1420000],
            [1200000, 1410000, 1620000]
        )
        m = x.merge_alns(merge_dist=100000, careful_merge=True)
        self.assertEqual(m.num_alns, 1)
        self.assertEqual((m.query_starts[0], m.query_ends[0]), (0, 620000))
        self.assertEqual((m.ref_starts[0], m.ref_ends[0]), (1000000, 1620000))

    def test_reverse_strand_query_gap(self):
        """ Reverse-strand alignments that are close on the reference but far apart on the query do not merge. """
        x = make_alns(
            "-",
            [700000, 0],
            [900000, 200000],
            [1000000, 1210000],
            [1200000, 1410000]
        )
        m = x.merge_alns(merge_dist=100000, careful_merge=True)
        self.assertEqual(m.num_alns, 2)

    def test_reverse_strand_orientation(self):
        """ The orientation follows the strand with most of the aligned bases. """
        # 600 kbp of reverse-strand alignments in 200 kbp pieces, then one 250 kbp forward-strand alignment.
        x = ContigAlignment(
            "q",
            1000000,
            [420000, 210000, 0, 700000],
            [620000, 410000, 200000, 950000],
            ["-", "-", "-", "+"],
            ["r"] * 4,
            [10000000] * 4,
            [1000000, 1210000, 1420000, 3000000],
            [1200000, 1410000, 1620000, 3250000],
            [200000, 200000, 200000, 250000],
            [200000, 200000, 200000, 250000],
            [0] * 4
        )
        m = x.merge_alns(merge_dist=100000, careful_merge=True)
        self.assertEqual(m.orientation, "-")
        self.assertGreater(m.orientation_confidence, 0.5)


if __name__ == "__main__":
    unittest.main()
