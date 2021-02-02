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

import re
import os
import bisect

import pysam


class RestrictionEnzymes:

    # Add restriction enzymes and corresponding sites here.
    enzymes_raw = [
        'HindIII',
        'Sau3AI',
        'MboI',
        'DpnII',
        'HinfI'
    ]

    # Use python regex syntax, not IUPAC ambiguity codes for wildcards
    # e.g. "[ATCG]" instead of "N" or "[CT]" instead of "Y"
    # Always use string literals i.e. r''
    sites = [
        r'AAGCTT',
        r'GATC',
        r'GATC',
        r'GATC',
        r'GA[ATCG]TC'
    ]

    # Keep everything upper case
    enzymes = [i.upper() for i in enzymes_raw]
    sites = [i.upper() for i in sites]

    assert len(enzymes) == len(sites)

    def __init__(self, res):
        """
        Blah
        :param res: List or set of restriction enzymes/sites
        """
        # The binary switch indicates the presence or absence of a restriction enzyme at a particular index
        self.switch = [0 for _ in range(len(self.enzymes))]

        # Remove empty elements and add enzymes/sites
        res = set(filter(None, res))
        if not res:
            raise ValueError("RestrictionEnzymes object must be instantiated with at least one restriction enzyme/site")
        self.add_from(res)

    def __str__(self):
        pairs = []
        for i, j in zip(self.get_enzymes_raw(), self.get_sites()):
            pairs.append("    {}: {}".format(i, j))

        return "\n".join(pairs)

    @staticmethod
    def get_info():
        """
        Print a help message showing all available enzymes/sites and other useful info.
        """
        msg = """
    Below is a list of all accepted restriction enzymes and
    their restriction sites:\n"""

        print(msg)
        pairs = []
        for i, j in zip(RestrictionEnzymes.enzymes_raw, RestrictionEnzymes.sites):
            pairs.append("        {}: {}".format(i, j))

        print("\n".join(pairs))

        msg = """
    For RagTag, use a comma separated list of enzymes or
    sites (or a mix). For example, for Arima Hi-C, use
    'Sau3AI,HinfI' or 'GATC,GA[ATCG]TC'."""
        print(msg)

        msg = """
    Note that for restriction sites, wildcards are
    represented with python regex syntax, not IUPAC
    ambiguity codes. e.g. '[ATCG]' instead of 'N'."""
        print(msg)

        msg = """
    Please contact the developers if you would like to add
    more enzymes/sites."""
        print(msg)

    def add(self, r):
        """
        Add a single restriction enzyme/site
        :param r:
        """
        if not isinstance(r, str):
            raise TypeError("r must be a string")

        if r:
            r = r.upper()
            if r in self.enzymes:
                self.switch[self.enzymes.index(r)] = 1
            elif r in self.sites:
                self.switch[self.sites.index(r)] = 1
            else:
                msg = "{} is not a valid restriction enzyme/site. Contact the developers if it should be.".format(r)
                raise ValueError(msg)

    def add_from(self, res):
        """
        Given a raw list/set of restriction enzymes/sites, set the master switch.
        :param res: list or set of restriction enzymes.sites
        """
        if not isinstance(res, set) and not isinstance(res, list):
            raise TypeError("res must be a list or set")

        for r in res:
            self.add(r)

    def get_all_enzymes(self):
        return self.enzymes

    def get_all_sites(self):
        return self.sites

    def get_enzymes(self):
        return [self.enzymes[i] for i in range(len(self.enzymes)) if self.switch[i]]

    def get_enzymes_raw(self):
        return [self.enzymes_raw[i] for i in range(len(self.enzymes_raw)) if self.switch[i]]

    def get_sites(self):
        return [self.sites[i] for i in range(len(self.sites)) if self.switch[i]]

    def make_regex(self):
        return "|".join(self.get_sites())


class RestrictionFragmentMap:
    """
    This class represents a restriction fragment map. It takes two arguments as input: a genome assembly file in fasta
    format and a RestrictionEnzymes objectt. The class performs an in silico restriction digest and provides
    additional helper functions.
    """

    def __init__(self, asm, res):
        if not isinstance(res, RestrictionEnzymes):
            raise TypeError("res must be a RestrictionEnzymes object")

        self.enzymes = res
        self.asm = os.path.abspath(asm)

        # Check that the assembly file exists
        if not os.path.isfile(self.asm):
            raise ValueError("Could not find file: %s" % self.asm)

        self.seq_lens = dict()
        self.site_pos = dict()
        self._find_sites()

    def _find_sites(self):
        """
        Store the index of all non-overlapping occurances of restriction sites
        """
        re_re = re.compile(self.enzymes.make_regex())
        fai = pysam.FastaFile(self.asm)

        # Iterate through each sequence in the assembly
        for ref in fai.references:
            self.seq_lens[ref] = fai.get_reference_length(ref)
            self.site_pos[ref] = []
            seq = fai.fetch(ref)

            # Find all non-overlapping sites within this sequence
            for match in re.finditer(re_re, seq):
                self.site_pos[ref].append(match.start(0))

    def count_sites_lte(self, seq, pos):
        if seq not in self.site_pos:
            raise ValueError("Invalid sequence: {}".format(seq))

        return bisect.bisect_right(self.site_pos[seq], pos)  # assumes the positions are sorted

    def count_sites_gt(self, seq, pos):
        if seq not in self.site_pos:
            raise ValueError("Invalide sequence: {}".format(seq))

        return len(self.site_pos[seq]) - bisect.bisect_right(self.site_pos[seq], pos)  # assumes sorted positions
