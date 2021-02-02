#!/usr/bin/env python

import sys

import pysam

from ragtag_utilities.utilities import reverse_complement

# Only one argument: FASTA file

x = pysam.FastaFile(sys.argv[1])

for i in x.references:
    print(">" + i)
    s1 = x.fetch(i)
    s2 = reverse_complement(s1)
    if s1 < s2:
        print(s1)
    else:
        print(s2)

x.close()