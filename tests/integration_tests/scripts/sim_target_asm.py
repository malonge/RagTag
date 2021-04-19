#!/usr/bin/env python

"""
This script creates test data for `ragtag patch`. First, we create a small query assembly. Then,
we induce changes (break sequences/add agaps etc.) in order to make a target assembly. `ragtag patch`
should be able to perfectly reconstruct the query assembly by patching the target assembly.
"""


import sys
import pysam

from ragtag_utilities.utilities import reverse_complement


# First argument is the reference used to make a query assembly

fasta_fn = sys.argv[1]
fai = pysam.FastaFile(fasta_fn)

# Write a small query assembly
with open("query.fa", "w") as f:
    f.write(">q1\n")
    f.write(fai.fetch("Chr1", 0, 5000000) + "\n")
    f.write(">q2\n")
    f.write(fai.fetch("Chr2", 0, 5000000) + "\n")

# Write a mutated target assembly
with open("target.fa", "w") as f:
    # Put a 500 bp gap in the first sequence
    f.write(">r1\n")
    f.write(fai.fetch("Chr1", 0, 1000000))
    f.write("N"*500)
    f.write(fai.fetch("Chr1", 1000500, 5000000) + "\n")

    # Break the second sequence into 5 sequences
    # The first should be reverse complemented
    # The second should be overlapping the first, + strand
    # The third should be reverse complemented
    # The fourth should be reverse complemnted, with some missing sequence between it and the third
    # The fifth should be + strand
    f.write(">r2a\n")
    f.write(fai.fetch("Chr2", 0, 1000000) + "\n")

    f.write(">r2b\n")
    f.write(fai.fetch("Chr2", 990000, 2000000) + "\n")

    f.write(">r2c\n")
    f.write(reverse_complement(fai.fetch("Chr2", 2000000, 3000000)) + "\n")

    f.write(">r2d\n")
    f.write(reverse_complement(fai.fetch("Chr2", 3001000, 4000000)) + "\n")

    f.write(">r2e\n")
    f.write(reverse_complement(fai.fetch("Chr2", 4000000, 5000000)) + "\n")

    # Write extra unscaffolded sequences
    f.write(">r3\n")
    f.write(fai.fetch("ChrM") + "\n")
    f.write(">r4\n")
    f.write(fai.fetch("ChrC"))