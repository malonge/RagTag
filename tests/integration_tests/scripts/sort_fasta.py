import sys

import pysam

# Only one argument: FASTA file

x = pysam.FastaFile(sys.argv[1])
    
for i in sorted(x.references):
    print(">" + i)
    print(x.fetch(i))
