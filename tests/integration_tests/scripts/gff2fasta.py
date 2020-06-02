import argparse

import pysam

from ragtag_utilities.utilities import reverse_complement

"""
Like bedtools getfasta, but use the gff ID attribute as the FASTA header and
always force strandedness. 
"""


def main():
    parser = argparse.ArgumentParser(description="Get fasta sequences from a GFF file")
    parser.add_argument("fasta", metavar="<sequences.fasta>", type=str, help="AGP v2.1 file")
    parser.add_argument("gff", metavar="<genes.gff>", type=str, help="FASTA file with component sequences to be scaffolded. must not be gzipped")
    
    args = parser.parse_args()
    fasta_file = args.fasta
    gff_file = args.gff
    
    x = pysam.FastaFile(fasta_file)
    
    # Read the gff file
    with open(gff_file, "r") as f:
        for line in f:
            if not line.startswith("#"):
                seqname, source, feature, start, end, score, strand, fname, attributes = line.rstrip().split("\t")
                start, end = int(start)-1, int(end)
                
                # Get the ID attribute
                gff_id = None
                tags = attributes.split(";")
                for j in tags:
                    if j.startswith("ID="):
                       gff_id = j[3:]
                       
                if gff_id is None:
                    raise ValueError("Need an ID attribute for each gff line.")
                    
                print(">" + gff_id)
                if strand == "+":
                    print(x.fetch(seqname, start, end))
                elif strand == "-":
                    print(reverse_complement(x.fetch(seqname, start, end)))
                else:
                    raise ValueError("Incorrect strand value")

if __name__ == "__main__":
    main()
