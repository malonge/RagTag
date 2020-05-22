#!/bin/bash

# run the ragoo2 pipeline with default settings

# position args:
## 1. objects
## 2. components
## 3. agp

function Usage() {
    echo "Usage: $0 objects.fa components.fa seqs.agp"
}

if [ $# -lt 3 ] ; then
    Usage
    exit 1
fi

OBJS=$1
COMPS=$2
AGP=$3

agp_validate -comp -out objs.fasta $OBJS $COMPS $AGP

# Check that the generated fasta file is the same as the ragoo objects file
seqtk seq -A $OBJS > 1.fasta
python3 scripts/sort_fasta.py 1.fasta > 1.s.fasta

sed 's/>lcl|/>/g' objs.fasta > 2.fasta
seqtk seq -A 2.fasta > 2.r.fasta
python3 scripts/sort_fasta.py 2.r.fasta > 2.s.fasta

echo "Comparing fasta files with 'cmp':"
cmp 1.s.fasta 2.s.fasta

rm objs.fasta
rm 1.fasta 1.s.fasta 1.fasta.fai
rm 2.fasta 2.r.fasta.fai 2.r.fasta 2.s.fasta
