#!/usr/bin/env bash

# check that two AGP files define the same objects, not enforcing strand

# position args:
## 1. assembly
## 2. first AGP file
## 3. second AGP file

Usage() {
    echo "Usage: $0 asm.fasta scf1.agp scf2.agp"
}

mecho() {
    NAME=`basename $0`
    echo "$NAME:" $1
}

if [ $# -lt 2 ] ; then
    Usage
    exit 1
fi

ASM=$1
AGP1=$2
AGP2=$3

# Prep the objects for the first AGP file
ragtag.py agp2fasta $AGP1 $ASM > 1.fasta
python3 scripts/choose_strand.py 1.fasta > 1.s.fasta
grep -v ">"  1.s.fasta | sort > 1.seq

# Prep the objects for the second AGP file
ragtag.py agp2fasta $AGP2 $ASM > 2.fasta
python3 scripts/choose_strand.py 2.fasta > 2.s.fasta
grep -v ">"  2.s.fasta | sort > 2.seq

echo ""
mecho "Comparing sequence files with 'cmp':"

cmp 1.seq 2.seq

mecho "done"

rm 1.fasta.fai 2.fasta.fai
rm 1.fasta 1.s.fasta 2.fasta 2.s.fasta
rm 1.seq 2.seq