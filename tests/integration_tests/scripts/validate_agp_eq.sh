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
ragtag.py agp2fasta $ASM $AGP1 > 1.fasta
ptyhon3 scripts/choose_strand.py 1.fasta > 1.s.fasta
python3 scripts/sort_fasta.py 1.s.fasta > 1.s.s.fasta

# Prep the objects for the second AGP file
ragtag.py agp2fasta $ASM $AGP2 > 2.fasta
ptyhon3 scripts/choose_strand.py 2.fasta > 2.s.fasta
python3 scripts/sort_fasta.py 2.s.fasta > 2.s.s.fasta

echo ""
mecho "Comparing fasta files with 'cmp':"

cmp 1.s.s.fasta 2.s.s.fasta

rm 1.fasta 1.s.fasta 1.s.s.fasta 2.fasta 2.s.fasta 2.s.s.fasta