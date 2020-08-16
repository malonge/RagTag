#!/usr/bin/env bash

# Validate that the GFF updating is accurate

# position args:
## 1. pre fasta
## 2. pre gff
## 3. post fasta
## 4. post gff

# Assumes fasta suffix is ".fasta"
# Assumes gff suffix is ".gff"

Usage() {
    echo "Usage: $0 pre.fa pre.gff post.fasta post.gff"
}

mecho() {
    NAME=`basename $0`
    echo "$NAME:" $1
}

if [ $# -lt 4 ] ; then
    Usage
    exit 1
fi

PREFILE=$1
PRE=${PREFILE%.*}

PREGFILE=$2
PREG=${PREGFILE%.*}

POSTFILE=$3
POST=${POSTFILE%.*}

POSTGFILE=$4
POSTG=${POSTGFILE%.*}

/Library/Frameworks/Python.framework/Versions/3.6/bin/python3 scripts/gff2fasta.py $PRE.fasta $PREG.gff > $PREG.fasta
/Library/Frameworks/Python.framework/Versions/3.6/bin/python3 scripts/sort_fasta.py $PREG.fasta > $PREG.s.fasta

/Library/Frameworks/Python.framework/Versions/3.6/bin/python3 scripts/gff2fasta.py $POST.fasta $POSTG.gff > $POSTG.fasta
/Library/Frameworks/Python.framework/Versions/3.6/bin/python3 scripts/sort_fasta.py $POSTG.fasta > $POSTG.s.fasta


mecho "Comparing gene fasta files with 'cmp':"
cmp $PREG.s.fasta $POSTG.s.fasta

rm $PREG.fasta $PREG.fasta.fai $PREG.s.fasta
rm $POSTG.fasta $POSTG.fasta.fai $POSTG.s.fasta
