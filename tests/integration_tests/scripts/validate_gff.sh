#!/bin/bash

# run the ragoo2 pipeline with default settings

# position args:
## 1. pre fasta
## 2. pre gff
## 3. post fasta
## 4. post gff

# Assumes fasta suffix is ".fasta"
# Assumes gff suffix is ".gff"

function Usage() {
    echo "Usage: $0 pre.fa pre.gff post.fasta post.gff"
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


python3 scripts/gff2fasta.py $PRE.fasta $PREG.gff > $PREG.fasta
python3 scripts/sort_fasta.py $PREG.fasta > $PREG.s.fasta

python3 scripts/gff2fasta.py $POST.fasta $POSTG.gff > $POSTG.fasta
python3 scripts/sort_fasta.py $POSTG.fasta > $POSTG.s.fasta


echo "Comparing gene fasta files with 'cmp':"
cmp $PREG.s.fasta $POSTG.s.fasta

rm $PREG.fasta $PREG.fasta.fai $PREG.s.fasta
rm $POSTG.fasta $POSTG.fasta.fai $POSTG.s.fasta
