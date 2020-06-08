#!/usr/bin/env bash

# run the ragtag pipeline with nucmer

# position args:
## 1. reference
## 2. query
## 3. gff
## 4. output dir

# Assumes query suffix is ".fasta"
# Assumes gff suffix is ".gff"

Usage() {
    echo "Usage: $0 ref.fa query.fa genes.gff output_dir"
}

mecho() {
    NAME=`basename $0`
    echo "$NAME:" $1
}

if [ $# -lt 3 ] ; then
    Usage
    exit 1
fi

REF=$1
QUERY=$2
QUERY_PREF=`basename $QUERY .fasta`
GENES=$3
GENES_PREF=`basename $GENES .gff`
OUTDIR=$4

# TODO add '-w' flag to all commands

ragtag.py correct --aligner nucmer --nucmer-params "-l 500 -c 1000" --debug -f 1000 -u --gff $GENES -o $OUTDIR $REF $QUERY
ragtag.py scaffold --aligner nucmer --nucmer-params "-l 500 -c 1000" --debug -f 1000 -u -C -o $OUTDIR $REF $OUTDIR/$QUERY_PREF.corrected.fasta
ragtag.py updategff -c $GENES $OUTDIR/ragtag.correction.agp > $OUTDIR/$GENES_PREF.corr.gff
ragtag.py updategff $OUTDIR/$GENES_PREF.corr.gff $OUTDIR/ragtag.scaffolds.agp > $OUTDIR/$GENES_PREF.scaf.gff
