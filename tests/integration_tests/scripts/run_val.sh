#!/usr/bin/env bash

# run the ragtag correction with read validation (error corrected long reads)

# position args:
## 1. reference
## 2. query
## 3. reads.fofn
## 4. output dir

# Assumes query suffix is ".fasta"

Usage() {
    echo "Usage: $0 ref.fa query.fa reads.fofn output_dir"
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
READS=$3
OUTDIR=$4

ragtag.py correct -F $READS -T corr --debug -t 2 -u -o $OUTDIR $REF $QUERY