#!/usr/bin/env bash

# check that two AGP files define the same objects, not enforcing strand

# position args:
## 1. assembly
## 2. first AGP file
## 3. second AGP file

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

OUTDIR=$1
ASM=$2
AGP="${*:3}"

ragtag.py merge --debug -u -o $OUTDIR $ASM $AGP
