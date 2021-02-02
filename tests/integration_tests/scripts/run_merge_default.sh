#!/usr/bin/env bash

# run ragtag merge with default settings

# position args:
## 1. output directory
## 2. assembly
## 3-n. AGP files


Usage() {
    echo "Usage: $0 output_dir asm.fa scf1.agp scf2.agp [...]"
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

ragtag.py merge --debug -u -l 0 -o $OUTDIR $ASM $AGP
