#!/usr/bin/env bash

# position args:
## 1. produced AGP file
## 2. static AGP file

Usage() {
    echo "Usage: $0 new.agp old.agp"
}

mecho() {
    NAME=`basename $0`
    echo "$NAME:" $1
}

if [ $# -lt 2 ] ; then
    Usage
    exit 1
fi

NEW=$1
OLD=$2

mecho "comparing AGP files with cmp:"

cmp $1 $2