#!/usr/bin/env python

import os
import argparse
from collections import defaultdict

from intervaltree import IntervalTree


def sub_update(gff_file, bed_file):
    # Make a dictionary associating each original sequence with an interval tree of component sequences
    trans = defaultdict(IntervalTree)
    with open(bed_file, "r") as f:
        for line in f:
            oh, start, end, stype, nh, strand = line.rstrip().split("\t")
            if not strand == "+":
                raise ValueError("The placement BED file is not formatted correctly. All strand should be positive.")
            if not stype == "s":
                raise ValueError("The placement BED file is not formatted correctly. All lines should be sequence (s).")
            trans[oh][int(start):int(end)] = nh

    # Iterate through the gff intervals and update them according to trans
    with open(gff_file, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("#"):
                print(line)  # Print this comment line
            else:
                fields = line.split("\t")
                h, s, e = fields[0], int(fields[3]), int(fields[4])
                s -= 1  # Keep everything zero-indexed
                ovlps = trans[h][s:e]
                if len(ovlps) > 1:
                    raise ValueError(
                        "%s:%d-%d in the gff file overlaps two sub sequences in the placement file. Make sure to run `ragoo correct` with '--gff'" % (h, s, e)
                    )
                if len(ovlps) < 1:
                    raise ValueError("The placement BED file is not formatted correctly.")

                # Get the data from the overlapping interval and print the new line
                o = list(ovlps)[0]
                new_s = s - o.begin
                new_e = e - o.begin
                fields[0] = o.data
                fields[3] = str(new_s + 1)  # back to one-based indexing for gff format
                fields[4] = str(new_e)
                print("\t".join(fields))


def sup_update(gff_file, bed_file):
    raise NotImplementedError()


def main():
    parser = argparse.ArgumentParser(description="Update gff invervals given a 'placement' BED file")
    parser.add_argument("gff", metavar="<genes.gff>", type=str, help="gff file. must not be gzipped")
    parser.add_argument("bed", metavar="<placement.bed>", type=str, help="placement BED file")
    parser.add_arguments("-c", action="store_true", default=False, type=str, help="placement BED file")

    args = parser.parse_args()
    gff_file = os.path.abspath(args.gff)
    bed_file = os.path.abspath(args.bed)
    is_sub = args.c

    if is_sub:
        sub_update(gff_file, bed_file)
    else:
        sup_update(gff_file, bed_file)


if __name__ == "__main__":
    main()