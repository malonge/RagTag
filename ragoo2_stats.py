#!/usr/bin/env python

import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description="Calculate scaffolding statistics")
    parser.add_argument("orderings", nargs='?', default="", metavar="<orderings.bed>", type=str, help="RaGOO ordering file")
    parser.add_argument("output_file", nargs='?', default="", metavar="<stats.txt>", type=str, help="output file name")

    args = parser.parse_args()

    if not args.orderings or not args.output_file:
        parser.print_help()
        sys.exit()

    orderings_file = args.orderings
    output_file = args.output_file

    placed_bp = 0
    placed_seq = 0
    unplaced_bp = 0
    unplaced_seq = 0

    with open(orderings_file, "r") as f:
        for line in f:
            l = line.rstrip().split("\t")
            if l[3] == "s":
                seq_len = int(l[2]) - int(l[1])
                if l[6] != "NA":
                    placed_bp += seq_len
                    placed_seq += 1
                else:
                    unplaced_bp += seq_len
                    unplaced_seq += 1

    with open(output_file, "w") as f:
        f.write("placed_sequences\tplaced_bp\tunplaced_sequences\tunplaced_bp\n")
        f.write("\t".join([str(placed_seq), str(placed_bp), str(unplaced_seq), str(unplaced_bp)]) + "\n")


if __name__ == "__main__":
    main()
