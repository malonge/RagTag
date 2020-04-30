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

    # Compile the AGP lines for the placed sequences
    uid = 1
    out_lines = []
    with open(orderings_file, "r") as f:
        for line in f:
            l = line.rstrip().split("\t")
            ref_header = l[0]
            ref_start, ref_end = int(l[1]), int(l[2])
            query_len = ref_end - ref_start
            ref_start = ref_start + 1  # AGP is 1-indexed
            seq_type = l[3]
            query_header = l[4]
            strand = l[5]

            if seq_type == "s":
                this_line = [ref_header, str(ref_start), str(ref_end), str(uid), "W", query_header, "1", str(query_len), strand]
            else:
                this_line = [ref_header, str(ref_start), str(ref_end), str(uid), "N", str(query_len), "scaffold", "yes", "align_genus"]

            out_lines.append("\t".join(this_line))
            uid += 1

    with open(output_file, "w") as f:
        f.write("## AGP-version 2.1\n")
        f.write("## AGP constructed by RaGOO2\n")
        f.write("\n".join(out_lines) + "\n")


if __name__ == "__main__":
    main()
