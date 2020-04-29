#!/usr/bin/env python

import argparse


def main():
    parser = argparse.ArgumentParser(description="Calculate scaffolding statistics")
    parser.add_argument("orderings", metavar="<orderings.bed>", type=str, help="RaGOO ordering file")
    parser.add_argument("unplaced_file", metavar="<unplaced.txt>", type=str, help="unplaced sequence file")
    parser.add_argument("output_file", metavar="<stats.txt>", type=str, help="output file name")
    parser.add_argument("gap_len", metavar="<100>", type=int, help="gap length used in RaGOO")
    parser.add_argument("-C", action="store_true", default=False, help="concatenate unplaced contigs and make 'chr0'")

    args = parser.parse_args()
    orderings_file = args.orderings
    unplaced_file = args.unplaced_file
    output_file = args.output_file
    gap_len = args.gap_len
    make_chr0 = args.C

    out_lines = []

    # Compile the AGP lines for the placed sequences
    uid = 1
    with open(orderings_file, "r") as f:
        for line in f:
            l = line.rstrip().split("\t")
            ref_header = l[0]

            ref_start = int(l[1])
            ref_end = int(l[2])
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

    # Compile the AGP lines for unplaced sequences
    with open(unplaced_file, "r") as f:
        if make_chr0:
            pos = 0
            for line in f:
                l = line.rstrip().split("\t")
                unplaced_header, unplaced_len = l[0], int(l[1])

                # The sequence line
                this_line = ["Chr0_RaGOO", str(pos + 1), str(pos + unplaced_len), str(uid), "W", unplaced_header, "1", str(unplaced_len), "+"]
                out_lines.append("\t".join(this_line))
                pos += unplaced_len
                uid += 1

                # Next the gap line
                this_line = ["Chr0_RaGOO", str(pos + 1), str(pos + gap_len), str(uid), "N", str(gap_len), "scaffold", "no", "na"]
                out_lines.append("\t".join(this_line))
                pos += gap_len
                uid += 1

            # Remove the last gap line
            out_lines = out_lines[:-1]
        else:
            for line in f:
                l = line.rstrip().split("\t")
                unplaced_header, unplaced_len = l[0], int(l[1])
                this_line = [unplaced_header, "1", str(unplaced_len), str(uid), "W", unplaced_header, "1", str(unplaced_len), "+"]
                uid += 1
                out_lines.append("\t".join(this_line))

    with open(output_file, "w") as f:
        f.write("## AGP-version 2.1\n")
        f.write("## AGP constructed by RaGOO2\n")
        f.write("\n".join(out_lines) + "\n")


if __name__ == "__main__":
    main()
