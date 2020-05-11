#!/usr/bin/env python

import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description="Calculate scaffolding statistics")
    parser.add_argument("agp", nargs='?', default="", metavar="<ragoo2.scaffolds.agp>", type=str, help="RaGOO2 scaffolding AGP file")
    parser.add_argument("confidence", nargs='?', default="", metavar="<ragoo2.confidence.txt>", type=str, help="RaGOO2 scaffolding confidence scores file")

    args = parser.parse_args()

    if not args.agp or not args.confidence:
        parser.print_help()
        sys.exit()

    agp_file = args.agp
    confidence_file = args.confidence

    placed_bp = 0
    placed_seq = 0
    unplaced_bp = 0
    unplaced_seq = 0
    gap_bp = 0
    gap_seq = 0

    allowed_seq_types = {"A", "D", "F", "G", "O", "P", "W"}
    allowed_gap_types = {"N", "U"}

    # Get the set of placed sequences from the confidence scores file
    placed_seqs = set()
    with open(confidence_file, "r") as f:
        f.readline()  # discard header
        for line in f:
            header, g_score, l_score, o_score = line.rstrip().split("\t")
            placed_seqs.add(header)

    # Iterate through the AGP file
    with open(agp_file, "r") as f:
        for line in f:
            if not line.startswith("#"):
                obj, obj_start, obj_end, pid, ctype, comp, comp_beg, comp_end, strand = line.rstrip().split("\t")
                if ctype in allowed_seq_types:
                    seq_len = int(comp_end) - (int(comp_beg) - 1)
                    if comp in placed_seqs:
                        placed_bp += seq_len
                        placed_seq += 1
                    else:
                        unplaced_bp += seq_len
                        unplaced_seq += 1
                elif ctype in allowed_gap_types:
                    seq_len = int(comp)
                    gap_bp += seq_len
                    gap_seq += 1
                else:
                    raise ValueError("AGP is not properly formatted. Invalid component type.")

    print("placed_sequences\tplaced_bp\tunplaced_sequences\tunplaced_bp\tgap_bp\tgap_sequences")
    print("\t".join([
        str(placed_seq),
        str(placed_bp),
        str(unplaced_seq),
        str(unplaced_bp),
        str(gap_bp),
        str(gap_seq)
    ]))


if __name__ == "__main__":
    main()
