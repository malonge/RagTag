#!/usr/bin/env python

"""
MIT License

Copyright (c) 2021 Michael Alonge <malonge11@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import os
import re
import sys
import argparse

import pysam

from ragtag_correct import make_gff_interval_tree

from ragtag_utilities.utilities import get_ragtag_version
from ragtag_utilities.utilities import log
from ragtag_utilities.AGPFile import AGPFile


def main():
    parser = argparse.ArgumentParser(description='Split sequences at gaps', usage="ragtag.py splitasm <asm.fa>")
    parser.add_argument("asm", metavar="<asm.fa>", default="", type=str, help="assembly fasta file (uncompressed or bgzipped)")
    parser.add_argument("-n", metavar="INT", type=int, default=0, help="minimum gap size [0]")
    parser.add_argument("-o", metavar="PATH", type=str, default="ragtag.splitasm.agp", help="output AGP file path [./ragtag.splitasm.agp]")
    parser.add_argument("--gff", metavar="<features.gff>", type=str, default="", help="don't break sequences within gff intervals [null]")
    parser.add_argument("--bed", metavar="<splits.bed>", type=str, default="", help="BED file of where to break sequences [null]")

    # Parse the command line arguments
    args = parser.parse_args()
    if not args.asm:
        parser.print_help()
        print("\n** The assembly FASTA file is required **")
        sys.exit()

    asm_fn = args.asm
    min_gap_size = args.n
    agp_fn = args.o
    gff_file = args.gff
    bed_file = args.bed
    # Initialize Tabix-indexed BED file
    if bed_file:
        log("INFO", "Breaking across positions from BED file")
        bed = pysam.TabixFile(bed_file)
        bed_contigs = [row.contig for row in bed.fetch(parser = pysam.asBed())]
        #print(bed_contigs)
        # Fetch will break if sequence is absent
        #for row in bed.fetch("seq00000000", parser = pysam.asBed()):
        #    print(row)
            #log("INFO", row.contig)
    if gff_file:
        gff_file = os.path.abspath(gff_file)
        it = make_gff_interval_tree(gff_file)
        log("INFO", "Avoiding breaks within GFF intervals")
    # Initialize the AGP file
    agp = AGPFile(agp_fn, mode="w")
    agp.add_pragma()
    agp.add_comment("# AGP created by RagTag {}".format(get_ragtag_version()))



    # Process the FASTA file
    new_header_idx = 0
    fai = pysam.FastaFile(asm_fn)
    for header in sorted(fai.references):
        seq = fai.fetch(header).upper()
        seq_len = fai.get_reference_length(header)
        gap_coords = [(i.start(), i.end()) for i in re.finditer(r'N+', seq) if i.end() - i.start() > min_gap_size]
        if bed_file and header in bed_contigs:
            splits_pos = [(row.start, row.end) for row in bed.fetch(header, parser=pysam.asBed()) if row.end <= seq_len]
        else:
            splits_pos = []
        #log("INFO", "Splitting %s at %s" % (header, splits_pos))
        #print("Splitting %s at %s" % (header, splits_pos))
        # Remove coordinates overlapping gff features
        if gff_file:
            #non_gff_breaks = dict()
            new_breaks = []
            for gap in gap_coords:
                if it[header][gap[0]] or it[header][gap[1]]:
                    log("INFO", "Avoiding breaking %s between %d-%d. This point intersects a feature in the gff file."
                                % (header, gap[0], gap[1]))
                else:
                    new_breaks.append(gap)
            if splits_pos:
                new_splits = []
                for split in splits_pos:
                    if it[header][split[1]]:
                        log("INFO", "Avoiding breaking %s at %d. This point intersects a feature in the gff file."
                                    % (header, split[1]))
                    else:
                        new_splits.append(split)
                if new_splits:
                    splits_pos = new_splits
            if new_breaks:
                #non_gff_breaks[header] = new_breaks
                gap_coords = new_breaks
            

        if not gap_coords and not splits_pos:
            new_header = "seq{0:08}".format(new_header_idx)
            new_header_idx += 1
            agp.add_seq_line(header, "1", seq_len, "1", "W", new_header, "1", seq_len, "+")
        else:
            gap_coords.append((seq_len, seq_len + 1))
            pid = 1
            if gap_coords[0][0]:
                # The sequence doesn't start with a gap
                new_header = "seq{0:08}".format(new_header_idx)
                # if there's a split coordinate before the end, that needs to be the "new" end
                # First, if/ELSE allows current behavior to be the same
                # Make any/all splits necessary before the *first* gap sequence
                # *later*: will need to check for split coords between every gap coord
                if splits_pos:
                    log("INFO", "Should break %s %s times at gaps %s" % (header, len(splits_pos), str(splits_pos)))
                    for i in range(0, len(splits_pos)):
                        log("INFO", "Breaking %s at gap %s" % (header, str(splits_pos)))
                        if splits_pos[i][1] < gap_coords[0][0]:
                            # Only first iteration starts at 1 (object; component always starts at 1)
                            if i == 0:
                                obj_start=1
                                obj_end=splits_pos[i][1]
                                cmp_end=obj_end
                                log("INFO", "Found first split in %s at %s" % (header, obj_end))
                            else:
                                obj_start=obj_end+1
                                obj_end=splits_pos[i][1]
                                cmp_end=obj_end - obj_start + 1
                                log("INFO", "Found another split in %s at %s producing seq of length %s" % (header, obj_end, cmp_end))
                            agp.add_seq_line(header, str(obj_start), str(obj_end), str(pid), "W", new_header, "1", str(cmp_end), "+")
                            new_header_idx += 1
                            new_header = "seq{0:08}".format(new_header_idx)
                            pid += 1
                    # Add in final component to object - from last break till end of the sequence
                    obj_start=obj_end+1
                    obj_end=seq_len
                    cmp_end=obj_end - obj_start + 1
                    agp.add_seq_line(header, str(obj_start), str(obj_end), str(pid), "W", new_header, "1", str(cmp_end), "+")
                else:  
                    agp.add_seq_line(header, "1", str(gap_coords[0][0]), str(pid), "W", new_header, "1", str(gap_coords[0][0]), "+")
                new_header_idx += 1
                pid += 1

            for i in range(1, len(gap_coords)):
                # Add the gap line
                gap_start, gap_end = gap_coords[i-1][0], gap_coords[i-1][1]
                gap_len = gap_end - gap_start
                agp.add_gap_line(header, str(gap_start + 1), str(gap_end), str(pid), "N", str(gap_len), "scaffold", "yes", "align_genus")
                pid += 1

                # Add the sequence line
                obj_start, obj_end = gap_coords[i-1][1], gap_coords[i][0]
                comp_len = obj_end - obj_start
                new_header = "seq{0:08}".format(new_header_idx)
                if gap_coords[i-1][1] != seq_len:
                    agp.add_seq_line(header, str(obj_start + 1), obj_end, pid, "W", new_header, "1", str(comp_len), "+")
                    new_header_idx += 1
                    pid += 1

        agp.write()

    # Iterate over the AGP file and print the sequences
    agp = AGPFile(agp_fn, mode="r")
    for line in agp.iterate_lines():
        if not line.is_gap:
            obj, comp, obj_beg, obj_end = line.obj, line.comp, line.obj_beg, line.obj_end
            print(">" + comp)
            print(fai.fetch(obj, obj_beg-1, obj_end))

    fai.close()


if __name__ == "__main__":
    main()
