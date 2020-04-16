#!/usr/bin/env python

import os
import argparse

import pysam

from ragoo2_utilities.utilities import log, run
from ragoo2_utilities.Aligner import Minimap2Aligner
from ragoo2_utilities.Aligner import Minimap2SAMAligner
from ragoo2_utilities.Aligner import NucmerAligner
from ragoo2_utilities.AlignmentReader import AlignmentReader
from ragoo2_utilities.ContigAlignment import ContigAlignment


def write_breaks(query_file, ctg_breaks, overwrite, out_path):
    """
    Write the intermediate file for contig breaks.
    This should be the same format as the intermediate output from 'ragoo_scaffold.py'. As a result,
    the same lift-over script could be used for either.
    """
    out_file = out_path + "correct.bed"

    # Check if the output file already exists
    if os.path.isfile(out_file):
        if not overwrite:
            log("retaining pre-existing file: " + out_file)
            return

    x = pysam.FastaFile(query_file)
    all_q_seqs = set(x.references)
    all_out_lines = []

    for q in all_q_seqs:
        if q not in ctg_breaks:
            # This query sequence was not broken
            all_out_lines.append(
                "\t".join([
                    q,
                    "0",
                    str(x.get_reference_length(q)),
                    "s",
                    q,
                    "+"
                ])
            )
        else:
            # This query sequence was broken
            sorted_breaks = sorted(ctg_breaks[q])
            start = 0
            for i in sorted_breaks:
                all_out_lines.append(
                    "\t".join([
                        q,
                        str(start),
                        str(i),
                        "s",
                        q + ":" + str(start) + "-" + str(i) + "(+)",
                        "+"
                    ])
                )
                start = i

            # Add one line for the last interval
            all_out_lines.append(
                "\t".join([
                    q,
                    str(start),
                    str(x.get_reference_length(q)),
                    "s",
                    q + ":" + str(start) + "-" + str(x.get_reference_length(q)) + "(+)",
                    "+"
                ])
            )

    log("Writing: " + out_file)
    with open(out_file, "w") as f:
        f.write("\n".join(all_out_lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description='Correct contigs according to alignments to a reference (v2.0.0)')
    parser.add_argument("reference", metavar="<reference.fasta>", type=str, help="reference fasta file. must not be gzipped.")
    parser.add_argument("query", metavar="<query.fasta>", type=str, help="query fasta file to be scaffolded. must not be gzipped.")
    parser.add_argument("-o", metavar="STR", type=str, default="ragoo_correct_output", help="output directory name [ragoo_correct_output]")
    parser.add_argument("--aligner", metavar="PATH", type=str, default="minimap2", help="Aligner ('nucmer' or 'minimap2') to use for correction. PATHs allowed [minimap2]")
    parser.add_argument("--mm2-params", metavar="STR", type=str, default="-k19 -w19 -t3", help="Space delimted parameters to pass directly to minimap2 ['-k19 -w19 -t3']")
    parser.add_argument("--nucmer-params", metavar="STR", type=str, default="-l 100 -c 500", help="Space delimted parameters to pass directly to nucmer ['-l 100 -c 500']")
    parser.add_argument("-f", metavar="INT", type=int, default=10000, help="minimum unique alignment length to use for scaffolding [10000]")
    parser.add_argument("-q", metavar="INT", type=int, default=-1, help="minimum mapping quality value for alignments. only pertains to minimap2 alignments [-1]")
    parser.add_argument("-d", metavar="INT", type=int, default=100000, help="merge contig alignments within this distance [100000]")
    parser.add_argument("-b", metavar="INT", type=int, default=5000, help="don't break contigs within this distance to the contigs ends")
    parser.add_argument("-e", metavar="<exclude.txt>", type=str, default="", help="single column text file of reference headers to ignore")
    parser.add_argument("-j", metavar="<skip.txt>", type=str, default="", help="List of contigs to automatically leave uncorrected")
    parser.add_argument("--inter", action="store_true", default=False, help="Only break misassemblies between reference sequences")
    parser.add_argument("--intra", action="store_true", default=False, help="Only break misassemblies within reference sequences")
    parser.add_argument("--gff", metavar="<features.gff>", type=str, default="", help="Avoid breaking query sequences within any specified gff interval")
    parser.add_argument("-w", action='store_true', default=False,help="overwrite pre-existing intermediate files. ragoo.fasta will always be overwritten")

    # If one uses the same output directory name for both correct and scaffold, could that be ok without having to worry
    # about re-writing files? e.g.
    # ragoo_correct.py -o test_out
    # ragoo_scaffold.py -o test_out
    # Would that be a problem?
    # Warning message for large gff intervals

    """
    NOTES ON COVERAGE VALIDATION
    
    1. use minimap2 to align reads and produce a SAM file, just as v1.1
    2. compress, sort and index the file with pysam
    (for 1 and 2, check if files already exist).
    for each putative region, pull out reads mapping within a window around the breakpoint (probably 10k by default)
    use AlignmentFile.pileup for this
    use truncate=True
    double check that 0-coverage bases in the region are reported
    
    for i in x.pileup("NC_003070.9", 0, 1000, truncate=True):
        print(i.n) # give the coverage of the first 1k bases.
    
    """

    args = parser.parse_args()
    reference_file = os.path.abspath(args.reference)
    query_file = os.path.abspath(args.query)
    output_path = args.o.replace("/", "").replace(".", "")
    min_ulen = args.f
    min_q = args.q
    merge_dist = args.d
    min_break_end_dist = args.b
    overwrite_files = args.w

    skip_file = args.j
    if skip_file:
        skip_file = os.path.abspath(args.j)

    exclude_file = args.e
    if exclude_file:
        exclude_file = os.path.abspath(args.e)

    # Get aligner arguments
    aligner_path = args.aligner
    aligner = aligner_path.split("/")[-1]
    if aligner.split("/")[-1] not in {'minimap2', 'nucmer'}:
        raise ValueError("Must specify either 'minimap2' or 'nucmer' (PATHs allowed) with '--aligner'.")
    mm2_params = args.mm2_params
    nucmer_params = args.nucmer_params

    # Make sure no quality filtering for nucmer alignments
    if aligner == "nucmer":
        min_q = -1

    # Check if intra/inter breaking is desired
    break_intra = True
    break_inter = True
    only_intra = args.intra
    only_inter = args.inter
    if only_intra and only_inter:
        raise ValueError("Must speficity either '--inter' or '--intra', not both.")

    if only_intra:
        break_inter = False
    if only_inter:
        break_intra = False

    # Get the skip and exclude sets
    query_blacklist = set()
    if skip_file:
        with open(skip_file, "r") as f:
            for line in f:
                query_blacklist.add(line.rstrip())

    ref_blacklist = set()
    if exclude_file:
        with open(exclude_file, "r") as f:
            for line in f:
                ref_blacklist.add(line.rstrip())

    # Get the current working directory and output path
    cwd = os.getcwd()
    output_path = cwd + "/" + output_path + "/"
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # Align the query to the reference
    log("Aligning the query to the reference")
    if aligner == "minimap2":
        al = Minimap2Aligner(reference_file, query_file, aligner_path, mm2_params, output_path + "c_query_against_ref",
                             in_overwrite=overwrite_files)
    else:
        al = NucmerAligner(reference_file, query_file, aligner_path, nucmer_params, output_path + "c_query_against_ref",
                           in_overwrite=overwrite_files)
    al.run_aligner()

    # If alignments are from Nucmer, need to convert from delta to paf
    if aligner == "nucmer":
        cmd = ["delta2paf.py", output_path + "c_query_against_ref.delta", ">", output_path + "c_query_against_ref.paf"]
        run(" ".join(cmd))

    # Read and organize the alignments
    log('Reading alignments')
    # Read all of the alignments into a temporary structure to save time on ContigAlignment instantiation
    tmp_ctg_alns = dict()
    aln_reader = AlignmentReader(output_path + "c_query_against_ref.paf")
    for aln_line in aln_reader.parse_alignments():
        # Check that the contig and reference in this alignment are allowed.
        if aln_line.query_header not in query_blacklist and aln_line.ref_header not in ref_blacklist:
            if aln_line.query_header not in tmp_ctg_alns:
                tmp_ctg_alns[aln_line.query_header] = [aln_line.query_header, aln_line.query_len,
                                                       [aln_line.ref_header], [aln_line.ref_len],
                                                       [aln_line.ref_start], [aln_line.ref_end],
                                                       [aln_line.query_start], [aln_line.query_end],
                                                       [aln_line.strand], [aln_line.aln_len],
                                                       [aln_line.mapq]]
            else:
                tmp_ctg_alns[aln_line.query_header][2].append(aln_line.ref_header)
                tmp_ctg_alns[aln_line.query_header][3].append(aln_line.ref_len)
                tmp_ctg_alns[aln_line.query_header][4].append(aln_line.ref_start)
                tmp_ctg_alns[aln_line.query_header][5].append(aln_line.ref_end)
                tmp_ctg_alns[aln_line.query_header][6].append(aln_line.query_start)
                tmp_ctg_alns[aln_line.query_header][7].append(aln_line.query_end)
                tmp_ctg_alns[aln_line.query_header][8].append(aln_line.strand)
                tmp_ctg_alns[aln_line.query_header][9].append(aln_line.aln_len)
                tmp_ctg_alns[aln_line.query_header][10].append(aln_line.mapq)

    # Put the tmp alignments into ContigAlignment objects
    ctg_alns = dict()
    for i in tmp_ctg_alns:
        ctg_alns[i] = ContigAlignment(
            tmp_ctg_alns[i][0],
            tmp_ctg_alns[i][1],
            tmp_ctg_alns[i][2],
            tmp_ctg_alns[i][3],
            tmp_ctg_alns[i][4],
            tmp_ctg_alns[i][5],
            tmp_ctg_alns[i][6],
            tmp_ctg_alns[i][7],
            tmp_ctg_alns[i][8],
            tmp_ctg_alns[i][9],
            tmp_ctg_alns[i][10]
        )

    # Filter and merge the alignments
    log("Filtering and merging alignments")
    for i in ctg_alns:
        ctg_alns[i] = ctg_alns[i].filter_mapq(min_q)
        if ctg_alns[i] is not None:
            ctg_alns[i] = ctg_alns[i].unique_anchor_filter(min_ulen)
            if ctg_alns[i] is not None:
                ctg_alns[i] = ctg_alns[i].merge_alns(merge_dist=merge_dist)

    # Get the putative breakpoints for each query sequence, if any
    ctg_breaks = dict()
    for i in ctg_alns:
        if ctg_alns[i] is not None:
            breaks = []
            intra_breaks, inter_breaks = ctg_alns[i].get_break_candidates(min_dist=min_break_end_dist)
            if break_intra:
                breaks = breaks + intra_breaks
            if break_inter:
                breaks = breaks + inter_breaks
            if breaks:
                ctg_breaks[i] = breaks

    write_breaks(query_file, ctg_breaks, overwrite_files, output_path)

    # Write the scaffolds
    log("Writing broken contigs")
    qf_name = query_file.split("/")[-1]
    qf_pref = qf_name[:qf_name.rfind(".")]
    cmd = [
        "break_query.py",
        output_path + "correct.bed",
        query_file,
        output_path + qf_pref + ".break.fasta"
    ]
    run(" ".join(cmd))


if __name__ == "__main__":
    main()
