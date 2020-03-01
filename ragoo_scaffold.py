import os
import argparse

from ragoo_utilities.utilities import log
from ragoo_utilities.Aligner import Minimap2Aligner
from ragoo_utilities.Aligner import NucmerAligner
from ragoo_utilities.AlignmentReader import AlignmentReader
from ragoo_utilities.ContigAlignment import ContigAlignment


def main():
    parser = argparse.ArgumentParser(description='Scaffold contigs according to alignments to a reference (v2.0.0)')
    parser.add_argument("reference", metavar="<reference.fasta>", type=str, help="reference fasta file. must not be gzipped.")
    parser.add_argument("query", metavar="<query.fasta>", type=str, help="query fasta file to be scaffolded. must not be gzipped.")
    parser.add_argument("-o", metavar="STR", type=str, default="ragoo_output", help="output directory name [ragoo_output]")
    parser.add_argument("--aligner", metavar="PATH", type=str, default="minimap2", help="Aligner ('nucmer' or 'minimap2') to use for scaffolding. PATHs allowed [minimap2]")
    parser.add_argument("--mm2-params", metavar="STR", type=str, default="-k19 -w19 -t1", help="Space delimted parameters to pass directly to minimap2 ['-k19 -w19 -t1']")
    parser.add_argument("--nucmer-params", metavar="STR", type=str, default="-l 100 -c 500", help="Space delimted parameters to pass directly to nucmer ['-l 100 -c 500']")
    parser.add_argument("-e", metavar="<exclude.txt>", type=str, default="", help="single column text file of reference headers to ignore")
    parser.add_argument("-j", metavar="<skip.txt>", type=str, default="", help="List of contigs to automatically leave unplaced")
    parser.add_argument("-g", metavar="INT", type=int, default=100, help="gap size for padding in pseudomolecules [100]")
    parser.add_argument("-l", metavar="INT", type=int, default=10000, help="minimum unique alignment length to use for scaffolding [10000]")
    parser.add_argument("-q", metavar="INT", type=int, default=40, help="minimum mapping quality value for alignments. only pertains to minimap2 alignments [0]")
    parser.add_argument("-i", metavar="FLOAT", type=float, default=0.2, help="minimum grouping confidence score needed to be localized [0.2]")
    parser.add_argument("-C", action='store_true', default=False, help="write unplaced contigs individually instead of making a chr0")
    parser.add_argument("-r", action='store_true', default=False, help="infer gap pad sizes from the reference. '-g' is used when adjacent contigs overlap")
    parser.add_argument("-w", action='store_true', default=False, help="overwrite pre-existing intermediate files.")

    # Get the command line arguments and ensure all paths are absolute.
    args = parser.parse_args()
    reference_file = os.path.abspath(args.reference)
    query_file = os.path.abspath(args.query)
    output_path = args.o.replace("/", "").replace(".", "")
    min_len = args.l
    min_q = args.q
    gap_size = args.g
    group_score_thresh = args.i
    make_chr0 = not args.C
    infer_gaps = args.r
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

    # Get the skip and exclude
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
    log("aligning the query to the reference")
    if aligner == "minimap2":
        al = Minimap2Aligner(reference_file, query_file, aligner_path, mm2_params, output_path + "query_against_ref", in_overwrite=overwrite_files)
    else:
        al = NucmerAligner(reference_file, query_file, aligner_path, nucmer_params, output_path + "query_against_ref", in_overwrite=overwrite_files)
    al.run_aligner()

    # Read and organize the alignments
    log('reading alignments')
    ctg_alns = dict()
    aln_reader = AlignmentReader(output_path + "query_against_ref", aligner)
    for aln_line in aln_reader.parse_alignments():

        ## Check that the contig and reference in this alignment are allowed.
        if aln_line.query_header not in query_blacklist and aln_line.ref_header not in ref_blacklist:
            if aln_line.query_header not in ctg_alns:
                ctg_alns[aln_line.query_header] = ContigAlignment(aln_line.query_header, aln_line.query_len, [], [], [], [], [], [], [], [])
            ctg_alns[aln_line.query_header] = ctg_alns[aln_line.query_header].add_alignment(aln_line.ref_header, aln_line.ref_len, aln_line.ref_start, aln_line.ref_end, aln_line.query_start, aln_line.query_end, aln_line.strand, aln_line.mapq)

    # Filter the alignments
    log("filtering alignments")
    if aligner == "minimap2":
        log('alignments are from minimap2. removing alignments with mapq < %r.' % min_q)
        for i in ctg_alns:
            ctg_alns[i] = ctg_alns[i].filter_mapq(min_q)
            ctg_alns[i] = ctg_alns[i].unique_anchor_filter(min_len)
    else:
        for i in ctg_alns:
            ctg_alns[i] = ctg_alns[i].unique_anchor_filter(min_len)


    # order and orient scaffolds


if __name__ == "__main__":
    main()

