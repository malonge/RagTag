import os
import argparse

from ragoo_utilities.utilities import log
from ragoo_utilities.Aligner import Minimap2Aligner
from ragoo_utilities.Aligner import NucmerAligner


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
    parser.add_argument("-l", metavar="INT", type=int, default=10000, help="minimum alignment length to use for scaffolding [10000]")
    parser.add_argument("-q", metavar="INT", type=int, default=0, help="minimum mapping quality value for alignments. only pertains to minimap2 alignments [0]")
    parser.add_argument("-i", metavar="FLOAT", type=float, default=0.2, help="minimum grouping confidence score needed to be localized [0.2]")
    parser.add_argument("-C", action='store_true', default=False, help="write unplaced contigs individually instead of making a chr0")
    parser.add_argument("-r", action='store_true', default=False, help="infer gap pad sizes from the reference. '-g' is used when adjacent contigs overlap")
    parser.add_argument("-w", action='store_true', default=False, help="overwrite pre-existing intermediate files.")

    # Get the command line arguments and ensure all paths are absolute.
    args = parser.parse_args()
    reference_file = os.path.abspath(args.reference)
    query_file = os.path.abspath(args.query)
    output_path = args.o.replace("/", "").replace(".", "")
    exclude_file = os.path.abspath(args.e)
    min_len = args.l
    gap_size = args.g
    group_score_thresh = args.i
    skip_file = os.path.abspath(args.j)
    make_chr0 = not args.C
    infer_gaps = args.r
    overwrite_files = args.w

    # Get aligner arguments
    aligner_path = args.aligner
    aligner = aligner_path.split("/")[-1]
    if aligner.split("/")[-1] not in {'minimap2', 'nucmer'}:
        raise ValueError("Must specify either 'minimap2' or 'nucmer' (PATHs allowed) with '--aligner'.")
    mm2_params = args.mm2_params
    nucmer_params = args.nucmer_params

    # Get the current working directory and output path
    cwd = os.getcwd()
    output_path = cwd + "/" + output_path + "/"
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # Align the query to the reference
    log("Aligning the query to the reference.")
    if aligner == "minimap2":
        al = Minimap2Aligner(reference_file, query_file, aligner_path, mm2_params, output_path + "query_against_ref", in_overwrite=overwrite_files)
    else:
        al = NucmerAligner(reference_file, query_file, aligner_path, nucmer_params, output_path + "query_against_ref", in_overwrite=overwrite_files)
    al.run_aligner()


if __name__ == "__main__":
    main()

