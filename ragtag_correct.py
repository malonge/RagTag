#!/usr/bin/env python

"""
MIT License

Copyright (c) 2020 Michael Alonge <malonge11@gmail.com>

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
import sys
import math
import argparse
from collections import defaultdict

import pysam
import numpy as np
from intervaltree import IntervalTree

from ragtag_utilities.utilities import log, run_oae, get_ragtag_version
from ragtag_utilities.AlignmentReader import PAFReader
from ragtag_utilities.ContigAlignment import ContigAlignment
from ragtag_utilities.AGPFile import AGPFile
from ragtag_utilities.Aligner import Minimap2Aligner
from ragtag_utilities.Aligner import Minimap2SAMAligner
from ragtag_utilities.Aligner import NucmerAligner


def read_genome_alignments(aln_file, query_blacklist, ref_blacklist):
    tmp_ctg_alns = dict()
    aln_reader = PAFReader(aln_file)
    for aln_line in aln_reader.parse_alignments():
        # Check that the contig and reference in this alignment are allowed.
        if aln_line.query_header not in query_blacklist and aln_line.ref_header not in ref_blacklist:
            if aln_line.query_header not in tmp_ctg_alns:
                tmp_ctg_alns[aln_line.query_header] = [aln_line.query_header, aln_line.query_len,
                                                       [aln_line.query_start], [aln_line.query_end], [aln_line.strand],
                                                       [aln_line.ref_header], [aln_line.ref_len],
                                                       [aln_line.ref_start], [aln_line.ref_end],
                                                       [aln_line.num_match], [aln_line.aln_len],
                                                       [aln_line.mapq]]
            else:
                tmp_ctg_alns[aln_line.query_header][2].append(aln_line.query_start)
                tmp_ctg_alns[aln_line.query_header][3].append(aln_line.query_end)
                tmp_ctg_alns[aln_line.query_header][4].append(aln_line.strand)
                tmp_ctg_alns[aln_line.query_header][5].append(aln_line.ref_header)
                tmp_ctg_alns[aln_line.query_header][6].append(aln_line.ref_len)
                tmp_ctg_alns[aln_line.query_header][7].append(aln_line.ref_start)
                tmp_ctg_alns[aln_line.query_header][8].append(aln_line.ref_end)
                tmp_ctg_alns[aln_line.query_header][9].append(aln_line.num_match)
                tmp_ctg_alns[aln_line.query_header][10].append(aln_line.aln_len)
                tmp_ctg_alns[aln_line.query_header][11].append(aln_line.mapq)

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
            tmp_ctg_alns[i][10],
            tmp_ctg_alns[i][11]
        )
    return ctg_alns


def get_median_read_coverage(output_path, num_threads, overwrite_files):
    """ Given the read alignments, use samtools stats to return an approximate median coverage value. """
    log("Calculating global read coverage")
    if os.path.isfile(output_path + "c_reads_against_query.s.bam.stats"):
        if not overwrite_files:
            log("retaining pre-existing file: " + output_path + "c_reads_against_query.s.bam.stats")
        else:
            log("overwriting pre-existing file: " + output_path + "c_reads_against_query.s.bam.stats")
            st = pysam.stats("-@", str(num_threads), output_path + "c_reads_against_query.s.bam")
            with open(output_path + "c_reads_against_query.s.bam.stats", "w") as f:
                f.write(st)
    else:
        st = pysam.stats("-@", str(num_threads), output_path + "c_reads_against_query.s.bam")
        with open(output_path + "c_reads_against_query.s.bam.stats", "w") as f:
            f.write(st)

    # Get the coverage histogram (for 1 to 1k)
    covs = []
    with open(output_path + "c_reads_against_query.s.bam.stats") as f:
        for line in f:
            if line.startswith("COV"):
                covs.append(int(line.split("\t")[3]))

    # Get the median from the histogram
    covs = np.asarray(covs, dtype=np.int32)

    # Remove the last value, which is a catch-all for coverages > 1k
    covs = covs[:-1]
    mid = sum(covs) // 2
    cs = 0
    for i in range(len(covs)):
        cs += covs[i]
        if cs >= mid:
            return i
    raise ValueError("Unable to calculate read coverage. Check SAM/BAM files and stats file.")


def run_samtools(output_path, num_threads, overwrite_files):
    """ Compress, sort and index alignments with pysam. """
    if os.path.isfile(output_path + "c_reads_against_query.s.bam"):
        if not overwrite_files:
            log("Retaining pre-existing file: " + output_path + "c_reads_against_query.s.bam")
        else:
            log("Overwriting pre-existing file: " + output_path + "c_reads_against_query.s.bam")
            pysam.view("-@", str(num_threads), "-b", "-o", output_path + "c_reads_against_query.bam", output_path + "c_reads_against_query.sam", catch_stdout=False)
            pysam.sort("-@", str(num_threads), "-o", output_path + "c_reads_against_query.s.bam", output_path + "c_reads_against_query.bam", catch_stdout=False)
    else:
        pysam.view("-@", str(num_threads), "-b", "-o", output_path + "c_reads_against_query.bam", output_path + "c_reads_against_query.sam", catch_stdout=False)
        pysam.sort("-@", str(num_threads), "-o", output_path + "c_reads_against_query.s.bam", output_path + "c_reads_against_query.bam", catch_stdout=False)

    log("Indexing read alignments")
    if os.path.isfile(output_path + "c_reads_against_query.s.bam.bai"):
        if not overwrite_files:
            log("Retaining pre-existing file: " + output_path + "c_reads_against_query.s.bam.bai")
        else:
            log("Overwriting pre-existing file: " + output_path + "c_reads_against_query.s.bam.bai")
            pysam.index(output_path + "c_reads_against_query.s.bam", catch_stdout=False)
    else:
        pysam.index(output_path + "c_reads_against_query.s.bam", catch_stdout=False)


def clean_breaks(val_breaks, d):
    """ Merge breakpoints that are within d bp of each other. """
    breaks = sorted(list(set(val_breaks)))
    i, j = 0, 1
    while j < len(breaks):
        if breaks[j] - breaks[i] < d:
            breaks.pop(j)
        else:
            i += 1
            j += 1
    return breaks


def validate_breaks(ctg_breaks, output_path, num_threads, overwrite_files, min_break_end_dist, max_cutoff, min_cutoff, window_size=10000, num_devs=3, clean_dist=1000, debug=False):
    """
    """
    # Get the median coverage over all bp
    glob_med = get_median_read_coverage(output_path, num_threads, overwrite_files)
    dev = round(math.sqrt(glob_med))

    if max_cutoff == -1:
        max_cutoff = glob_med + (num_devs*dev)

    if min_cutoff == -1:
        min_cutoff = max(0, (glob_med - (num_devs*dev)))

    log("The global median read coverage is %dX" % glob_med)
    log("The max and min coverage thresholds are %dX and %dX, respectively" % (max_cutoff, min_cutoff))

    # Go through each break point and query the coverage within the vicinity of the breakpoint.
    bam = pysam.AlignmentFile(output_path + "c_reads_against_query.s.bam")
    validated_ctg_breaks = dict()
    for ctg in ctg_breaks:
        val_breaks = []

        # Iterate over each breakpoint for this query sequence
        for b in ctg_breaks[ctg]:
            # Don't extend the validation window too close to the contig ends (defined by min_break_end_dist)
            min_range = max(min_break_end_dist, b - (window_size//2))
            max_range = min((bam.get_reference_length(ctg) - min_break_end_dist), b + (window_size // 2))

            if min_range >= max_range:
                continue

            region = "%s:%d-%d" % (ctg, min_range, max_range-1)
            depth_out = pysam.samtools.depth("-aa", "-r", region, output_path + "c_reads_against_query.s.bam")
            covs = np.asarray(
                [j.split("\t")[2] for j in [i for i in depth_out.rstrip().split("\n")]],
                dtype=np.int32
            )
            assert len(covs) == max_range - min_range

            # Given the coverage in vicinity of the breakpoint, find the max and min coverage.
            cov_min, cov_max = np.min(covs), np.max(covs)
            too_high = True if cov_max >= max_cutoff else False
            too_low = True if cov_min <= min_cutoff else False
            new_break = None
            status = "not validated"
            if too_low and too_high:
                val_breaks.append(np.argmin(covs) + min_range)
                new_break = np.argmin(covs) + min_range
                status = "low and high cov"
            elif too_low:
                val_breaks.append(np.argmin(covs) + min_range)
                new_break = np.argmin(covs) + min_range
                status = "low cov"
            elif too_high:
                val_breaks.append(np.argmax(covs) + min_range)
                new_break = np.argmax(covs) + min_range
                status = "high cov"

            if debug:
                log("query: %s, original break: %s, window start: %d, window end: %d, status: %s, new_break: %s, cov max: %d, cov min: %d" %(ctg, b, min_range, max_range, status, str(new_break), cov_max, cov_min))

        validated_ctg_breaks[ctg] = clean_breaks(val_breaks, clean_dist)

    return validated_ctg_breaks


def make_gff_interval_tree(gff_file):
    # Dictionary storing an interval tree for each sequence header
    t = defaultdict(IntervalTree)

    # Iterate over the gff file
    with open(gff_file, "r") as f:
        for line in f:
            if not line.startswith("#"):
                fields = line.split("\t")
                h, start, end = fields[0], int(fields[3]), int(fields[4])
                start = start - 1  # make everything zero-indexed
                assert start < end

                if end - start > 100000:
                    coords = "%s:%d-%d" %(h, start+1, end)
                    log("WARNING: large interval in this gff file (%s). This could disproportionately invalidate putative query breakpoints." % coords)
                t[h][start:end] = (start, end)

    return t


def write_breaks(out_file, query_file, ctg_breaks, overwrite, remove_suffix):
    """ Write the intermediate file for contig breaks in AGP v2.1 format."""
    # Check if the output file already exists
    if os.path.isfile(out_file):
        if not overwrite:
            log("Retaining pre-existing file: " + out_file)
            return

        else:
            log("Overwriting pre-existing file: " + out_file)

    fai = pysam.FastaFile(query_file)
    all_q_seqs = sorted(fai.references)
    agp = AGPFile(out_file, mode="w")

    agp.add_pragma()
    agp.add_comment("# AGP created by RagTag {}".format(get_ragtag_version()))

    for q in all_q_seqs:

        # Check if this sequence was broken during misassembly correction
        if q not in ctg_breaks:

            # Add suffix to query header, unless otherwise requested
            unchanged_comp_header = q
            if not remove_suffix:
                unchanged_comp_header = q + ":0" + "-" + str(fai.get_reference_length(q)) + "(+)"

            agp.add_seq_line(
                    q,
                    "1",
                    str(fai.get_reference_length(q)),
                    "1",
                    "W",
                    unchanged_comp_header,
                    "1",
                    str(fai.get_reference_length(q)),
                    "+"
            )
        else:  # This query sequence was broken
            pid = 1
            sorted_breaks = sorted(ctg_breaks[q])
            start = 0
            for i in sorted_breaks:
                agp.add_seq_line(
                        q,
                        str(start+1),
                        str(i),
                        str(pid),
                        "W",
                        q + ":" + str(start) + "-" + str(i) + "(+)",
                        "1",
                        str(i-start),
                        "+"
                )
                start = i
                pid += 1

            # Add one line for the last interval
            agp.add_seq_line(
                    q,
                    str(start+1),
                    str(fai.get_reference_length(q)),
                    str(pid),
                    "W",
                    q + ":" + str(start) + "-" + str(fai.get_reference_length(q)) + "(+)",
                    "1",
                    str(fai.get_reference_length(q)-start),
                    "+"
            )

    log("Writing: " + out_file)
    agp.write()
    fai.close()


def main():
    parser = argparse.ArgumentParser(description='Reference-guided misassembly correction', usage="ragtag.py correct <reference.fa> <query.fa>")

    parser.add_argument("reference", metavar="<reference.fa>", nargs='?', default="", type=str, help="reference fasta file (uncompressed or bgzipped)")
    parser.add_argument("query", metavar="<query.fa>", nargs='?', default="", type=str, help="query fasta file (uncompressed or bgzipped)")

    cor_options = parser.add_argument_group("correction options")
    cor_options.add_argument("-f", metavar="INT", type=int, default=1000, help="minimum unique alignment length [1000]")
    cor_options.add_argument("--remove-small", action="store_true", default=False, help="remove unique alignments shorter than -f")
    cor_options.add_argument("-q", metavar="INT", type=int, default=10, help="minimum mapq (NA for Nucmer alignments) [10]")
    cor_options.add_argument("-d", metavar="INT", type=int, default=100000, help="alignment merge distance [100000]")
    cor_options.add_argument("-b", metavar="INT", type=int, default=5000, help="minimum break distance from contig ends [5000]")
    cor_options.add_argument("-e", metavar="<exclude.txt>", type=str, default="", help="list of reference headers to ignore [null]")
    cor_options.add_argument("-j", metavar="<skip.txt>", type=str, default="", help="list of query headers to leave uncorrected [null]")
    cor_options.add_argument("--inter", action="store_true", default=False, help="only break misassemblies between reference sequences")
    cor_options.add_argument("--intra", action="store_true", default=False, help="only break misassemblies within reference sequences")
    cor_options.add_argument("--gff", metavar="<features.gff>", type=str, default="", help="don't break sequences within gff intervals [null]")

    io_options = parser.add_argument_group("input/output options")
    io_options.add_argument("-o", metavar="PATH", type=str, default="ragtag_output", help="output directory [./ragtag_output]")
    io_options.add_argument("-w", action='store_true', default=False, help="overwrite intermediate files")
    io_options.add_argument("-u", action='store_true', default=False, help="add suffix to unaltered sequence headers")
    io_options.add_argument("--debug", action='store_true', default=False, help=argparse.SUPPRESS)

    aln_options = parser.add_argument_group("mapping options")
    mm2_default = "-x asm5"
    aln_options.add_argument("-t", metavar="INT", type=int, default=1, help="number of minimap2 threads [1]")
    aln_options.add_argument("--aligner", metavar="PATH", type=str, default="minimap2", help="whole genome aligner executable ('nucmer' or 'minimap2') [minimap2]")
    aln_options.add_argument("--mm2-params", metavar="STR", type=str, default=mm2_default, help="space delimited minimap2 whole genome alignment parameters ['%s']" % mm2_default)
    aln_options.add_argument("--nucmer-params", metavar="STR", type=str, default="-l 100 -c 500", help="space delimted nucmer whole genome alignment parameters ['-l 100 -c 500']")

    val_options = parser.add_argument_group("validation options")
    val_options.add_argument("--read-aligner", metavar="PATH", type=str, default="minimap2", help="read aligner executable (only 'minimap2' is allowed) [minimap2]")
    val_options.add_argument("-R", metavar="<reads.fasta>", type=str, default="", help="validation reads (uncompressed or gzipped) [null]")
    val_options.add_argument("-F", metavar="<reads.fofn>", type=str, default="", help="same as '-R', but a list of files [null]")
    val_options.add_argument("-T", metavar="STR", type=str, default="", help="read type. 'sr' and 'corr' accepted for short reads and error corrected long-reads, respectively [null]")
    val_options.add_argument("-v", metavar="INT", type=int, default=10000, help="coverage validation window size [10000]")
    val_options.add_argument("--max-cov", metavar="INT", type=int, default=-1, help="break sequences at regions at or above this coverage level [AUTO]")
    val_options.add_argument("--min-cov", metavar="INT", type=int, default=-1, help="break sequences at regions at or below this coverage level [AUTO]")
    val_options.add_argument("-m", metavar="INT", type=int, default=1000, help=argparse.SUPPRESS)  # Merge breakpoints within this distance after validation

    args = parser.parse_args()

    if not args.reference or not args.query:
        parser.print_help()
        print("\n** The reference and query FASTA files are required **")
        sys.exit()

    log("RagTag " + get_ragtag_version())
    log("CMD: ragtag.py correct " + " ".join(sys.argv[1:]))

    reference_file = os.path.abspath(args.reference)
    query_file = os.path.abspath(args.query)

    # Check that the reference/query file exists
    if not os.path.isfile(reference_file):
        raise ValueError("Could not find file: %s" % reference_file)

    if not os.path.isfile(query_file):
        raise ValueError("Could not find file: %s" % query_file)

    num_threads = args.t
    min_ulen = args.f
    keep_small_uniques = not args.remove_small
    merge_dist = args.d
    min_break_dist = args.m
    min_break_end_dist = args.b
    val_window_size = args.v

    # I/O options
    output_path = args.o
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = os.path.abspath(output_path) + "/"

    # Setup a log file for external RagTag scripts
    ragtag_log = output_path + "ragtag.correct.err"
    open(ragtag_log, "w").close()  # Wipe the log file

    overwrite_files = args.w
    remove_suffix = not args.u
    if remove_suffix:
        log("WARNING: Without '-u' invoked, some component/object AGP pairs might share the same ID. Some external programs/databases don't like this. To ensure valid AGP format, use '-u'.")

    gff_file = args.gff
    if gff_file:
        gff_file = os.path.abspath(gff_file)

    # Skip/exclude options
    query_blacklist = set()
    skip_file = args.j
    if skip_file:
        skip_file = os.path.abspath(args.j)
        with open(skip_file, "r") as f:
            for line in f:
                query_blacklist.add(line.rstrip())

    ref_blacklist = set()
    exclude_file = args.e
    if exclude_file:
        exclude_file = os.path.abspath(args.e)
        with open(exclude_file, "r") as f:
            for line in f:
                ref_blacklist.add(line.rstrip())

    # Get aligner arguments
    genome_aligner_path = args.aligner
    genome_aligner = genome_aligner_path.split("/")[-1]
    if genome_aligner.split("/")[-1] not in {'minimap2', 'nucmer'}:
        raise ValueError("Must specify either 'minimap2' or 'nucmer' (PATHs allowed) with '--aligner'.")

    mm2_params = args.mm2_params
    nucmer_params = args.nucmer_params

    # Mapq filtering params
    min_mapq = args.q
    if genome_aligner == "nucmer":
        min_mapq = 0

    # Add the number of mm2 threads if the mm2 params haven't been overridden.
    if mm2_params == mm2_default:
        mm2_params += " -t " + str(num_threads)

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

    # read-alignment parameters
    val_reads = args.R
    val_reads_fofn = args.F
    val_reads_tech = args.T
    read_aligner_path = args.read_aligner
    read_aligner = read_aligner_path.split("/")[-1]
    if read_aligner != "minimap2":
        raise ValueError("Only minimap2 can be used for read alignments. got: %s" % read_aligner)

    # If the genome aligner is minimap2, we can just use that path for read alignment
    if genome_aligner == 'minimap2':
        read_aligner_path = genome_aligner_path

    # Make sure that if -R or -F, -T has been specified.
    if val_reads or val_reads_fofn:
        if not val_reads_tech:
            raise ValueError("'-T' must be provided when using -R or -F.")

    # Make a list of read sequences.
    read_files = []
    if val_reads_fofn:
        with open(val_reads_fofn, "r") as f:
            for line in f:
                read_files.append(os.path.abspath(line.rstrip()))
    elif val_reads:
        read_files.append(os.path.abspath(val_reads))

    # Coverage thresholds
    max_cov = args.max_cov
    min_cov = args.min_cov

    if max_cov < 0:
        if max_cov != -1:
            raise ValueError("--max-cov must be >=0")

    if min_cov < 0:
        if min_cov != -1:
            raise ValueError("--min-cov must be >=0")

    # Debugging options
    debug_mode = args.debug
    debug_non_fltrd_file = output_path + "ragtag.correction.debug.unfiltered.paf"
    debug_fltrd_file = output_path + "ragtag.correction.debug.filtered.paf"
    debug_merged_file = output_path + "ragtag.correction.debug.merged.paf"
    debug_query_info_file = output_path + "ragtag.correction.debug.query.info.txt"

    # Align the query to the reference.
    log("Mapping the query genome to the reference genome")
    if genome_aligner == "minimap2":
        al = Minimap2Aligner(reference_file, [query_file], genome_aligner_path, mm2_params, output_path + "c_query_against_ref", in_overwrite=overwrite_files)
    else:
        al = NucmerAligner(reference_file, [query_file], genome_aligner_path, nucmer_params, output_path + "c_query_against_ref", in_overwrite=overwrite_files)
    al.run_aligner()

    # If alignments are from Nucmer, convert from delta to paf.
    if genome_aligner == "nucmer":
        cmd = ["ragtag_delta2paf.py", output_path + "c_query_against_ref.delta"]
        run_oae(cmd, output_path + "c_query_against_ref.paf", ragtag_log)

    # Read and organize the alignments.
    log('Reading whole genome alignments')
    # ctg_alns = dict :: key=query header, value=ContigAlignment object
    ctg_alns = read_genome_alignments(output_path + "c_query_against_ref.paf", query_blacklist, ref_blacklist)

    # Filter and merge the alignments.
    if debug_mode:
        # create new empty copies of debugging output files
        open(debug_non_fltrd_file, "w").close()
        open(debug_fltrd_file, "w").close()
        open(debug_merged_file, "w").close()
        open(debug_query_info_file, "w").close()

    log("Filtering and merging alignments")
    for i in ctg_alns:

        # Write unfiltered alignments
        if debug_mode:
            with open(debug_non_fltrd_file, "a") as f:
                f.write(str(ctg_alns[i]))

        ctg_alns[i] = ctg_alns[i].unique_anchor_filter(min_ulen, keep_small=keep_small_uniques)
        if ctg_alns[i] is not None:
            ctg_alns[i] = ctg_alns[i].filter_mapq(min_mapq)
            if ctg_alns[i] is not None:

                # Write filtered alignments
                if debug_mode:
                    with open(debug_fltrd_file, "a") as f:
                        f.write(str(ctg_alns[i]))

                ctg_alns[i] = ctg_alns[i].merge_alns(merge_dist=merge_dist)

    # Get the putative breakpoints for each query sequence, if any.
    ctg_breaks = dict()
    for i in ctg_alns:
        if ctg_alns[i] is not None:

            # Write merged alignments and confidence scores
            if debug_mode:
                with open(debug_merged_file, "a") as f:
                    f.write(str(ctg_alns[i]))

                with open(debug_query_info_file, "a") as f:
                    f.write("\t".join([
                        i,
                        ctg_alns[i].best_ref_header,
                        str(ctg_alns[i].grouping_confidence),
                        str(ctg_alns[i].location_confidence),
                        str(ctg_alns[i].orientation_confidence),
                    ]) + "\n")

            breaks = []
            intra_breaks, inter_breaks = ctg_alns[i].get_break_candidates(min_dist=min_break_end_dist)
            if break_intra:
                breaks = breaks + intra_breaks
            if break_inter:
                breaks = breaks + inter_breaks
            if breaks:
                ctg_breaks[i] = breaks

    # If desired, validate the putative breakpoints by observing read coverage.
    if read_files:
        log("Validating putative query breakpoints via read alignment.")
        log("Aligning reads to query sequences.")
        if not os.path.isfile(output_path + "c_reads_against_query.s.bam"):
            if val_reads_tech == "sr":
                al = Minimap2SAMAligner(query_file, read_files, read_aligner_path, "-ax sr -t " + str(num_threads),
                                        output_path + "c_reads_against_query", in_overwrite=overwrite_files)
            elif val_reads_tech == "corr":
                al = Minimap2SAMAligner(query_file, read_files, read_aligner_path, "-ax asm5 -t " + str(num_threads),
                                        output_path + "c_reads_against_query", in_overwrite=overwrite_files)
            else:
                raise ValueError("'-T' must be either 'sr' or 'corr'.")
            al.run_aligner()
        else:
            log("Retaining pre-existing read alignments: " + output_path + "c_reads_against_query.s.bam")

        # Compress, sort and index the alignments.
        log("Compressing, sorting, and indexing read alignments")
        run_samtools(output_path, num_threads, overwrite_files)

        # Validate the breakpoints
        log("Validating putative query breakpoints")
        
        # Give at least 10k/1k from ctg ends for coverage to accumulate for corr and sr, respectively.
        val_min_break_end_dist = min_break_end_dist
        if val_reads_tech == "corr":
            val_min_break_end_dist = max(10000, min_break_end_dist)
        if val_reads_tech == "sr":
            val_min_break_end_dist = max(1000, min_break_end_dist)
            
        # Validate the breakpoints
        ctg_breaks = validate_breaks(ctg_breaks, output_path, num_threads, overwrite_files, val_min_break_end_dist, max_cov, min_cov, window_size=val_window_size, clean_dist=min_break_dist, debug=debug_mode)

    # Check if we need to avoid gff intervals
    if gff_file:
        log("Avoiding breaks within GFF intervals")
        it = make_gff_interval_tree(gff_file)
        non_gff_breaks = dict()
        for ctg in ctg_breaks:
            new_breaks = []
            for i in ctg_breaks[ctg]:
                if it[ctg][i]:
                    log("Avoiding breaking %s at %d. This point intersects a feature in the gff file." % (ctg, i))
                else:
                    new_breaks.append(i)
            if new_breaks:
                non_gff_breaks[ctg] = new_breaks
        ctg_breaks = non_gff_breaks

    # Write the summary of query sequence breaks in AGP format
    agp_file = output_path + "ragtag.correction.agp"
    write_breaks(agp_file, query_file, ctg_breaks, True, remove_suffix)

    # Write the scaffolds.
    log("Writing broken contigs")
    qf_name = query_file.split("/")[-1]
    qf_pref = qf_name[:qf_name.rfind(".")]
    cmd = [
        "ragtag_break_query.py",
        agp_file,
        query_file
    ]
    run_oae(cmd, output_path + qf_pref + ".corrected.fasta", ragtag_log)

    log("Goodbye")


if __name__ == "__main__":
    main()
