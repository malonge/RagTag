#!/usr/bin/env python

import os
import math
import argparse
from collections import defaultdict

import pysam
import numpy as np
from intervaltree import IntervalTree
import matplotlib.pyplot as plt

from ragoo2_utilities.utilities import log, run
from ragoo2_utilities.Aligner import Minimap2Aligner
from ragoo2_utilities.Aligner import Minimap2SAMAligner
from ragoo2_utilities.Aligner import NucmerAligner
from ragoo2_utilities.AlignmentReader import AlignmentReader
from ragoo2_utilities.ContigAlignment import ContigAlignment


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
    raise ValueError()


def run_samtools(output_path, num_threads, overwrite_files):
    """ Compress, sort and index alignments with pysam. """
    log("Compressing and sorting read alignments")
    if os.path.isfile(output_path + "c_reads_against_query.s.bam"):
        if not overwrite_files:
            log("retaining pre-existing file: " + output_path + "c_reads_against_query.s.bam")
        else:
            log("overwriting pre-existing file: " + output_path + "c_reads_against_query.s.bam")
            pysam.view("-@", str(num_threads), "-b", "-o", output_path + "c_reads_against_query.bam", output_path + "c_reads_against_query.sam", catch_stdout=False)
            pysam.sort("-@", str(num_threads), "-o", output_path + "c_reads_against_query.s.bam", output_path + "c_reads_against_query.bam", catch_stdout=False)
    else:
        pysam.view("-@", str(num_threads), "-b", "-o", output_path + "c_reads_against_query.bam", output_path + "c_reads_against_query.sam", catch_stdout=False)
        pysam.sort("-@", str(num_threads), "-o", output_path + "c_reads_against_query.s.bam", output_path + "c_reads_against_query.bam", catch_stdout=False)

    log("Indexing read alignments")
    if os.path.isfile(output_path + "c_reads_against_query.s.bam.bai"):
        if not overwrite_files:
            log("retaining pre-existing file: " + output_path + "c_reads_against_query.s.bam.bai")
        else:
            log("overwriting pre-existing file: " + output_path + "c_reads_against_query.s.bam.bai")
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


def validate_breaks(ctg_breaks, output_path, num_threads, overwrite_files, window_size=10000, num_devs=3, clean_dist=10000):
    """
    """
    # Get the median coverage over all bp
    glob_med = get_median_read_coverage(output_path, num_threads, overwrite_files)
    dev = int(math.sqrt(glob_med))
    max_cutoff = glob_med + (num_devs*dev)
    min_cutoff = max(0, (glob_med - (num_devs*dev)))
    log("The global median read coverage is %d" % glob_med)
    #print("max and min cutoff: %f, %f" % (max_cutoff, min_cutoff))

    # Go through each break point and query the coverage within the vicinity of the breakpoint.
    bam = pysam.AlignmentFile(output_path + "c_reads_against_query.s.bam")
    validated_ctg_breaks = dict()
    for ctg in ctg_breaks:
        val_breaks = []

        # Iterate over each breakpoint for this query sequence
        for b in ctg_breaks[ctg]:
            min_range = max(0, b - (window_size//2))
            max_range = min(bam.get_reference_length(ctg), b + (window_size // 2))
            region = "%s:%d-%d" % (ctg, min_range, max_range-1)
            depth_out = pysam.samtools.depth("-aa", "-r", region, output_path + "c_reads_against_query.s.bam")
            covs = np.asarray(
                [j.split("\t")[2] for j in [i for i in depth_out.rstrip().split("\n")]],
                dtype=np.int32
            )
            assert len(covs) == max_range - min_range

            # Given the coverage in vicinity of the breakpoint, find the max and min coverage.
            cov_min, cov_max = np.min(covs), np.max(covs)
            too_high = True if cov_max > max_cutoff else False
            too_low = True if cov_min < min_cutoff else False
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

            # Temp logging
            #print("Contig: %s, break: %s, status: %s, new_break: %s, cov max: %d, cov min: %d" %(ctg, b, status, str(new_break), cov_max, cov_min))

            ####################
            ## TEMP PLOTTING ###
            ####################


            #plt.figure(figsize=(10, 10))
            #plt.plot(covs)
            #plt.axvline(b-min_range, color="r", linestyle="--")
            #plt.axvline(new_break-min_range, color="g", linestyle="--")
            #plt.savefig(output_path + "plots/" + ctg + "_" + str(b) + ".png")
            #plt.close()


            ####################
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
                    log("WARNING: there are large intervals in this gff file (>= 100 kbp). This could disproportionately invalidate putative query breakpoints.")
                t[h][start:end] = (start, end)

    return t


def write_breaks(out_file, query_file, ctg_breaks, overwrite, out_path):
    """
    Write the intermediate file for contig breaks.
    This should be the same format as the intermediate output from 'ragoo_scaffold.py'. As a result,
    the same lift-over script could be used for either.
    """
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
    parser = argparse.ArgumentParser(description='Correct contigs according to alignments to a reference')
    parser.add_argument("reference", metavar="<reference.fasta>", type=str, help="reference fasta file. must not be gzipped.")
    parser.add_argument("query", metavar="<query.fasta>", type=str, help="query fasta file to be scaffolded. must not be gzipped.")
    parser.add_argument("-o", metavar="STR", type=str, default="ragoo_output", help="output directory name [ragoo_correct_output]")
    parser.add_argument("-t", metavar="INT", type=int, default=1, help="number of threads to use when running minimap2 [1]")
    parser.add_argument("--genome-aligner", metavar="PATH", type=str, default="minimap2", help="whole genome aligner ('nucmer' or 'minimap2') to use for correction. PATHs allowed [minimap2]")
    parser.add_argument("--mm2-params", metavar="STR", type=str, default="-k19 -w19", help="Space delimted parameters to pass directly to minimap2 for whole genome alignment ['-k19 -w19 -t3']")
    parser.add_argument("--nucmer-params", metavar="STR", type=str, default="-l 100 -c 500", help="Space delimted parameters to pass directly to nucmer for whole genome alignment ['-l 100 -c 500']")
    parser.add_argument("-f", metavar="INT", type=int, default=10000, help="minimum unique alignment length to use for scaffolding [10000]")
    parser.add_argument("-q", metavar="INT", type=int, default=-1, help="minimum mapping quality value for alignments. only pertains to minimap2 alignments [-1]")
    parser.add_argument("-d", metavar="INT", type=int, default=100000, help="merge contig alignments within this distance [100000]")
    parser.add_argument("-b", metavar="INT", type=int, default=5000, help="don't break contigs within this distance to the contigs ends [5000]")
    parser.add_argument("-v", metavar="INT", type=int, default=10000, help="for each putative breakpoint, query the coverage in a window this size around the breakpoint [10000]")
    parser.add_argument("-m", metavar="INT", type=int, default=10000, help="merge query breakpoints within this distance of each other [10000]")
    parser.add_argument("-e", metavar="<exclude.txt>", type=str, default="", help="single column text file of reference headers to ignore")
    parser.add_argument("-j", metavar="<skip.txt>", type=str, default="", help="List of contigs to automatically leave uncorrected")
    parser.add_argument("--inter", action="store_true", default=False, help="Only break misassemblies between reference sequences")
    parser.add_argument("--intra", action="store_true", default=False, help="Only break misassemblies within reference sequences")
    parser.add_argument("--gff", metavar="<features.gff>", type=str, default="", help="Avoid breaking query sequences within any specified gff interval")
    parser.add_argument("-w", action='store_true', default=False,help="overwrite pre-existing intermediate files. ragoo.fasta will always be overwritten")
    parser.add_argument("--read-aligner", metavar="PATH", type=str, default="minimap2", help="read aligner (only 'minimap2' is allowed). use this to specify minimap2 path if using Nucmer as the genome aligner.")
    parser.add_argument("-R", metavar="<reads.fasta>", type=str, default="", help="Align provided reads to the contigs to validate misassembly correction breakpoints. gzipped fastq or fasta allowed.")
    parser.add_argument("-F", metavar="<reads.fofn>", type=str, default="", help="same as '-R', but a list of files.")
    parser.add_argument("-T", metavar="sr", type=str, default="", help="type of reads provided by '-R'. 'sr' and 'corr' accepted for short reads and error corrected long reads respectively.")

    # If one uses the same output directory name for both correct and scaffold, could that be ok without having to worry
    # about re-writing files? e.g.
    # ragoo_correct.py -o test_out
    # ragoo_scaffold.py -o test_out
    # Would that be a problem?

    args = parser.parse_args()
    reference_file = os.path.abspath(args.reference)
    query_file = os.path.abspath(args.query)
    output_path = args.o.replace("/", "").replace(".", "")
    num_threads = args.t
    min_ulen = args.f
    min_q = args.q
    merge_dist = args.d
    min_break_dist = args.m
    min_break_end_dist = args.b
    val_window_size = args.v
    overwrite_files = args.w

    gff_file = args.gff
    if gff_file:
        gff_file = os.path.abspath(gff_file)

    skip_file = args.j
    if skip_file:
        skip_file = os.path.abspath(args.j)

    exclude_file = args.e
    if exclude_file:
        exclude_file = os.path.abspath(args.e)

    # Get aligner arguments
    genome_aligner_path = args.genome_aligner
    genome_aligner = genome_aligner_path.split("/")[-1]
    if genome_aligner.split("/")[-1] not in {'minimap2', 'nucmer'}:
        raise ValueError("Must specify either 'minimap2' or 'nucmer' (PATHs allowed) with '--aligner'.")
    mm2_params = args.mm2_params
    nucmer_params = args.nucmer_params

    # Add the number of mm2 threads if the mm2 params haven't been overridden.
    if mm2_params == "-k19 -w19":
        mm2_params += " -t" + str(num_threads)

    # Make sure quality filtering is disabled for nucmer alignments
    if genome_aligner == "nucmer":
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

    # read alignment parameters
    corr_reads = args.R
    corr_reads_fofn = args.F
    corr_reads_tech = args.T
    read_aligner_path = args.read_aligner
    read_aligner = read_aligner_path.split("/")[-1]
    if read_aligner != "minimap2":
        raise ValueError("Only minimap2 can be used for read alignments. got: %s" % read_aligner)

    # If the genome aligner is minimap2, we can just use that path for read alignment
    if genome_aligner == 'minimap2':
        read_aligner_path = genome_aligner_path

    # Make sure that if -R or -F, -T has been specified.
    if corr_reads or corr_reads_fofn:
        if not corr_reads_tech:
            raise ValueError("'-T' must be provided when using -R or -F.")

    # Make a list of read sequences.
    read_files = []
    if corr_reads_fofn:
        with open(corr_reads_fofn, "r") as f:
            for line in f:
                read_files.append(line.rstrip())
    elif corr_reads:
        read_files.append(corr_reads)

    # Get the skip and exclude sets.
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

    # Get the current working directory and output path.
    cwd = os.getcwd()
    output_path = cwd + "/" + output_path + "/"
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # Align the query to the reference.
    log("Aligning the query genome to the reference")
    if genome_aligner == "minimap2":
        al = Minimap2Aligner(reference_file, query_file, genome_aligner_path, mm2_params, output_path + "c_query_against_ref", in_overwrite=overwrite_files)
    else:
        al = NucmerAligner(reference_file, query_file, genome_aligner_path, nucmer_params, output_path + "c_query_against_ref", in_overwrite=overwrite_files)
    al.run_aligner()

    # If alignments are from Nucmer, convert from delta to paf.
    if genome_aligner == "nucmer":
        cmd = ["delta2paf.py", output_path + "c_query_against_ref.delta", ">", output_path + "c_query_against_ref.paf"]
        run(" ".join(cmd))

    # Read and organize the alignments.
    log('Reading whole genome alignments')
    # Read all of the alignments into a temporary structure to save time on ContigAlignment instantiation.
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

    # Put the tmp alignments into ContigAlignment objects.
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

    # Filter and merge the alignments.
    log("Filtering and merging alignments")
    for i in ctg_alns:
        ctg_alns[i] = ctg_alns[i].filter_mapq(min_q)
        if ctg_alns[i] is not None:
            ctg_alns[i] = ctg_alns[i].unique_anchor_filter(min_ulen)
            if ctg_alns[i] is not None:
                ctg_alns[i] = ctg_alns[i].merge_alns(merge_dist=merge_dist)

    # Get the putative breakpoints for each query sequence, if any.
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

    # If desired, validate the putative breakpoints by observing read coverage.
    if read_files:
        log("Validating putative query breakpoints via read alignment.")
        if not os.path.isfile(output_path + "c_reads_against_query.s.bam"):
            log("Aligning reads to query sequences.")
            if corr_reads_tech == "sr":
                al = Minimap2SAMAligner(query_file, " ".join(read_files), read_aligner_path, "-ax sr -t " + str(num_threads),
                                        output_path + "c_reads_against_query", in_overwrite=overwrite_files)
            elif corr_reads_tech == "corr":
                al = Minimap2SAMAligner(query_file, " ".join(read_files), read_aligner_path, "-ax asm5 -t " + str(num_threads),
                                        output_path + "c_reads_against_query", in_overwrite=overwrite_files)
            else:
                raise ValueError("'-T' must be either 'sr' or 'corr'.")
            al.run_aligner()

        # Compress, sort and index the alignments.
        run_samtools(output_path, num_threads, overwrite_files)

        # Validate the breakpoints
        log("Validating putative query breakpoints")
        ctg_breaks = validate_breaks(ctg_breaks, output_path, num_threads, overwrite_files, window_size=val_window_size, clean_dist=min_break_dist)

    # Check if we need to avoid gff intervals
    if gff_file:
        log("Avoiding breaks within GFF intervals")
        it = make_gff_interval_tree(gff_file)
        non_gff_breaks = dict()
        for ctg in ctg_breaks:
            new_breaks = []
            for i in ctg_breaks[ctg]:
                if it[ctg][i]:
                    log("Avoiding breaking %s at %d. This point intersects a feature in the gff file.)" % (ctg, i))
                else:
                    new_breaks.append(i)
            if new_breaks:
                non_gff_breaks[ctg] = new_breaks
        ctg_breaks = non_gff_breaks

    # Write the summary of query sequence breaks in BED format
    bed_file = output_path + "correction.placement.bed"
    write_breaks(bed_file, query_file, ctg_breaks, overwrite_files, output_path)

    # Write the scaffolds.
    log("Writing broken contigs")
    qf_name = query_file.split("/")[-1]
    qf_pref = qf_name[:qf_name.rfind(".")]
    cmd = [
        "break_query.py",
        bed_file,
        query_file,
        output_path + qf_pref + ".break.fasta"
    ]
    run(" ".join(cmd))


if __name__ == "__main__":
    main()
