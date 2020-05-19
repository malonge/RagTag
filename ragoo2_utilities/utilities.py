import time
import subprocess
import operator

""" A collection of various helper functions"""

complements = str.maketrans("ACGTNURYSWKMBVDHacgtnuryswkmbvdh", "TGCANAYRSWMKVBHDtgcanayrswmkvbhd")


def reverse_complement(seq):
    """
    Reverse complement a nucleotide sequence.
    :param seq: Sequence to be reverse complemented
    :return: A reverse complemented sequence
    """
    return seq.translate(complements)[::-1]


def run(cmd):
    """ Run command and report status. """
    if not isinstance(cmd, list):
        raise TypeError("'run' expects a list")

    log('Running : %s' % " ".join(cmd))
    if subprocess.call(cmd) != 0:
        raise RuntimeError('Failed : %s' % " ".join(cmd))
    log('Finished running : %s' % " ".join(cmd))


def run_oe(cmd, out, err):
    """ Run command and redirect stdout/stderr. """

    if not isinstance(out, str) or not isinstance(err, str):
        raise TypeError("out/err should be file names (strings)")

    f_out = open(out, "w")
    f_err = open(err, "w")

    log('Running : %s > %s 2> %s' % (" ".join(cmd), out, err))
    if subprocess.call(cmd, stdout=f_out, stderr=f_err) != 0:
        raise RuntimeError('Failed : %s > %s 2> %s' % (" ".join(cmd), out, err))

    log('Finished running : %s > %s 2> %s' % (" ".join(cmd), out, err))

    f_out.close()
    f_err.close()


def run_o(cmd, out):
    """ Run command and redirect stdout but not stderr. """

    if not isinstance(out, str):
        raise TypeError("out should be a file name (string)")

    f_out = open(out, "w")

    log('Running : %s > %s' % (" ".join(cmd), out))
    if subprocess.call(cmd, stdout=f_out) != 0:
        raise RuntimeError('Failed : %s > %s' % (" ".join(cmd), out))

    log('Finished running : %s > %s' % (" ".join(cmd), out))

    f_out.close()


def log(message):
    """ Log messages to standard output. """
    print(time.ctime() + ' --- ' + message, flush=True)


def binary_search(query, numbers, left, right):
    """
    The contents of this method are either influenced by or directly copied from "Assemblytics_uniq_anchor.py"
    written by Maria Nattestad. The original script can be found here:

    https://github.com/MariaNattestad/Assemblytics

    And the publication associated with Maria's work is here:

    Nattestad, Maria, and Michael C. Schatz. "Assemblytics: a
    web analytics tool for the detection of variants from an
    assembly." Bioinformatics 32.19 (2016): 3021-3023.
    """
    #  Returns index of the matching element or the first element to the right

    if left >= right:
        return right
    mid = (right + left) // 2

    if query == numbers[mid]:
        return mid
    elif query < numbers[mid]:
        return binary_search(query, numbers, left, mid)
    else:
        return binary_search(query, numbers, mid + 1, right)


def summarize_planesweep(lines, unique_length_required, keep_small_uniques=False):
    """
    The contents of this method are either influenced by or directly copied from "Assemblytics_uniq_anchor.py"
    written by Maria Nattestad. The original script can be found here:

    https://github.com/MariaNattestad/Assemblytics

    And the publication associated with Maria's work is here:

    Nattestad, Maria, and Michael C. Schatz. "Assemblytics: a
    web analytics tool for the detection of variants from an
    assembly." Bioinformatics 32.19 (2016): 3021-3023.
    """
    alignments_to_keep = []

    # If no alignments:
    if len(lines) == 0:
        return []

    # If only one alignment:
    if len(lines) == 1:
        if keep_small_uniques == True or abs(lines[0][1] - lines[0][0]) >= unique_length_required:
            return [0]
        else:
            return []

    starts_and_stops = []
    for query_min, query_max in lines:
        starts_and_stops.append((query_min, "start"))
        starts_and_stops.append((query_max, "stop"))

    sorted_starts_and_stops = sorted(starts_and_stops, key=operator.itemgetter(0))

    current_coverage = 0
    last_position = -1
    sorted_unique_intervals_left = []
    sorted_unique_intervals_right = []
    for pos, change in sorted_starts_and_stops:

        if current_coverage == 1:
            sorted_unique_intervals_left.append(last_position)
            sorted_unique_intervals_right.append(pos)

        if change == "start":
            current_coverage += 1
        else:
            current_coverage -= 1
        last_position = pos

    linecounter = 0
    for query_min, query_max in lines:

        i = binary_search(query_min, sorted_unique_intervals_left, 0, len(sorted_unique_intervals_left))

        exact_match = False
        if sorted_unique_intervals_left[i] == query_min and sorted_unique_intervals_right[i] == query_max:
            exact_match = True
        sum_uniq = 0
        while i < len(sorted_unique_intervals_left) and sorted_unique_intervals_left[i] >= query_min and \
                        sorted_unique_intervals_right[i] <= query_max:
            sum_uniq += sorted_unique_intervals_right[i] - sorted_unique_intervals_left[i]
            i += 1

        if sum_uniq >= unique_length_required:
            alignments_to_keep.append(linecounter)
        elif keep_small_uniques == True and exact_match == True:
            alignments_to_keep.append(linecounter)

        linecounter += 1

    return alignments_to_keep
