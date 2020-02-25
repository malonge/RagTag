import time
import subprocess

""" A collection of various helper functions"""

complements = str.maketrans("ACGTNURYSWKMBVDHacgtnuryswkmbvdh", "TGCANAYRSWMKVBHDtgcanayrswmkvbhd")


def reverse_complement(seq):
    """
    Reverse complement a nucleotide sequence.
    :param seq: Sequence to be reverse complemented
    :return: A reverse complemented sequence
    """
    return seq.translate(complements)[::-1]


def run(cmnd):
    """ Run command and report status. """
    log('Running : %s' % cmnd)
    if subprocess.call(cmnd, shell=True, executable='/bin/bash') != 0:
        raise RuntimeError('Failed : %s ' % cmnd)
    log('Finished running : %s' % cmnd)


def log(message):
    """ Log messages to standard output. """
    print(time.ctime() + ' --- ' + message, flush=True)

