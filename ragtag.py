#!/usr/bin/env python

import sys
import subprocess

from ragtag_utilities.utilities import get_ragtag_version


def main():
    VERSION = get_ragtag_version()
    CITATION = """
Alonge, Michael, et al. "RaGOO: fast and accurate reference-guided scaffolding of draft genomes."
Genome biology 20.1 (2019): 1-17.
    """

    description = """
RagTag: Reference-guided scaffolding and misassembly correction.
Version: %s

usage: ragtag.py <command> [options]
    
    commands:
      scaffold        scaffold contigs
      correct         correct contig misassemblies
      updategff       update gff intervals

    options:
      -c, --citation  
      -v, --version""" % VERSION

    arg_len = len(sys.argv)
    if arg_len == 1:
        print(description)

    if arg_len > 1:
        cmd = sys.argv[1]

        if cmd == "-h" or cmd == "--help":
            print(description)

        elif cmd == "-v" or cmd == "--version":
            print(VERSION)

        elif cmd == "-c" or cmd == "--citation":
            print(CITATION)

        elif cmd == "scaffold":
            subcmd = ["ragtag_scaffold.py"] + sys.argv[2:]
            subprocess.call(subcmd)

        elif cmd == "correct":
            subcmd = ["ragtag_correct.py"] + sys.argv[2:]
            subprocess.call(subcmd)

        elif cmd == "updategff":
            subcmd = ["ragtag_update_gff.py"] + sys.argv[2:]
            subprocess.call(subcmd)

        else:
            print(description)
            print("\n** unrecognized command: %s" % cmd)


if __name__ == "__main__":
    main()
