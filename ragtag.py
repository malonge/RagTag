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
    
    assembly improvement:
      correct         misassembly correction 
      scaffold        synteny scaffolding
      merge           scaffold merging
      
    file utilities:
      agp2fasta       build a FASTA file from an AGP file
      agpcheck        check for valid AGP file format
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

        elif cmd == "correct":
            subcmd = ["ragtag_correct.py"] + sys.argv[2:]
            subprocess.call(subcmd)

        elif cmd == "scaffold":
            subcmd = ["ragtag_scaffold.py"] + sys.argv[2:]
            subprocess.call(subcmd)

        elif cmd == "merge":
            subcmd = ["ragtag_merge.py"] + sys.argv[2:]
            subprocess.call(subcmd)

        elif cmd == "agp2fasta":
            subcmd = ["ragtag_agp2fasta.py"] + sys.argv[2:]
            subprocess.call(subcmd)

        elif cmd == "agpcheck":
            subcmd = ["ragtag_agpcheck.py"] + sys.argv[2:]
            subprocess.call(subcmd)

        elif cmd == "updategff":
            subcmd = ["ragtag_update_gff.py"] + sys.argv[2:]
            subprocess.call(subcmd)

        else:
            print(description)
            print("\n** unrecognized command: %s **" % cmd)


if __name__ == "__main__":
    main()
