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
import argparse

from ragtag_utilities.utilities import log
from ragtag_utilities.AGPFile import AGPFile


def main():
    parser = argparse.ArgumentParser(description="Check AGP v2.1 files for validity.", usage="ragtag.py agpcheck <asm1.agp> [<asm2.agp> ... <asmN.agp>]")
    parser.add_argument("agp", metavar="<asm1.agp> [<asm2.agp> ... <asmN.agp>]", nargs='+', default=[], type=str, help="AGP v2.1 files")

    DISCLAIMER = """
    DISCLAIMER:
    This utility performs most (but not all) checks necessary to validate an
    AGP v2.1 file: https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
    
    Please additionally use the NCBI AGP validator for robust
    validation: https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Validation/
    """

    args = parser.parse_args()

    print(DISCLAIMER)
    agp_file_list = [os.path.abspath(i) for i in args.agp]
    for agp_file in agp_file_list:
        print()
        log("Checking {} ...".format(agp_file))
        agp = AGPFile(agp_file, mode="r")
        for _ in agp.iterate_lines():
            pass
        log("Check for {} is complete with no errors.".format(agp_file))


if __name__ == "__main__":
    main()
