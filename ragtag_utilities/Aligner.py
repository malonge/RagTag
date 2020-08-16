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

import abc
import os
import shutil

from ragtag_utilities.utilities import run_oe, run_e, log


class Aligner:
    """
    This class is an Abstract Base Class and should be inherited by subclasses with
    the goal of running pairwise sequence aligners on the command line. Subclasses should be
    aligner-specific, specifying the parameters/output and performing all validity
    checks.

    The base class specifies generic information that is sufficient to run most pairwise aligners.
    The generic information is as follows:

        1. A file handle pointing to a reference sequence
        2. A file handle pointing to a set of query sequences
        3. Command line parameters
        4. A file handle (or sometimes prefix) for the output alignments
        5. A file handle for logging (i.e. stderr from the aligner)

    Inheriting classes should perform the following tasks:
        1. Check that all the file handles/executables are valid
        2. Parse the command line arguments
        3. Ensure a valid command can be constructed
        4. Check if the output alignments already exist
        5. Execute the command (run the aligner via subprocess)

    A small number of these tasks are generic (or mostly generic) and can be specified in this base class.
    For example, the base class can check that the aligner executable (e.g. 'minimap2') exists. However,
    the majority of these tasks are specific to the aligner, and should be specified in inheriting classes
    by overriding base class methods or by creating new methods.
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, in_ref_file, in_query_files, in_aligner, in_params, in_out_file, in_overwrite=False):
        self.r_file = in_ref_file
        self.q_files = in_query_files
        self.aligner = in_aligner
        self.params_string = in_params
        self.params = self._split_params(in_params)
        self.outfile_prefix = in_out_file
        self.overwrite = in_overwrite

        # Aligner specific parameters
        self.out_file = None
        self.aligner_exec = None
        self.out_log = None

        self._update_attrs()

    @abc.abstractmethod
    def _update_attrs(self):
        """ Update class attributes for a specific aligner. """
        pass

    @staticmethod
    def _split_params(a):
        """ Split aligner parameters into a suitable format for 'subprocess.call'. """
        return a.split(" ")

    @abc.abstractmethod
    def params_are_valid(self):
        """
        Check to ensure that not illegal (for RagTag) parameters are being passed to the aligner
        :return: True if no illegal parameters
        """
        pass

    @abc.abstractmethod
    def compile_command(self):
        """
        Compile the space delimited command to be executed.
        :return: The command as a list of strings
        """
        pass

    def exec_is_valid(self):
        """
        Check if the aligner executable is valid.
        :return: True if the executable is valid. Raise appropriate error otherwise.
        """
        ex = self.aligner.split("/")[-1]
        if ex != self.aligner_exec:
            raise ValueError("The aligner executable must be `" + self.aligner_exec + "`. Got `" + ex + "'")

        if not shutil.which(self.aligner):
            raise ValueError("The provided aligner executable is not valid: " + self.aligner)

        return True

    def output_exists(self):
        """ Check if the output file already exists. """
        return os.path.isfile(self.out_file)

    def run_aligner(self):
        """ Run the aligner. """
        if all([self.params_are_valid(), self.exec_is_valid()]):
            if not self.output_exists():
                run_oe(self.compile_command(), self.out_file, self.out_log)
            else:
                if self.overwrite:
                    log("overwriting pre-existing file: " + self.out_file)
                    run_oe(self.compile_command(), self.out_file, self.out_log)
                else:
                    log("Retaining pre-existing file: " + self.out_file)


class NucmerAligner(Aligner):
    """ The 'Aligner' subclass specific to whole genome alignment with Nucmer """

    def _update_attrs(self):
        """ Update class attributes for a specific aligner. """
        self.aligner_exec = "nucmer"
        self.out_file = self.outfile_prefix + ".delta"
        self.out_log = self.out_file + ".log"

    def params_are_valid(self):
        """
        Do a basic check to make sure the nucmer parameters are valid.
        I won't check that every parameter is valid, but will check anything that can
        cause a problem for RagTag later on.
        :return: True if the parameters are valid. Raises appropriate errors otherwise
        """
        if "-p" in self.params_string:
            raise ValueError("Please don't specify '-p' when using nucmer. RagTag names its own alignment files.")

        return True

    def compile_command(self):
        return [
            self.aligner,
            *self.params,
            '-p',
            self.outfile_prefix,
            self.r_file,
            *self.q_files
        ]

    def run_aligner(self):
        """ Run the aligner. """
        if all([self.params_are_valid(), self.exec_is_valid()]):
            if not self.output_exists():
                run_e(self.compile_command(), self.out_log)
            else:
                if self.overwrite:
                    log("overwriting pre-existing file: " + self.out_file)
                    run_e(self.compile_command(), self.out_log)
                else:
                    log("Retaining pre-existing file: " + self.out_file)


class Minimap2Aligner(Aligner):
    """ The 'Aligner' subclass specific to Minimap2 alignments in PAF format. """

    def _update_attrs(self):
        """ Update class attributes for a specific aligner. """
        self.aligner_exec = "minimap2"
        self.out_file = self.outfile_prefix + ".paf"
        self.out_log = self.out_file + ".log"

    def params_are_valid(self):
        """
        Do a basic check to make sure the minimap2 parameters are valid.
        I won't check that every parameter is valid, but will check anything that can
        cause a problem for RagTag later on.
        :return: True if the parameters are valid. Raises appropriate errors otherwise
        """
        all_flags = "".join([i for i in self.params_string.split(" ") if i.startswith("-")])
        if "a" in all_flags:
            raise ValueError("Alignments must not be in SAM format (-a).")

        if "c" in all_flags:
            log("WARNING: computing base-alignments (-c) will slow down Minimap2 alignment.")

        return True

    def compile_command(self):
        """
        Compile the space delimited command to be executed.
        :return: The command as a string
        """
        return [
            self.aligner,
            *self.params,
            self.r_file,
            *self.q_files
        ]


class Minimap2SAMAligner(Aligner):
    """ The 'Aligner' subclass specific to Minimap2 alignments in SAM format. """

    def _update_attrs(self):
        """ Update class attributes for a specific aligner. """
        self.aligner_exec = "minimap2"
        self.out_file = self.outfile_prefix + ".sam"
        self.out_log = self.out_file + ".log"

    def params_are_valid(self):
        """
        Do a basic check to make sure the minimap2 parameters are valid.
        I won't check that every parameter is valid, but will check anything that can
        cause a problem for RagTag later on.
        :return: True if the parameters are valid. Raises appropriate errors otherwise
        """
        all_flags = "".join([i for i in self.params_string.split(" ") if i.startswith("-")])
        if "a" not in all_flags:
            raise ValueError("Alignments must be in SAM format (-a).")

        return True

    def compile_command(self):
        """
        Compile the space delimited command to be executed.
        :return: The command as a string
        """
        return [
            self.aligner,
            *self.params,
            self.r_file,
            *self.q_files
        ]
