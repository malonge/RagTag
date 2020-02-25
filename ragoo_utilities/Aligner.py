import abc
import os

from ragoo_utilities.utilities import run, log


class Aligner:
    """
    A description of the base class
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, in_ref_file, in_query_file, in_aligner, in_params, in_out_file, in_overwrite=False):
        self.r_file = in_ref_file
        self.q_file = in_query_file
        self.aligner = in_aligner
        self.params_string = in_params
        self.outfile_prefix = in_out_file
        self.overwrite = in_overwrite
        self.out_file = None

    @abc.abstractmethod
    def params_are_valid(self):
        """
        Check to ensure that not illegal (for RaGOO) parameters are being passed to the aligner
        :return: True if no illegal parameters
        """
        pass

    @abc.abstractmethod
    def compile_command(self):
        """
        Compile the space delimited command to be executed.
        :return: The command as a string
        """
        pass

    @abc.abstractclassmethod
    def output_exists(self):
        """
        Check if the output file already exists
        :return: True if output exists, False otherwise
        """
        pass


"""
class NucmerAligner(Aligner):
    
    Description of this inheriting class.

    def params_are_valid(self):
        if "-p" in self.params_string:
            raise ValueError("Please don't specify '-p' when using nucmer. RaGOO names its own alignment files.")

    def compile_command(self):
        return " ".join([
            self.aligner,
            self.params_string,
            '-p ' + self.outfile_prefix,
            self.r_file,
            self.q_file,
            ">",
            "/dev/null",
            "2>",
            self.outfile_prefix + ".err"
        ])

"""


class Minimap2Aligner(Aligner):
    """
    Description of this inheriting class.
    """

    def __init__(self, in_ref_file, in_query_file, in_aligner, in_params, in_out_file, in_overwrite=False):
        super().__init__(in_ref_file, in_query_file, in_aligner, in_params, in_out_file, in_overwrite)
        self.out_file = self.outfile_prefix + ".paf"
        self.out_log = self.out_file + ".log"

    def params_are_valid(self):
        """
        Do a basic check to make sure the Minimap2 parameters are valid.
        I won't check that every parameter is valid, but will check anything that can
        cause a problem for RaGOO later on.
        :return: True if the parameters are valid. Raises appropriate errors otherwise
        """
        all_flags = "".join([i for i in self.params_string.split(" ") if i.startswith("-")])
        if "a" in all_flags:
            raise ValueError("Please don't output Minimap2 alignments in SAM format (-a). RaGOO needs PAF format.")

        if "c" in all_flags:
            log("WARNING: getting base alignments (-c) will slow down Minimap2 alignment.")

        return True

    def compile_command(self):
        """
        Compile the space delimited command to be executed.
        :return: The command as a string
        """
        return " ".join([
            self.aligner,
            self.params_string,
            self.r_file,
            self.q_file,
            ">",
            self.out_file,
            "2>",
            self.out_log
        ])

    def output_exists(self):
        return os.path.isfile(self.out_file)

    def run_aligner(self):
        """ Run the aligner. """
        if self.params_are_valid():
            if not self.output_exists():
                run(self.compile_command())
            else:
                if self.overwrite:
                    run(self.compile_command())
                else:
                    log("retaining pre-existing file: " + self.out_file)