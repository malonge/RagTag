import abc
import os
import shutil

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

        # Aligner specific parameters
        self.out_file = None
        self.aligner_exec = None
        self.out_log = None

        self._update_attrs()

    @abc.abstractmethod
    def _update_attrs(self):
        """ Update class attributes for a specific aligner. """
        pass

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
                run(self.compile_command())
            else:
                if self.overwrite:
                    log("overwriting pre-existing file: " + self.out_file)
                    run(self.compile_command())
                else:
                    log("retaining pre-existing file: " + self.out_file)


class NucmerAligner(Aligner):
    """ Description of this inheriting class. """

    def _update_attrs(self):
        """ Update class attributes for a specific aligner. """
        self.aligner_exec = "nucmer"
        self.out_file = self.outfile_prefix + ".delta"
        self.out_log = self.out_file + ".log"

    def params_are_valid(self):
        """
        Do a basic check to make sure the nucmer parameters are valid.
        I won't check that every parameter is valid, but will check anything that can
        cause a problem for RaGOO later on.
        :return: True if the parameters are valid. Raises appropriate errors otherwise
        """
        if "-p" in self.params_string:
            raise ValueError("Please don't specify '-p' when using nucmer. RaGOO names its own alignment files.")

        return True

    def compile_command(self):
        return " ".join([
            self.aligner,
            self.params_string,
            '-p ' + self.outfile_prefix,
            self.r_file,
            self.q_file,
            "2>",
            self.out_log
        ])


class Minimap2Aligner(Aligner):
    """
    Description of this inheriting class.
    """

    def _update_attrs(self):
        """ Update class attributes for a specific aligner. """
        self.aligner_exec = "minimap2"
        self.out_file = self.outfile_prefix + ".paf"
        self.out_log = self.out_file + ".log"

    def params_are_valid(self):
        """
        Do a basic check to make sure the minimap2 parameters are valid.
        I won't check that every parameter is valid, but will check anything that can
        cause a problem for RaGOO later on.
        :return: True if the parameters are valid. Raises appropriate errors otherwise
        """
        all_flags = "".join([i for i in self.params_string.split(" ") if i.startswith("-")])
        if "a" in all_flags:
            raise ValueError("Please don't output Minimap2 alignments in SAM format (-a). RaGOO needs PAF format.")

        if "c" in all_flags:
            log("WARNING: computing base alignments (-c) will slow down Minimap2 alignment.")

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


class Minimap2SAMAligner(Aligner):
    """
    The same as Minimap2Aligner, but outputs in SAM rather than PAF format.
    """

    def _update_attrs(self):
        """ Update class attributes for a specific aligner. """
        self.aligner_exec = "minimap2"
        self.out_file = self.outfile_prefix + ".sam"
        self.out_log = self.out_file + ".log"

    def params_are_valid(self):
        """
        Do a basic check to make sure the minimap2 parameters are valid.
        I won't check that every parameter is valid, but will check anything that can
        cause a problem for RaGOO later on.
        :return: True if the parameters are valid. Raises appropriate errors otherwise
        """
        all_flags = "".join([i for i in self.params_string.split(" ") if i.startswith("-")])
        if "a" not in all_flags:
            raise ValueError("Minimap2 alignments must be in SAM format (-a). RaGOO needs PAF format.")

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