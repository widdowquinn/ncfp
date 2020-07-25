# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) University of Strathclyde 2019-2020
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G1 1XQ
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2017-2019 The James Hutton Institute
# Copyright (c) 2019-2020 University of Strathclyde
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
"""Test command-line parsing for ncfp program"""

import logging
import unittest

from pathlib import Path

import pytest

from ncbi_cds_from_protein.scripts import ncfp
from ncbi_cds_from_protein import NCFPException


class TestCLIParsing(unittest.TestCase):
    """Class defining tests of ncfp CLI parsing."""

    def setUp(self):
        """Set attributes for tests."""
        self.indir = Path("tests/test_input/sequences")
        self.outdir = Path("tests/test_output/parsertests")
        self.email = "ncfptest@dev.null"
        # Lists of command-line arguments for each test
        self.argsdict = {
            "validate_and_log": [
                str(self.indir / "input_uniprot_stockholm_small.fasta"),
                str(self.outdir),
                self.email,
                "--stockholm",
                "--disabletqdm",
                "-l",
                str(self.outdir / "download_logfile.log"),
            ],
            "bad_infile": [str(self.indir / "notexist.fasta"), str(self.outdir), self.email, "--disabletqdm",],
            "local_cache": [
                str(self.indir / "human.fasta"),
                str(self.outdir),
                self.email,
                "--disabletqdm",
                "-d",
                str(self.outdir / "cache"),
                "-c",
                "humancache",
            ],
        }

        # Null logger for testing
        self.logger = logging.getLogger("TestCLIParsing logger")
        self.logger.addHandler(logging.NullHandler())

    def test_download_and_log(self):
        """ncfp downloads coding sequences and logs output from CLI."""
        ncfp.run_main(self.argsdict["validate_and_log"])

    def test_bad_infile(self):
        """ncfp stops if CLI input file does not exist."""
        with pytest.raises(NCFPException):
            ncfp.run_main(self.argsdict["bad_infile"])

    def test_create_and_keep_cache(self):
        """ncfp creates named cache from CLI and keeps it when rerunning."""
        self.logger.info("Creating local cache")
        ncfp.run_main(self.argsdict["local_cache"])

        self.logger.info("Reusing local cache")
        ncfp.run_main(self.argsdict["local_cache"] + ["--keepcache"])
