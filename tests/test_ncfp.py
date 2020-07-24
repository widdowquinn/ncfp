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
"""Test commandline entry for the ncfp program."""

import logging
import os
import time
import unittest

from argparse import Namespace
from pathlib import Path

from ncbi_cds_from_protein.scripts import ncfp


class TestBasicUse(unittest.TestCase):
    """Class defining tests of basic ncfp usage."""

    def setUp(self):
        """Set attributes for tests."""
        # Paths to test file directories
        self.datadir = Path("tests/test_input/sequences")
        self.outdir = Path("tests/test_output")
        self.targetdir = Path("tests/test_targets")

        # Null logger instance
        self.logger = logging.getLogger("TestBasicUse logger")
        self.logger.addHandler(logging.NullHandler())

        # Default command-line namespace
        self.base_namespace = Namespace(
            infname=self.datadir / "input_ncbi.fasta",
            outdirname=self.outdir / "basic_ncbi",
            email="ncfp@dev.null",
            infmt="ncbi",
            stockholm=False,
            cachedir=Path(".ncfp_cache"),
            cachestem=time.strftime("%Y-%m-%d-%H-%m-%S"),
            batchsize=100,
            retries=10,
            limit=None,
            filestem="ncfp",
            keepcache=False,
            skippedfname="skipped.fasta",
            logfile=None,
            verbose=False,
            disabletqdm=True,
            debug=False,
        )

    def test_basic_ncbi(self):
        """ncfp collects correct coding sequences for basic NCBI input."""
        # Modify default arguments
        outdirname = "basic_ncbi"  # common dirname for output,target
        namespace = self.base_namespace
        namespace.infname = self.datadir / "input_ncbi.fasta"
        namespace.outdirname = self.outdir / outdirname

        # Run ersatz command-line
        ncfp.run_main(namespace)

        # Compare output
        self.check_files(outdirname, ("ncfp_aa.fasta", "ncfp_nt.fasta"))

    def test_basic_uniprot(self):
        """ncfp collects correct coding sequences for basic UniProt input."""
        # Modify default arguments
        outdirname = "basic_uniprot"  # common dirname for output,target
        namespace = self.base_namespace
        namespace.infname = self.datadir / "input_uniprot.fasta"
        namespace.outdirname = self.outdir / outdirname
        namespace.infmt = "uniprot"

        # Run ersatz command-line
        ncfp.run_main(namespace)

        # Compare output
        self.check_files(outdirname, ("ncfp_aa.fasta", "ncfp_nt.fasta", "skipped.fasta"))

    def test_basic_stockholm(self):
        """ncfp collects correct coding sequences for basic UniProt/Stockholm input."""
        # Modify default arguments
        outdirname = "small_stockholm"  # common dirname for output,target
        namespace = self.base_namespace
        namespace.infname = self.datadir / "input_uniprot_stockholm_small.fasta"
        namespace.outdirname = self.outdir / outdirname
        namespace.infmt = "uniprot"
        namespace.stockholm = True

        # Run ersatz command-line
        ncfp.run_main(namespace)

        # Compare output
        self.check_files(outdirname, ("ncfp_aa.fasta", "ncfp_nt.fasta"))

    def check_files(self, outdirname, fnames):
        """Test whether output and target files are identical"""
        for fname in fnames:
            with (self.outdir / outdirname / fname).open() as fh1:
                with (self.targetdir / outdirname / fname).open() as fh2:
                    assert fh1.read() == fh2.read()
