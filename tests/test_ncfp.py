#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_cmdlines.py

Test commandline entry for the ncfp program.

This test suite is intended to be run from the repository root using:

nosetests -v

For each test, command-line options are defined in a Namespace,
and passed as the sole argument to the run_main() function in
ncfp.py.

Paths are defined relative to the repository root, where the tests
are expected to be run.

(c) The James Hutton Institute 2017
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD6 9LH,
Scotland,
UK

The MIT License

Copyright (c) 2017 The James Hutton Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import logging
import os
import time
import unittest

from argparse import Namespace

from nose.tools import (assert_equal, raises)

from ncbi_cds_from_protein.scripts import ncfp


class TestBasicUse(unittest.TestCase):
    """Class defining tests of basic ncfp command lines."""

    def setUp(self):
        """Set attributes for tests."""
        # Paths to test file directories
        self.datadir = os.path.join('tests', 'test_input', 'sequences')
        self.outdir = os.path.join('tests', 'test_output')
        self.targetdir = os.path.join('tests', 'test_targets')

        # Null logger instance
        self.logger = logging.getLogger('TestBasicUse logger')
        self.logger.addHandler(logging.NullHandler())

        # Default command-line namespace
        self.base_namespace = Namespace(
            infname=os.path.join(self.datadir, 'input_ncbi.fasta'),
            outdirname=os.path.join(self.outdir, 'basic_ncbi'),
            email="ncfp@dev.null",
            uniprot=False,
            stockholm=False,
            cachedir='.ncfp_cache',
            cachestem=time.strftime("%Y-%m-%d-%H-%m-%S"),
            batchsize=100,
            retries=10,
            limit=None,
            filestem="ncfp",
            keepcache=False,
            skippedfname="skipped.fasta",
            logfile=None,
            verbose=False,
            disabletqdm=True
        )

    def test_basic_ncbi(self):
        """ncfp collects correct coding sequences for basic NCBI input."""
        # Modify default arguments
        outdirname = 'basic_ncbi'        # common dirname for output,target
        namespace = self.base_namespace
        namespace.infname = os.path.join(
            self.datadir, 'input_ncbi.fasta')
        namespace.outdirname = os.path.join(self.outdir, outdirname)

        # Run ersatz command-line
        ncfp.run_main(namespace)

        # Compare output
        for fname in ('ncfp_aa.fasta', 'ncfp_nt.fasta'):
            with open(os.path.join(namespace.outdirname, fname)) as fh1:
                with open(os.path.join(self.targetdir, outdirname, fname)) as fh2:
                    assert_equal(fh1.read(), fh2.read())

    def test_basic_uniprot(self):
        """ncfp collects correct coding sequences for basic UniProt input."""
        # Modify default arguments
        outdirname = 'basic_uniprot'        # common dirname for output,target
        namespace = self.base_namespace
        namespace.infname = os.path.join(self.datadir, 'input_uniprot.fasta')
        namespace.outdirname = os.path.join(self.outdir, outdirname)
        namespace.uniprot = True

        # Run ersatz command-line
        ncfp.run_main(namespace, self.logger)

        # Compare output
        for fname in ('ncfp_aa.fasta', 'ncfp_nt.fasta', 'skipped.fasta'):
            with open(os.path.join(namespace.outdirname, fname)) as fh1:
                with open(os.path.join(self.targetdir, outdirname, fname)) as fh2:
                    assert_equal(fh1.read(), fh2.read())

    def test_basic_stockholm(self):
        """ncfp collects correct coding sequences for basic UniProt/Stockholm input."""
        # Modify default arguments
        outdirname = 'small_stockholm'        # common dirname for output,target
        namespace = self.base_namespace
        namespace.infname = os.path.join(
            self.datadir, 'input_uniprot_stockholm_small.fasta')
        namespace.outdirname = os.path.join(self.outdir, outdirname)
        namespace.uniprot = True
        namespace.stockholm = True

        # Run ersatz command-line
        ncfp.run_main(namespace)

        # Compare output
        for fname in ('ncfp_aa.fasta', 'ncfp_nt.fasta'):
            with open(os.path.join(namespace.outdirname, fname)) as fh1:
                with open(os.path.join(self.targetdir, outdirname, fname)) as fh2:
                    assert_equal(fh1.read(), fh2.read())
