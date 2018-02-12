#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_cli_parsing.py

Test command-line parsing for ncfp program

This test suite is intended to be run from the repository root using:

nosetests -v

(c) The James Hutton Institute 2018
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
import unittest

from nose.tools import (raises,)

from ncbi_cds_from_protein.scripts import (ncfp,)
from ncbi_cds_from_protein.ncfp_tools import (NCFPException,)


class TestCLIParsing(unittest.TestCase):
    """Class defining tests of ncfp cli parsing."""

    def setUp(self):
        """Set attributes for tests."""
        # Lists of command-line arguments for each test
        self.argsdict = {'validate_and_log': [os.path.join('tests', 'test_input',
                                                           'sequences', 'input_uniprot_stockholm_small.fasta'),
                                              os.path.join('tests', 'test_output',
                                                           'parsertests'),
                                              'ncfp@dev.null',
                                              '--uniprot', '--stockholm',
                                              '--disabletqdm',
                                              '-l', os.path.join('tests', 'test_output', 'parsertests',
                                                                 'download_logfile.log'), ],
                         'bad_infile': [os.path.join('tests', 'test_input',
                                                     'sequences', 'notexist.fasta'),
                                        os.path.join('tests', 'test_output',
                                                     'parsertests'),
                                        'ncfp@dev.null', ],
                         }

        # Null logger for testing
        self.logger = logging.getLogger("TestCLIParsing logger")
        self.logger.addHandler(logging.NullHandler())

    def test_download_and_log(self):
        """ncfp downloads coding sequences and logs output from CLI."""
        ncfp.run_main(self.argsdict['validate_and_log'])

    @raises(NCFPException)
    def test_bad_infile(self):
        """ncfp baulks at invalid input file from CLI."""
        ncfp.run_main(self.argsdict['bad_infile'], logger=self.logger)
