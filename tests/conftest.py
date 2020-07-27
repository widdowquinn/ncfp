# -*- coding: utf-8 -*-
# (c) The University of Strathclude 2019-2020
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute of Pharmaceutical and Biomedical Sciences
# The University of Strathclyde
# 161 Cathedral Street
# Glasgow
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# (c) The University of Strathclude 2019-2020
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
"""Pytest configuration file."""

from pathlib import Path

import pytest


# Path to tests, contains tests and data subdirectories
# This conftest.py file should be found in the top directory of the tests
# module. The fixture data should be in a subdirectory named fixtures
TESTSPATH = Path(__file__).parents[0]
FIXTUREPATH = TESTSPATH / "fixtures"
TARGETPATH = FIXTUREPATH / "targets"


@pytest.fixture
def dir_sequences():
    """Path to directory containing sequence file fixtures."""
    yield FIXTUREPATH / "sequences"


@pytest.fixture
def email_address():
    """Package-specific email address, passed to NCBI."""
    yield "ncfptest@dev.null"


@pytest.fixture
def path_human(dir_sequences):
    """Path to FASTA file of human sequences."""
    yield dir_sequences / "human.fasta"


@pytest.fixture
def path_ncbi(dir_sequences):
    """Path to NCBI sequences."""
    yield dir_sequences / "input_ncbi.fasta"


@pytest.fixture
def path_ncbi_targets():
    """Path to NCBI target outputs."""
    yield TARGETPATH / "basic_ncbi"


@pytest.fixture
def path_notexist(dir_sequences):
    """Path to nonexistent FASTA file."""
    yield dir_sequences / "notexist.fasta"


@pytest.fixture
def path_stockholm(dir_sequences):
    """Path to Stockholm sequences."""
    yield dir_sequences / "input_uniprot_stockholm.fasta"


@pytest.fixture
def path_stockholm_targets():
    """Path to Stockholm target outputs."""
    yield TARGETPATH / "basic_stockholm"


@pytest.fixture
def path_uniprot(dir_sequences):
    """Path to UniProt sequences."""
    yield dir_sequences / "input_uniprot.fasta"


@pytest.fixture
def path_uniprot_stockholm_small(dir_sequences):
    """Path to small FASTA file of UniProt sequences with Stockholm format."""
    yield dir_sequences / "input_uniprot_stockholm_small.fasta"


@pytest.fixture
def path_uniprot_targets():
    """Path to UniProt target outputs."""
    yield TARGETPATH / "basic_uniprot"


@pytest.fixture
def path_uniprot_stockholm_small_targets(dir_sequences):
    """Path to targets for small FASTA file of UniProt sequences with Stockholm format."""
    yield TARGETPATH / "small_stockholm"
