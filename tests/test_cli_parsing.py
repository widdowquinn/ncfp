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


import pytest

from ncbi_cds_from_protein.scripts import ncfp
from ncbi_cds_from_protein import NCFPException


@pytest.fixture
def args_bad_infile(email_address, path_notexist, tmp_path):
    """Cmd-line arguments for passing a nonexistent input file."""
    return [path_notexist, tmp_path, email_address, "--disabletqdm"]


@pytest.fixture
def args_create_reuse_cache(email_address, path_human, tmp_path):
    """Cmd-line arguments for creating and reusing a local cache."""
    return [
        path_human,
        tmp_path,
        email_address,
        "--disabletqdm",
        "-d",
        tmp_path / "cache",
        "-c",
        "humancache",
    ]


@pytest.fixture
def args_validate_and_log(
    email_address, path_uniprot_stockholm_small, tmp_path,
):
    """Cmd-line arguments for downloading sequences and logging output.

    Validates that ncfp runs without error, does not check output
    """
    return [
        path_uniprot_stockholm_small,
        tmp_path,
        email_address,
        "--stockholm",
        "--disabletqdm",
        "-l",
        tmp_path / "validate_and_log.log",
    ]


def test_bad_infile(args_bad_infile):
    """ncfp stops if CLI input file does not exist.

    Validates that ncfp throws correct error, does not check output
    """
    with pytest.raises(NCFPException):
        ncfp.run_main(args_bad_infile)


def test_create_and_keep_cache(args_create_reuse_cache):
    """ncfp creates named cache from CLI and keeps it when rerunning.

    This calls ncfp twice and expects both calls run without error. Output
    is not checked.
    """
    # Create local cache
    ncfp.run_main(args_create_reuse_cache)
    # Reuse created cache
    ncfp.run_main(args_create_reuse_cache + ["--keepcache"])


def test_download_and_log(args_validate_and_log):
    """ncfp downloads coding sequences and logs output from CLI.

    Validates that ncfp runs without error, does not check output
    """
    ncfp.run_main(args_validate_and_log)
