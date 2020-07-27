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

import time

from argparse import Namespace

import pytest

from ncbi_cds_from_protein.scripts import ncfp
from utils import check_files, modify_namespace


@pytest.fixture
def namespace_base(email_address, path_ncbi, tmp_path):
    """Cmd-line arguments for passing a nonexistent input file."""
    yield Namespace(
        infname=path_ncbi,
        outdirname=tmp_path,
        email=email_address,
        stockholm=False,
        cachedir=tmp_path / ".ncfp_cache",
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


def test_basic_ncbi(namespace_base, path_ncbi, path_ncbi_targets, tmp_path):
    """ncfp collects correct coding sequences for basic NCBI input."""
    # Modify default arguments
    infile = path_ncbi
    outdir = tmp_path / "basic_ncbi"
    args = modify_namespace(namespace_base, infname=infile, outdirname=outdir)

    # Run ersatz command-line
    ncfp.run_main(args)

    # Compare output
    check_files(outdir, path_ncbi_targets, ("ncfp_aa.fasta", "ncfp_nt.fasta"))


def test_basic_uniprot(namespace_base, path_uniprot, path_uniprot_targets, tmp_path):
    """ncfp collects correct coding sequences for basic UniProt input."""
    # Modify default arguments
    infile = path_uniprot
    outdir = tmp_path / "basic_uniprot"
    args = modify_namespace(namespace_base, infname=infile, outdirname=outdir)

    # Run ersatz command-line
    ncfp.run_main(args)

    # Compare output
    check_files(
        outdir,
        path_uniprot_targets,
        ("ncfp_aa.fasta", "ncfp_nt.fasta", "skipped.fasta"),
    )


def test_basic_stockholm(
    namespace_base, path_stockholm, path_stockholm_targets, tmp_path
):
    """ncfp collects correct coding sequences for basic UniProt/Stockholm input."""
    # Modify default arguments
    infile = path_stockholm
    outdir = tmp_path / "small_stockholm"
    args = modify_namespace(
        namespace_base, infname=infile, outdirname=outdir, stockholm=True
    )

    # Run ersatz command-line
    ncfp.run_main(args)

    # Compare output
    check_files(
        outdir,
        path_stockholm_targets,
        ("ncfp_aa.fasta", "ncfp_nt.fasta", "skipped.fasta"),
    )


def test_small_stockholm(
    namespace_base,
    path_uniprot_stockholm_small,
    path_uniprot_stockholm_small_targets,
    tmp_path,
):
    """ncfp collects correct coding sequences for small UniProt/Stockholm input."""
    # Modify default arguments
    infile = path_uniprot_stockholm_small
    outdir = tmp_path / "small_stockholm"
    args = modify_namespace(
        namespace_base, infname=infile, outdirname=outdir, stockholm=True
    )

    # Run ersatz command-line
    ncfp.run_main(args)

    # Compare output
    check_files(
        outdir, path_uniprot_stockholm_small_targets, ("ncfp_aa.fasta", "ncfp_nt.fasta")
    )
