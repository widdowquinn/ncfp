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

import io

import pytest

from unittest.mock import Mock

from Bio import Entrez
from bioservices import UniProt

from ncbi_cds_from_protein.scripts import ncfp
from ncbi_cds_from_protein import NCFPException


@pytest.fixture
def args_bad_infile(email_address, path_notexist, tmp_path):
    """Cmd-line arguments for passing a nonexistent input file."""
    yield [path_notexist, tmp_path, email_address, "--disabletqdm"]


@pytest.fixture
def args_create_reuse_cache(email_address, path_human, tmp_path):
    """Cmd-line arguments for creating and reusing a local cache."""
    yield [
        path_human,
        tmp_path,
        email_address,
        "--disabletqdm",
        "-d",
        tmp_path / "cache",
        "-c",
        "humancache",
        "--debug",
    ]


@pytest.fixture
def args_validate_and_log(
    email_address,
    path_uniprot_stockholm_small,
    tmp_path,
):
    """Cmd-line arguments for downloading sequences and logging output.

    Validates that ncfp runs without error, does not check output
    """
    yield [
        path_uniprot_stockholm_small,
        tmp_path,
        email_address,
        "--stockholm",
        "--disabletqdm",
        "-l",
        tmp_path / "validate_and_log.log",
    ]


@pytest.fixture
def mock_uniprot_download_and_log(monkeypatch):
    """Mock remote service call to UniProt for download_and_log test.

    Returns the result expected from l125 of sequences.py, for each query of
    the path_uniprot_stockholm_small dataset

    u_service.search(match.group(0), columns="database(EMBL)")  # type: ignore
    """

    def mock_search(*args, **kwargs):
        """Mock call to UniProt.search() method.

        This output specific to the download_and_log() test
        """
        return "Cross-reference (EMBL)\nCM000618;\n"

    monkeypatch.setattr(UniProt, "search", mock_search)


@pytest.fixture
def mock_entrez_create_and_keep_cache(monkeypatch):
    """Mock remote service call to NCBI/Entrez for create_and_keep_cache test.

    Returns the result expected from NCBI with a call for each of the proteins
    in fixtures/sequences/human.fasta
    """

    responses = Mock()
    responses.side_effect = [
        io.BytesIO(_)
        for _ in (
            (
                b'<?xml version="1.0" encoding="UTF-8" ?>\n<!DOCTYPE eLinkResult '
                b'PUBLIC "-//NLM//DTD elink 20101123//EN" '
                b'"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20101123/elink.dtd">\n<eLinkResult>\n\n  '
                b"<LinkSet>\n    <DbFrom>protein</DbFrom>\n    <IdList>\n      <Id>10835069</Id>\n    "
                b"</IdList>\n    <LinkSetDb>\n      <DbTo>nuccore</DbTo>\n      "
                b"<LinkName>protein_nuccore</LinkName>\n      \n        "
                b"<Link>\n\t\t\t\t<Id>568815587</Id>\n\t\t\t</Link>\n        "
                b"<Link>\n\t\t\t\t<Id>197116381</Id>\n\t\t\t</Link>\n      \n    </LinkSetDb>\n  "
                b"</LinkSet>\n</eLinkResult>\n"
            ),
            (
                b'<?xml version="1.0" encoding="UTF-8" ?>\n<!DOCTYPE eLinkResult '
                b'PUBLIC "-//NLM//DTD elink 20101123//EN" '
                b'"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20101123/elink.dtd">\n<eLinkResult>\n\n  '
                b"<LinkSet>\n    <DbFrom>protein</DbFrom>\n    <IdList>\n      <Id>283135214</Id>\n    "
                b"</IdList>\n    <LinkSetDb>\n      <DbTo>nuccore</DbTo>\n      "
                b"<LinkName>protein_nuccore</LinkName>\n      \n        "
                b"<Link>\n\t\t\t\t<Id>1677498684</Id>\n\t\t\t</Link>\n        "
                b"<Link>\n\t\t\t\t<Id>568815587</Id>\n\t\t\t</Link>\n      \n    </LinkSetDb>\n  "
                b"</LinkSet>\n</eLinkResult>\n"
            ),
            (
                b'<?xml version="1.0" encoding="UTF-8" ?>\n<!DOCTYPE eLinkResult '
                b'PUBLIC "-//NLM//DTD elink 20101123//EN" '
                b'"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20101123/elink.dtd">\n<eLinkResult>\n\n  '
                b"<LinkSet>\n    <DbFrom>protein</DbFrom>\n    <IdList>\n      <Id>530397002</Id>\n    "
                b"</IdList>\n    <LinkSetDb>\n      <DbTo>nuccore</DbTo>\n      "
                b"<LinkName>protein_nuccore</LinkName>\n      \n        "
                b"<Link>\n\t\t\t\t<Id>767968522</Id>\n\t\t\t</Link>\n        "
                b"<Link>\n\t\t\t\t<Id>568815587</Id>\n\t\t\t</Link>\n      \n    </LinkSetDb>\n  "
                b"</LinkSet>\n</eLinkResult>\n"
            ),
            (
                b'<?xml version="1.0" encoding="UTF-8" ?>\n<!DOCTYPE eLinkResult '
                b'PUBLIC "-//NLM//DTD elink 20101123//EN" '
                b'"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20101123/elink.dtd">\n<eLinkResult>\n\n  '
                b"<LinkSet>\n    <DbFrom>protein</DbFrom>\n    <IdList>\n      <Id>283135242</Id>\n    "
                b"</IdList>\n    <LinkSetDb>\n      <DbTo>nuccore</DbTo>\n      "
                b"<LinkName>protein_nuccore</LinkName>\n      \n        "
                b"<Link>\n\t\t\t\t<Id>1677500256</Id>\n\t\t\t</Link>\n        "
                b"<Link>\n\t\t\t\t<Id>568815587</Id>\n\t\t\t</Link>\n      \n    "
                b"</LinkSetDb>\n  </LinkSet>\n</eLinkResult>\n"
            ),
        )
    ]

    monkeypatch.setattr(Entrez, "elink", responses)


def test_bad_infile(args_bad_infile):
    """ncfp stops if CLI input file does not exist.

    Validates that ncfp throws correct error, does not check output
    """
    with pytest.raises(NCFPException):
        ncfp.run_main(args_bad_infile)


def test_create_and_keep_cache(
    args_create_reuse_cache, mock_entrez_create_and_keep_cache
):
    """ncfp creates named cache from CLI and keeps it when rerunning.

    This calls ncfp twice and expects both calls run without error. Output
    is not checked.
    """
    # Create local cache
    ncfp.run_main(args_create_reuse_cache)
    # Reuse created cache
    ncfp.run_main(args_create_reuse_cache + ["--keepcache"])


def test_download_and_log(args_validate_and_log, mock_uniprot_download_and_log):
    """ncfp downloads coding sequences and logs output from CLI.

    Validates that ncfp runs without error, does not check output. Calls to UniProt
    are mocked for use with cloud CI.
    """
    ncfp.run_main(args_validate_and_log)
