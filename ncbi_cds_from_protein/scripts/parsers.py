# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) University of Strathclyde 2019-2022
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
# Copyright (c) 2019-2022 University of Strathclyde
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
"""Provides command-line/subcommand parsers for the ncfp script."""

import sys
import time

from pathlib import Path

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


# Process command-line for ncbi_cds_from_protein script
def parse_cmdline(args=None):
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="ncfp", formatter_class=ArgumentDefaultsHelpFormatter)

    # Add compulsory arguments
    parser.add_argument(
        action="store",
        dest="infname",
        default=None,
        help="input sequence file",
        type=Path,
    )
    parser.add_argument(
        action="store",
        dest="outdirname",
        default=Path("ncfp_output"),
        help="output directory for sequence files",
        type=Path,
    )
    parser.add_argument(
        action="store",
        dest="email",
        default=None,
        help="email address for NCBI/Entrez",
        type=str,
    )

    # Add options
    parser.add_argument(
        "-s",
        "--stockholm",
        dest="stockholm",
        action="store_true",
        default=False,
        help="parse Stockholm format domain regions",
    )
    parser.add_argument(
        "-d",
        "--cachedir",
        dest="cachedir",
        action="store",
        default=Path(".ncfp_cache"),
        type=Path,
        help="directory for cached data",
    )
    parser.add_argument(
        "-c",
        "--cachestem",
        dest="cachestem",
        action="store",
        default=time.strftime("%Y-%m-%d-%H-%m-%S"),
        type=str,
        help="suffix for cache filestems",
    )
    parser.add_argument(
        "-b",
        "--batchsize",
        dest="batchsize",
        action="store",
        default=100,
        type=int,
        help="batch size for EPost submissions",
    )
    parser.add_argument(
        "-r",
        "--retries",
        dest="retries",
        action="store",
        default=10,
        type=int,
        help="maximum number of Entrez retries",
    )
    parser.add_argument(
        "--limit",
        dest="limit",
        action="store",
        default=None,
        type=int,
        help="maximum number of sequences to process " + "(for testing)",
    )
    parser.add_argument(
        "--filestem",
        dest="filestem",
        action="store",
        default="ncfp",
        type=str,
        help="stem for output sequence files",
    )
    parser.add_argument(
        "--keepcache",
        dest="keepcache",
        action="store_true",
        default=False,
        help="keep download cache (for testing)",
    )
    parser.add_argument(
        "--skippedfile",
        dest="skippedfname",
        action="store",
        default="skipped.fasta",
        type=str,
        help="filename for skipped sequences",
    )
    parser.add_argument(
        "--use_protein_ids",
        dest="use_protein_ids",
        action="store_true",
        default=False,
        help=("Use protein_id field when recovering matched sequence information"),
    )
    parser.add_argument(
        "--unify_seqid",
        dest="unify_seqid",
        action="store_true",
        default=False,
        help=("Use the input sequence seqID for the recovered sequence"),
    )
    parser.add_argument(
        "--allow_alternative_start_codon",
        dest="alternative_start_codon",
        action="store_true",
        default=False,
        help=("Allow for an alternative start codon in the first sequence position"),
    )
    parser.add_argument(
        "-l",
        "--logfile",
        dest="logfile",
        action="store",
        default=None,
        type=Path,
        help="path to logfile",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="report verbosely",
    )
    parser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        default=False,
        help="report debug-level information",
    )
    parser.add_argument(
        "--disabletqdm",
        dest="disabletqdm",
        action="store_true",
        default=False,
        help="disable progress bar (for testing)",
    )

    # Parse arguments
    if args is None:
        args = sys.argv[1:]
    else:
        args = map(str, args)  # Ensure that args look like what we get from CLI
    return parser.parse_args(args)
