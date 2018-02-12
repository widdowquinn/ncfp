# -*- coding: utf-8 -*-
"""Provides command-line/subcommand parsers for the ncfp script.

(c) The James Hutton Institute 2016-2017
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

Copyright (c) 2016-2017 The James Hutton Institute

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

import sys
import time

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


# Process command-line for ncbi_cds_from_protein script
def parse_cmdline(args=None):
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="ncfp",
                                 formatter_class=ArgumentDefaultsHelpFormatter)

    # Add compulsory arguments
    parser.add_argument(action='store', dest='infname',
                        default=None,
                        help='input sequence file')
    parser.add_argument(action='store', dest='outdirname',
                        default=None,
                        help='output directory for sequence files')
    parser.add_argument(action='store', dest='email',
                        default=None,
                        help='email address for NCBI/Entrez')

    # Add options
    parser.add_argument('-u', '--uniprot', dest='uniprot',
                        action='store_true', default=False,
                        help='treat input sequences as UniProt FASTA')
    parser.add_argument('-s', '--stockholm', dest='stockholm',
                        action='store_true', default=False,
                        help='parse Stockholm format sequence regions')
    parser.add_argument('-d', '--cachedir', dest='cachedir',
                        action='store', default='.ncfp_cache',
                        help='directory for cached data')
    parser.add_argument('-c', '--cachestem', dest='cachestem',
                        action='store',
                        default=time.strftime("%Y-%m-%d-%H-%m-%S"),
                        help='suffix for cache filestems')
    parser.add_argument('-b', '--batchsize', dest='batchsize',
                        action='store', default=100, type=int,
                        help='batch size for EPost submissions')
    parser.add_argument('-r', '--retries', dest='retries',
                        action='store', default=10, type=int,
                        help='maximum number of Entrez retries')
    parser.add_argument('--limit', dest='limit',
                        action='store', default=None, type=int,
                        help='maximum number of sequences to process ' +
                        '(for testing)')
    parser.add_argument('--filestem', dest='filestem',
                        action='store', default='ncfp',
                        help='stem for output sequence files')
    parser.add_argument('--keepcache', dest='keepcache',
                        action='store_true', default=False,
                        help='keep download cache (for testing)')
    parser.add_argument('--skippedfile', dest='skippedfname',
                        action='store', default='skipped.fasta',
                        help='path to file with skipped sequences')
    parser.add_argument('-l', '--logfile', dest='logfile',
                        action='store', default=None,
                        help='path to logfile')
    parser.add_argument('-v', '--verbose', dest='verbose',
                        action='store_true', default=False,
                        help='report verbose progress')
    parser.add_argument('--disabletqdm', dest='disabletqdm',
                        action='store_true', default=False,
                        help='disable progress bar (for testing)')

    # Parse arguments
    if args is None:
        args = sys.argv[1:]
    return parser.parse_args(args)
