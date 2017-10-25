#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Implements the ncbi_cds_from_protein script for getting nt sequences

(c) The James Hutton Institute 2017
Author: Leighton Pritchard

Contact: leighton.pritchard@hutton.ac.uk
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
import sys
import time

from collections import namedtuple

from Bio import SeqIO

from .parsers import parse_cmdline
from .. import __version__
from ..ncfp_tools import (last_exception, NCFPException)
from ..sequences import add_seqrecord_query
from ..caches import load_elink_cache
from ..entrez import (set_entrez_email, )


# Paths cache files for script
Cachepaths = namedtuple("Cachepaths", "elink acc gb gbfull")


# Process input sequences
def process_input_sequences(args, logger):
    """Process input sequences."""
    if args.infname is None or args.infname == '-':
        instream = sys.stdin
        logger.info("Reading sequences from stdin")
    else:
        if not os.path.isfile(args.infname):
            msg = "Input sequence file %s does not exist" % args.infname
            logger.error(msg)
            raise NCFPException(msg)
        try:
            instream = open(args.infname, 'r')
        except OSError:
            logger.error("Could not open input file %s", args.infname)
            logger.error(last_exception())
            sys.exit(1)
        logger.info("Reading sequences from %s", args.infname)
    try:
        records = list(SeqIO.parse(instream, 'fasta'))
    except IOError:
        logger.error("Could not parse sequence file %s", args.infname)
        logger.error(last_exception())
        sys.exit(1)
    logger.info("%d sequence records read successfully from %s",
                len(records), args.infname)
    return records


def build_logger(args):
    """Return a logger for this script.

    Instantiates a logger for the script, and adds basic info.
    """
    logger = logging.getLogger('ncfp: {}'.format(time.asctime))
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)

    # Verbose output?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # If a logfile was specified, use it
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, 'w')
        except OSError:
            logger.error('Could not open %s for logging', args.logfile)
            logger.error(last_exception())
            sys.exit(1)
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)

    # Report arguments
    args.cmdline = ' '.join(sys.argv)
    logger.info('Processed arguments: %s', args)
    logger.info('command-line: %s', args.cmdline)

    return logger


# Main script function
def run_main(namespace=None):
    """Run main process for ncfp script."""
    # Parse command-line if no namespace provided
    if namespace is None:
        args = parse_cmdline()
    else:
        args = namespace

    # Catch execution with no arguments
    if len(sys.argv) == 1:
        sys.stderr.write("ncbi_cds_from_protein " +
                         "version: {0}\n".format(__version__))
        return 0

    # Set up logging
    time0 = time.time()
    logger = build_logger(args)

    # Get input sequences and add query string
    logger.info("Parsing sequence input...")
    seqrecords = process_input_sequences(args, logger)

    # Set up cache filenames
    logger.info("Setting up data caches...")
    caches = Cachepaths(os.path.join(args.cachedir,
                                     'elink_{}'.format(args.cachestem)),
                        os.path.join(args.cachedir,
                                     'gb_{}'.format(args.cachestem)),
                        os.path.join(args.cachedir,
                                     'gbfull_{}'.format(args.cachestem)),
                        os.path.join(args.cachedir,
                                     'acc_{}'.format(args.cachestem)))
    logger.info(caches)

    # Rather than subclass the Biopython SeqRecord class, we add
    # the ncfp attribute. This will be a dictionary that holds key:values
    # for:
    # header_id - the accession pulled from the FASTA header
    # nt_query  - query term for NCBI nucleotide database
    # aa_query  - query term for NCBI protein database
    # These accessions are taken from the FASTA header and, if we can't
    # parse that appropriately, we can't search - so we skip those
    # sequences
    if args.uniprot:
        fmt = 'uniprot'
    else:
        fmt = 'ncbi'
    logger.info("Processing input sequences as %s format", fmt)
    qrecords = []
    qskipped = []
    for record in seqrecords:
        qrecord = add_seqrecord_query(record, fmt)
        if qrecord.ncfp['header_id'] is None:
            qskipped.append(qrecord)
        else:
            qrecords.append(qrecord)
    if len(qskipped):
        logger.warning("Skipped %d sequences (no header ID)",
                       len(qskipped))
        SeqIO.write(qskipped, args.skippedfname, 'fasta')
        logger.warning("Skipped sequences were written to %s",
                       args.skippedfname)
    logger.info("%d sequences taken forward with query",
                len(qrecords))

    # Set email address at Entrez
    set_entrez_email(args.email)

    
#    # Get Entrez EUtils protein_nuccore links for sequences
#    # Check local cache first
#    if not os.path.isfile(elink_cachename):
#        logger.info("No ELink cache... not checking")
#    else:
#        cachedata = load_elink_cache(elink_cachename)
#        logger.info("Checking ELink cache at %s", elink_cachename)
#        logger.info("Loaded data for %d queries", len(cachedata))
#    # Get EUtils links
#    for record in qrecords:
#        logger.info("Downloading protein_nuccore record for %s", record.id)
#        matches = elink_fetch_with_retries(record, "protein",
#                                           "protein_nuccore",
#                                           args.retries)
    
    # Report success
    logger.info('Completed. Time taken: %.3f',
                (time.time() - time0))
    return 0
