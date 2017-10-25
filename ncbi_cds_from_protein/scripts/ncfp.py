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
import tqdm

from collections import namedtuple

from Bio import SeqIO

from .parsers import parse_cmdline
from .. import __version__
from ..ncfp_tools import (last_exception, NCFPException)
from ..sequences import add_seqrecord_query
from ..caches import (initialise_caches, load_cache, write_cache)
from ..entrez import (set_entrez_email, search_nt_ids, fetch_gb_headers)


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
    cachepaths = initialise_caches(args.cachedir, args.cachestem)
    logger.info(cachepaths)

    # Initialise caches
    

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

    # NCBI protein accessions can't be queried directly against the
    # nucleotide database, so we must perform an ELink search to
    # connect protein entries to nuccore entries, and populate the
    # .ncfp['nt_query'] attribute

    # UniProt sequence accessions should have a matching NCBI
    # nucleotide entry, that we need to identify with an ESearch
    # and put in the .ncfp['nt_query'] slot
    logger.info("Identifying nucleotide accessions...")
    acc_cache = load_cache(cachepaths.acc)
    logger.info("Accession cache %s contains %d entries",
                cachepaths.acc, len(acc_cache))
    qrecords, acc_cache, accfail = search_nt_ids(qrecords, acc_cache, args.retries)
    logger.info("Writing accession cache with %d entries to %s",
                len(acc_cache), cachepaths.acc)
    write_cache(acc_cache, cachepaths.acc)
    logger.info("NCBI nucleotide accessions retrieved for %d records (%d failed)",
                len(qrecords), len(accfail))

    # At this point, all records should have a .ncfp['nt_acc']
    # attribute, and be queryable against the NCBI nuccore db.
    # First, we associate each of the records with the accession for
    # its corresponding GenBank file.
    
    # Then, we retrieve GenBank headers, for inspection and
    # identification of the most useful nucleotide record.
    logger.info("Downloading GenBank headers for each record")
    gb_cache = load_cache(cachepaths.gb)
    logger.info("GenBank header cache %s contains %d entries",
                cachepaths.gb, len(gb_cache))
    gb_cache = fetch_gb_headers(qrecords, gb_cache,
                                args.batchsize, args.retries)
    logger.info("Writing GenBank header cache with %d entries to %s",
                len(gb_cache), cachepaths.gb)
    write_cache(gb_cache, cachepaths.gb)

    # Next we recover the complete GenBank records for useful
    # sequences
    pass

    # Report success
    logger.info('Completed. Time taken: %.3f',
                (time.time() - time0))
    return 0
