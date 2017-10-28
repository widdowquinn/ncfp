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

import os
import re
import sys
import time

from io import StringIO

from Bio import SeqIO

from .parsers import parse_cmdline
from .logger import build_logger
from .. import __version__
from ..ncfp_tools import (last_exception, NCFPException)
from ..sequences import (process_sequences, re_uniprot_gn,
                         extract_feature_by_locus_tag,
                         extract_feature_cds)
from ..caches import (initialise_dbcache, find_record_cds)
from ..entrez import (set_entrez_email, search_nt_ids,
                      update_gb_accessions,
                      fetch_gb_headers,
                      fetch_shortest_genbank)


# Process input sequences
def load_input_sequences(args, logger):
    """Load input FASTA sequences."""
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
    logger = build_logger('ncfp', args)

    # Set email address at Entrez
    set_entrez_email(args.email)

    # Make sure we can write to the output directory
    try:
        os.makedirs(args.outdirname, exist_ok=True)
    except OSError:
        logger.error("Could not use/create output directory %s (exiting)",
                     args.outdirname)
        logger.error(last_exception)
        raise SystemExit(1)
    
    # Initialise cache
    cachepath = os.path.join(args.cachedir,
                             'ncfpcache_%s.sqlite3' %
                             args.cachestem)
    # Use the old SQLite3 database if --keepcache set
    if (args.keepcache and os.path.isfile(cachepath)):
        logger.info("Not overwriting old cache at %s...",
                    cachepath)
    else:
        logger.info("Setting up SQLite3 database cache at %s...",
                    cachepath)
        initialise_dbcache(cachepath)

    # Get input sequences
    logger.info("Parsing sequence input...")
    seqrecords = load_input_sequences(args, logger)

    # Rather than subclass the Biopython SeqRecord class, we use the cache
    # database and query it to obtain relevant query terms and NCBI database
    # matches.
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
    qrecords, qskipped = process_sequences(seqrecords, cachepath, fmt)
    if len(qskipped):
        logger.warning("Skipped %d sequences (no query term found)",
                       len(qskipped))
        SeqIO.write(qskipped, args.skippedfname, 'fasta')
        logger.warning("Skipped sequences were written to %s",
                       args.skippedfname)
    logger.info("%d sequences taken forward with query",
                len(qrecords))
    if len(qrecords) == 0:
        logger.warning("No new input sequences were found! (in cache?)")

    # NCBI protein accessions can't be queried directly against the
    # nucleotide database, so we must perform an ELink search to
    # connect protein entries to nuccore entries, and populate the
    # .ncfp['nt_query'] attribute
    if not args.uniprot:
        raise NotImplementedError("Only Uniprot inputs for now")

    # UniProt sequence accessions should have at least one matching
    # NCBI nucleotide entry, that we need to identify with an ESearch
    # and put in the cache.
    logger.info("Identifying nucleotide accessions...")
    addedrows, countfail = search_nt_ids(qrecords, cachepath, args.retries)
    logger.info("Added %d new UIDs to cache", len(addedrows))
    if countfail:
        logger.warning("NCBI nucleotide accession search failed for " +
                       "%d records", countfail)
    if len(addedrows) == 0 and countfail == 0:
        logger.warning(
            "No nucleotide accession downloads were required! (in cache?)")

    # At this point, we want to retrieve all records that are
    # queryable against the NCBI nt database.
    # First, we associate GenBank accessions with a UID. This can be done
    # without reference to the records, using only the cache.
    logger.info("Collecting GenBank accessions...")
    updatedrows, countfail = update_gb_accessions(cachepath, args.retries)
    logger.info("Updated GenBank accessions for %d UIDs", len(updatedrows))
    if countfail:
        logger.warning("Unable to update GenBank accessions for %d UIDs",
                       countfail)
    if len(updatedrows) == 0 and countfail == 0:
        logger.warning(
            "No GenBank accession downloads were required! (in cache?)")

    # Next we recover GenBank headers and extract useful information -
    # sequence length, taxonomy, and so on.
    logger.info("Fetching GenBank headers...")
    addedrows, countfail = fetch_gb_headers(cachepath,
                                            args.retries, args.batchsize)
    logger.info("Fetched GenBank headers for %d UIDs", len(addedrows))
    if countfail:
        logger.warning("Unable to update GenBank headers for %d UIDs",
                       countfail)
    if len(addedrows) == 0 and countfail == 0:
        logger.warning(
            "No GenBank header downloads were required! (in cache?)")

    # Next we recover the shortest complete GenBank record for each input
    # sequence
    logger.info("Fetching shortest complete GenBank records...")
    addedrows, countfail = fetch_shortest_genbank(cachepath,
                                                  args.retries, args.batchsize)
    logger.info("Fetched GenBank records for %d UIDs", len(addedrows))
    if countfail:
        logger.warning("Unable to get complete GenBank files for %d UIDs",
                       countfail)
    if len(addedrows) == 0 and countfail == 0:
        logger.warning(
            "No complete GenBank downloads were required! (in cache?)")

    # Now that all the required GenBank nucleotide information is in the
    # local cache, we extract the CDS for each of the input sequences
    logger.info("Extracting CDS for each input sequence...")
    nt_sequences = []  # Holds extracted nucleotide sequences
    for record in seqrecords:
        result = find_record_cds(cachepath, record.id)
        if len(result) == 0:
            logger.warning("No record found for sequence input %s", record.id)
        elif len(result) > 1:
            logger.error("More than one record returned for %s (exiting)",
                         record.id)
            raise SystemExit(1)
        else:
            gbrecord = SeqIO.read(StringIO(result[0][-1]), 'gb')
            logger.info("Sequence %s matches GenBank entry %s",
                        record.id, gbrecord.id)
            # For Uniprot sequences, we need to extract the gene name
            if args.uniprot:
                match = re.search(re_uniprot_gn, record.description)
                gene_name = match.group(0)
            else:
                raise NotImplementedError
            # Get the matching CDS
            logger.info("Searching for CDS: %s", gene_name)
            feature = extract_feature_by_locus_tag(gbrecord, gene_name)
            if feature is None:
                logger.info("Could not identify CDS feature for %s", record.id)
            else:
                logger.info("\tSequence %s matches CDS feature %s", record.id,
                            feature.qualifiers['protein_id'][0])
                logger.info("\tExtracting coding sequence...")
                ntseq, aaseq = extract_feature_cds(feature, gbrecord,
                                                   args.stockholm)
                if aaseq.seq == record.seq:
                    logger.info("\t\tTranslated sequence matches input sequence")
                    nt_sequences.append((record, ntseq))
                else:
                    logger.warning("\t\tTranslated sequence does not match input sequence!")

    # Write matching nucleotide and protein sequences to file
    logger.info("Matched %d/%d records", len(nt_sequences),
                len(seqrecords))
    for (record, cds) in nt_sequences:
        logger.info("\t%-40s to CDS: %s", record.id, cds.id)
        
    # Write input sequences that were matched
    aafilename = os.path.join(args.outdirname,
                              '_'.join([args.filestem, 'aa.fasta']))
    logger.info("Writing matched input sequences to %s", aafilename)
    SeqIO.write([aaseq for (aaseq, ntseq) in nt_sequences],
                aafilename, 'fasta')

    # Write coding sequences
    ntfilename = os.path.join(args.outdirname,
                              '_'.join([args.filestem, 'nt.fasta']))
    logger.info("Writing matched output sequences to %s", ntfilename)
    SeqIO.write([ntseq for (aaseq, ntseq) in nt_sequences],
                ntfilename, 'fasta')

    # Report success
    logger.info('Completed. Time taken: %.3f',
                (time.time() - time0))
    return 0
