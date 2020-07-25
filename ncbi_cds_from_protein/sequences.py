#!/usr/bin/env python3
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
"""Functions for handling sequence data"""

import logging
import re
import sqlite3

from pathlib import Path

from Bio.SeqRecord import SeqRecord
from bioservices import UniProt
from tqdm import tqdm

from .caches import add_input_sequence, has_query

# regexes for parsing out Uniprot
re_uniprot_head = re.compile(r".*\|.*\|.*")
re_uniprot_gn = re.compile(r"(?<=GN=)[^\s]+")


def guess_seqtype(record):
    """Return best guess at SeqRecord origin on basis of header.

    :param record:  SeqRecord describing sequence.

    Guesses at the most likely sequence origin on the basis of sequence ID/header,
    as follows.

    If the description has a GN=* field, we guess UniProt.
    If the ID starts with "UPI" we guess UniParc
    Otherwise, we guess NCBI
    """
    logger = logging.getLogger(__name__)
    logger.debug("Guessing sequence type for %s...", record.id)

    if record.id.startswith("UPI"):
        logger.debug("...guessed UniParc")
        return "UniParc"
    else:  # Test for UniProt, see https://www.uniprot.org/help/fasta-headers
        match = re.search(re_uniprot_head, record.description)  # test for UniProt header

    if match is None:
        logger.debug("...guessed NCBI")
        return "NCBI"
    else:
        logger.debug("...guessed UniProt")
        return "UniProt"


# Process collection of SeqRecords into cache and skipped/kept
def process_sequences(records, cachepath: Path, disabletqdm: bool = True):
    """Triage SeqRecords into those that can/cannot be used

    This function also caches all inputs into the SQLite cache
    at cachepath

    :param records:  collection of SeqRecords
    :param cachepath: path to local sequence cache
    :param disabletqdm:  turn off tqdm progress bar
    """
    logger = logging.getLogger(__name__)
    logger.info("Processing sequences...")

    kept, skipped = [], []
    for record in tqdm(records, desc="Process input sequences", disable=disabletqdm):
        seqtype = guess_seqtype(record)
        if seqtype == "UniParc":
            logger.warning("Record %s looks like a UniParc cluster (skipping)", record.id)
            skipped.append(record)
            continue
        elif seqtype == "UniProt":
            match = re.search(re_uniprot_gn, record.description)
            if match is None:  #  No GN field
                logger.warning("Uniprot record %s has no GN field (skipping)", record.id)
                skipped.append(record)
                continue
            service = UniProt()
            result = service.search(match.group(0), columns="database(EMBL)")  # type: ignore
            qstring = result.split("\n")[1].strip()[:-1]
            try:  # Uniprot sequences are added to cache as (accession, NULL, nt_query)
                add_input_sequence(cachepath, record.id, None, qstring)
            except sqlite3.IntegrityError:  # Sequence exists
                continue
        elif seqtype == "NCBI":
            try:  # NCBI sequences are added to cache as (accession, aa_query, NULL)
                add_input_sequence(cachepath, record.id, record.id, None)
            except sqlite3.IntegrityError:  # Sequence exists
                continue
        # If the record has no query terms, skip it
        if has_query(cachepath, record.id):
            kept.append(record)
        else:
            skipped.append(record)
    return kept, skipped


# Extract a gene feature by locus tag
def extract_feature_by_locus_tag(record, tag, ftype="CDS"):
    """Returns the gene feature with passed tag from passed seqrecord.

    record      - Biopython SeqRecord
    tag         - locus tag to search for
    ftype       - feature types to search
    """
    for feature in [ftr for ftr in record.features if ftr.type == ftype]:
        try:
            if tag in feature.qualifiers["locus_tag"]:
                return feature
        except KeyError:
            continue
    return None


# Extract a gene feature by protein_id
def extract_feature_by_protein_id(record, tag, ftype="CDS"):
    """Returns the gene feature with passed tag from passed seqrecord.

    record      - Biopython SeqRecord
    tag         - locus tag to search for
    ftype       - feature types to search
    """
    for feature in [ftr for ftr in record.features if ftr.type == ftype]:
        try:
            if tag in feature.qualifiers["protein_id"]:
                return feature
        except KeyError:
            continue
    return None


# Extract the coding sequence from a feature
def extract_feature_cds(feature, record, stockholm):
    """Returns SeqRecord with CDS and translation of GenBank record feature.

    feature      - SeqFeature object
    record       - SeqRecord object
    stockholm    - List of (start, end) aa positions
    """
    # Extract nucleotide coding sequence
    ntseq = feature.extract(record.seq)

    # Account for offset start codon
    if "codon_start" in feature.qualifiers:
        startpos = int(feature.qualifiers["codon_start"][0]) - 1
    else:
        startpos = 0
    ntseq = ntseq[startpos:]

    # If Stockholm (start, end) headers were provided, trime sequences
    if len(stockholm) != 0:
        start, end = stockholm[0], stockholm[1]
        ntseq = ntseq[(start - 1) * 3 : (end * 3)]

    # Generate conceptual translation
    aaseq = ntseq.translate()
    if aaseq[-1] == "*":
        aaseq = aaseq[:-1]

    # Create SeqRecords of CDS and conceptual translation
    if "locus_tag" in feature.qualifiers:
        ntrecord = SeqRecord(seq=ntseq, description="coding sequence", id=feature.qualifiers["locus_tag"][0])
        aarecord = SeqRecord(seq=aaseq, description="conceptual translation", id=feature.qualifiers["locus_tag"][0])
    else:
        ntrecord = SeqRecord(seq=ntseq, description="coding sequence", id=feature.qualifiers["protein_id"][0])
        aarecord = SeqRecord(seq=aaseq, description="conceptual translation", id=feature.qualifiers["protein_id"][0])

    return ntrecord, aarecord
