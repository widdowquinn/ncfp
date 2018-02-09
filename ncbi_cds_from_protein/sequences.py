#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions for handling sequence data

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

import re
import sqlite3

from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from .caches import (add_input_sequence, has_query)

# regexes for parsing out Uniprot
re_uniprot_gn = re.compile('(?<=GN=)[^\s]+')


# Process collection of SeqRecords into cache and skipped/kept
def process_sequences(records, cachepath, fmt='ncbi', disabletqdm=True):
    """Triage SeqRecords into those that can/cannot be used

    This function also caches all inputs into the SQLite cache
    at cachepath

    records      - collection of SeqRecords
    cachepath    - path to SQLite3 cache database
    fmt          - sequence format: ncbi or uniprot
    """
    kept, skipped = [], []
    for record in tqdm(records, desc="Process input sequences",
                       disable=disabletqdm):
        if fmt == 'uniprot':
            match = re.search(re_uniprot_gn, record.description)
            if match is None:
                qstring = None
            else:
                qstring = match.group(0)
            # Uniprot sequences are added to cache as
            # (accession, NULL, nt_query)
            try:
                add_input_sequence(cachepath, record.id, None, qstring)
            except sqlite3.IntegrityError:  # Sequence exists
                continue
        else:
            # NCBI sequences are added to cache as
            # (accession, aa_query, NULL)
            try:
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
def extract_feature_by_locus_tag(record, tag, ftype='CDS'):
    """Returns the gene feature with passed tag from passed seqrecord.

    record      - Biopython SeqRecord
    tag         - locus tag to search for
    ftype       - feature types to search
    """
    for feature in [ftr for ftr in record.features if
                    ftr.type == ftype]:
        try:
            if tag in feature.qualifiers['locus_tag']:
                return feature
        except KeyError:
            continue
    return None


# Extract a gene feature by protein_id
def extract_feature_by_protein_id(record, tag, ftype='CDS'):
    """Returns the gene feature with passed tag from passed seqrecord.

    record      - Biopython SeqRecord
    tag         - locus tag to search for
    ftype       - feature types to search
    """
    for feature in [ftr for ftr in record.features if
                    ftr.type == ftype]:
        try:
            if tag in feature.qualifiers['protein_id']:
                return feature
        except KeyError:
            continue
    return None


# Extract the coding sequence from a feature
def extract_feature_cds(feature, record, stockholm):
    """Returns SeqRecord with CDS and translation of GenBank record feature.

    feature      - SeqFeature object
    record       - SeqRecord object
    stockholm    - If not False, expects (start, end) aa positions
    """
    # Extract nucleotide coding sequence
    ntseq = feature.extract(record.seq)

    # Account for offset start codon
    if 'codon_start' in feature.qualifiers:
        startpos = int(feature.qualifiers['codon_start'][0]) - 1
    else:
        startpos = 0
    ntseq = ntseq[startpos:]

    # If Stockholm (start, end) headers were provided, trime sequences
    if stockholm:
        start, end = stockholm[0], stockholm[1]
        ntseq = ntseq[(start - 1) * 3:(end * 3)]

    # Generate conceptual translation
    aaseq = ntseq.translate()
    if aaseq[-1] == '*':
        aaseq = aaseq[:-1]

    # Create SeqRecords of CDS and conceptual translation
    if 'locus_tag' in feature.qualifiers:
        ntrecord = SeqRecord(seq=ntseq, description="coding sequence",
                             id=feature.qualifiers['locus_tag'][0])
        aarecord = SeqRecord(seq=aaseq, description="conceptual translation",
                             id=feature.qualifiers['locus_tag'][0])
    else:
        ntrecord = SeqRecord(seq=ntseq, description="coding sequence",
                             id=feature.qualifiers['protein_id'][0])
        aarecord = SeqRecord(seq=aaseq, description="conceptual translation",
                             id=feature.qualifiers['protein_id'][0])

    return ntrecord, aarecord
