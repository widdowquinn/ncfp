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


# regexes for parsing out query strings
re_uniprot = re.compile('(?<=GN=)[^\s]+')


# Add a .query_string attribute to a Biopython SeqRecord object,
def add_seqrecord_query(record, fmt="ncbi"):
    """Adds .ncfp attribute to a SeqRecord

    The .ncfp attribute is a dictionary that holds key:values

    header_id - the accession pulled from the FASTA header
    nt_query  - query term for NCBI nucleotide database
    aa_query  - query term for NCBI protein database

    The header ID is dependent on the input sequence source.
    NCBI proteins have the accession as their identifier. UniProt
    proteins have the originating gene sequence accession as
    a 'GN=' gene name field in the description.
    """
    # Identify query string
    if fmt == 'uniprot':
        match = re.search(re_uniprot, record.description)
        if match is None:
            qstring = None
        else:
            qstring = match.group(0)
    else:
        qstring = record.id

    # Add query string to record
    record.ncfp = {'header_id': qstring}

    return record
