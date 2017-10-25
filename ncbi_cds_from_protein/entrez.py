#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions to interact with Entrez

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

from io import StringIO

from Bio import Entrez, SeqIO
from tqdm import tqdm

from .ncfp_tools import last_exception

# EXCEPTIONS
# ==========
class NCFPELinkException(Exception):
    """Exception for ELink qeries."""

    def __init__(self, msg="Error in ncfp ELink query"):
        """Instantiate class."""
        Exception.__init__(self, msg)


class NCFPEFetchException(Exception):
    """Exception for EFetch queries."""

    def __init__(self, msg="Error in ncfp EFetch query"):
        """Instantiate class."""
        Exception.__init__(self, msg)


class NCFPESearchException(Exception):
    """Exception for ESearch queries."""

    def __init__(self, msg="Error in ncfp ESearch query"):
        """Instantiate class."""
        Exception.__init__(self, msg)


class NCFPMaxretryException(Exception):
    """Exception raised when maximum retries are reached"""

    def __init__(self, msg="Maximum retries exceeded"):
        """Instantiate class."""
        Exception.__init__(self, msg)


def fetch_gb_headers(records, cache, batchsize, retries):
    """Returns NCBI GenBank headers for passed records

    records    - collection of SeqRecords
    cache      - gb_cache for GenBank headers
    batchsize  - the batch size for EPost history
    retries    - number of Entrez retries
    """
    # Collate all unique nucleotide accessions
    accessions = set()
    for record in tqdm(records, desc="Collating accessions"):
        for acc in record.ncfp['nt_acc']:
            accessions.add(acc)

    # Check accessions against those in the cache. We only
    # retrieve accessions not in the cache
    fetchlist = list(accessions.difference(set(cache.keys())))

    # Batch fetchlist for separate EPost submissions
    histories = []
    for batch in [fetchlist[i:i+batchsize] for i in
                  range(0, len(fetchlist), batchsize)]:
        histories.append(epost_history_with_retries(batch,
                                                    'nucleotide',
                                                    retries))

    # Fetch headers with each history
    for history in tqdm(histories, desc="Fetching GenBank headers"):
        gbheader = efetch_history_with_retries(history, 'nucleotide',
                                               'gb', 'text', retries)
        # Parse the gbheader data and get the record ID, length
        # and annotations
        print(accessions)
        for record in SeqIO.parse(gbheader, 'gb'):
            print(record.id, len(record), record.annotations['accessions'])
    return cache


# Query NCBI singly with each record, to recover nucleotide accessions
def search_nt_ids(records, cache, retries):
    """Returns records with NCBI nucleotide ID as an attribute.

    records - collection of SeqRecords
    cache   - accession cache (dictionary)
    retries - number of Entrez retries

    Passed records must have the record.ncfp attribute. This will
    be extended with the ['nt_accession'] key, whose value is the
    NCBI accession for that sequence.

    If the record's ID is in the cache, the ESearch is not
    performed - the cache is presumed to be up to date.
    """
    successrecords = []
    failrecords = []
    for record in tqdm(records, desc="Search NT IDs"):
        if record.ncfp['nt_query'] not in cache:
            try:
                result = esearch_with_retries(record.ncfp['nt_query'],
                                              'nucleotide', retries)
            except NCFPESearchException:
                failrecords.append(record)
            if int(result['Count']) == 0:
                failrecords.append(record)
            else:
                successrecords.append(record)
                record.ncfp['nt_acc'] = result['IdList']
                cache[record.ncfp['nt_query']] = result['IdList']
        else:
            record.ncfp['nt_acc'] = cache[record.ncfp['nt_query']]
            successrecords.append(record)
    return successrecords, cache, failrecords


# Run an ESearch on a single ID
def esearch_with_retries(query_id, dbname, maxretries):
    """Entrez ESearch for single query ID.

    query_id    - query term for search
    dbname      - NCBI target database name
    maxretries  - maximum number of attempts to make

    Returns the parsed record resulting from the ESearch,
    after trying the ESearch up to a maximum number of times.
    """
    tries = 0
    while tries < maxretries:
        try:
            handle = Entrez.esearch(db=dbname, term=query_id)
            return Entrez.read(handle)
        except:
            tries += 1
    raise NCFPMaxretryException("Query ID %s ESearch failed\n%s" %
                                (query_id, last_exception()))


# Run an EFetch on a single ID
def efetch_with_retries(query_id, dbname, rettype, retmode, maxretries):
    """Entrez EFetch for a single query ID.

    query_id      - query term for search
    dbname        - NCBI target database name
    rettype       - requested return type
    retmode       - format of data to be returned
    maxretries    - maximum number of attempts to make

    Returns a handle to a completely buffered string, after a read()
    operation, sanity check, and retokenising.
    """
    tries = 0
    while tries < maxretries:
        try:
            data = Entrez.efetch(db=dbname, rettype=rettype,
                                 retmode=retmode,
                                 id=query_id).read()
            if rettype in ['gb', 'gbwithparts'] and retmode == 'text':
                assert data.startswith('LOCUS')
            # Return data string as stream
            return StringIO(data)
        except:
            tries += 1
    raise NCFPMaxretryException("Query ID %s EFetch failed\n" %
                                (query_id, last_exception()))


def epost_history_with_retries(batch, dbname, maxretries):
    """Returns EPost history key for query batch

    batch      - collection of query terms
    dbname     - NCBI target database
    maxretries - maximum number of query attempts
    """
    tries = 0
    while tries < maxretries:
        try:
            result = Entrez.epost(id=','.join(batch), db=dbname)
            return Entrez.read(result)
        except:
            tries += 1
    raise NCFPMaxretryException("EPost history failed\n%s" % last_exception())


def efetch_history_with_retries(history, dbname, rettype, retmode, maxretries):
    """Return batched EFetch as tokenised string

    history         - NCBI EPost history key
    dbname          - NCBI target database
    rettype         - data return type
    retmode         - data return mode
    maxretries      - maximum retrieval attempts
    """
    tries = 0
    while tries < maxretries:
        try:
            data = Entrez.efetch(db=dbname, rettype=rettype, retmode=retmode,
                                 webenv=history['WebEnv'],
                                 query_key=history['QueryKey']).read()
            if rettype in ['gb', 'gbwithparts'] and retmode == 'text':
                assert data.startswith('LOCUS')
            return StringIO(data)
        except:
            tries += 1
    raise NCFPMaxretryException("EFetch history failed\n%s" % last_exception())


def set_entrez_email(email):
    """Sets Entrez email address."""
    Entrez.email = email
