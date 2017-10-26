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

from Bio import (Entrez, SeqIO)
from tqdm import tqdm

from .caches import (has_nt_query, get_nt_query,
                     has_ncbi_uid, add_ncbi_uids,
                     get_nt_noacc_uids, update_nt_uid_acc,
                     get_nt_uids, add_gb_headers,
                     get_nogbhead_nt_uids)
from .ncfp_tools import last_exception

# EXCEPTIONS
# ==========


class NCFPELinkException(Exception):
    """Exception for ELink qeries."""

    def __init__(self, msg="Error in ncfp ELink query"):
        """Instantiate class."""
        Exception.__init__(self, msg)


class NCFPEPostException(Exception):
    """Exception for EPost qeries."""

    def __init__(self, msg="Error in ncfp EPost query"):
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


def fetch_gb_headers(cachepath, retries, batchsize):
    """Update cache with NCBI GenBank headers for passed records

    cachepath - path to cache
    retries   - number of Entrez retries
    batchsize - number of uids per EPost query

    Gets list of UIDs with no existing cached GenBank headers, and
    batch EFetches the GenBank headers.
    """
    addedrows = []
    failcount = 0
    # Create batches and get EPost keys for each batch
    epost_histories = []
    nogbhead_uids = get_nogbhead_nt_uids(cachepath)
    for batch in [nogbhead_uids[idx:idx + batchsize] for idx in
                  range(0, len(nogbhead_uids), batchsize)]:
        epost_histories.append(epost_history_with_retries(batch,
                                                          'nucleotide',
                                                          retries))
    for history in tqdm(epost_histories,
                        desc="Fetching GenBank headers"):
        print(history)
        try:
            records = SeqIO.parse(efetch_history_with_retries(history,
                                                              'nucleotide',
                                                              'gb', 'text',
                                                              retries), 'gb')
            print(records)
            for record in records:
                taxonomy = ' '.join(record.annotations['taxonomy'])
                addedrows.append(add_gb_headers(cachepath,
                                                record.id, len(record),
                                                record.annotations['organism'],
                                                taxonomy,
                                                record.annotations['date']))
        except NCFPMaxretryException:
            failcount += 1
    return addedrows, failcount


# Query NCBI singly with each record, to recover nucleotide accessions
def search_nt_ids(records, cachepath, retries):
    """Queries NCBI nucleotide database and populates cache

    records - collection of SeqRecords
    cache   - path to cache
    retries - number of Entrez retries

    If the record's ID is in the cache, the ESearch is not
    performed - the cache is presumed to be up to date.
    """
    addedrows = []  # Holds list of added rows in nt_uid_acc
    noresult = 0    # Count of records with no result
    for record in tqdm(records, desc="Search NT IDs"):
        if (not has_ncbi_uid(cachepath, record.id) and
                has_nt_query(cachepath, record.id)):
            result = esearch_with_retries(get_nt_query(cachepath,
                                                       record.id),
                                          'nucleotide', retries)
            if len(result['IdList']):
                addedrows.extend(add_ncbi_uids(cachepath, record.id,
                                               result['IdList']))
            else:
                noresult += 1
    return addedrows, noresult


# Update existing cache nt_uid_acc table with accessions from NCBI
def update_gb_accessions(cachepath, retries):
    """Updates the cache table with GenBank accession for each UID.

    cachepath     - path to cache database
    retries       - number of Entres retries

    For each UID in nt_uid_acc where there is no GenBank accession,
    obtain the GenBank accession and update the row.
    """
    updatedrows = []
    noupdate = 0
    for uid in tqdm(get_nt_noacc_uids(cachepath),
                    desc="Fetch UID accessions"):
        result = efetch_with_retries(uid, 'nucleotide', 'acc',
                                     'text', retries).read().strip()
        if result is None:
            noupdate += 1
        else:
            updatedrows.extend(update_nt_uid_acc(cachepath, uid, result))
    return updatedrows, noupdate


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
    raise NCFPMaxretryException("Query ID %s ESearch failed" % query_id)


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
            print(last_exception())
            tries += 1
    raise NCFPMaxretryException("Query ID %s EFetch failed" % query_id)


# Batch EPost to get history for a set of query IDs
def epost_history_with_retries(qids, dbname, maxretries):
    """Entrez EPost for a batch of queryids.

    qids            - Collection of query IDs
    dbname          - target NCBI database
    maxretries      - maximum download attempts

    Returns the generated history, as parsed by Entrez.read()
    """
    tries = 0
    while tries < maxretries:
        try:
            history = Entrez.read(Entrez.epost(dbname,
                                               id=','.join(qids)))
            return history
        except:
            tries += 1
    raise NCFPMaxretryException("Batch EPost failed.")


def efetch_history_with_retries(history, dbname, rettype, retmode, maxretries):
    """Run Entrez EFetch on a passed EPost history

    history          - EPost history
    dbname           - target NCBI database
    rettype          - return datatype
    retmode          - return data mode
    maxretries       - maximum download attempts

    Returns the result data as a tokenised string.
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
    raise NCFPMaxretryException("Failed to recover batch EFetch")


def set_entrez_email(email):
    """Sets Entrez email address."""
    Entrez.email = email
