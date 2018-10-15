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

from .caches import (
    has_nt_query,
    get_nt_query,
    has_aa_query,
    has_ncbi_uid,
    add_ncbi_uids,
    get_nt_noacc_uids,
    update_nt_uid_acc,
    add_gbheaders,
    add_gbfull,
    get_nogbhead_nt_uids,
    get_nogbfull_nt_acc,
    find_shortest_genbank,
)
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


def fetch_gb_headers(cachepath, retries, batchsize, disabletqdm=True):
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
    for batch in [
        nogbhead_uids[idx : idx + batchsize]
        for idx in range(0, len(nogbhead_uids), batchsize)
    ]:
        epost_histories.append(epost_history_with_retries(batch, "nucleotide", retries))
    for history in tqdm(
        epost_histories, desc="Fetching GenBank headers", disable=disabletqdm
    ):
        try:
            records = SeqIO.parse(
                efetch_history_with_retries(
                    history, "nucleotide", "gb", "text", retries
                ),
                "gb",
            )
            for record in records:
                taxonomy = " ".join(record.annotations["taxonomy"])
                addedrows.append(
                    add_gbheaders(
                        cachepath,
                        record.id,
                        len(record),
                        record.annotations["organism"],
                        taxonomy,
                        record.annotations["date"],
                    )
                )
        except NCFPMaxretryException:
            failcount += 1

    return addedrows, (failcount * batchsize)


def fetch_shortest_genbank(cachepath, retries, batchsize, disabletqdm=True):
    """Update cache with shortest full GenBank record for each input.

    cachepath     - path to cache database
    retries       - number of times to retry Entrez fetch
    batchsize     - number of GenBank records to fetch each time

    Checks the sequence length for each GenBank record associated with
    an input sequence, and records the shortest one. This list is used
    to batch fetch full GenBank records from Entrez. Every input
    sequence ends up with a corresponding nucleotide sequence.
    """
    addedrows = []
    failcount = 0

    # Get set of shortest GenBank accessions that cover input sequences
    # Identify those that need to be downloaded as full records
    shortids = find_shortest_genbank(cachepath)
    nogbfull_acc = get_nogbfull_nt_acc(cachepath)
    fetchaccs = list(shortids.intersection(nogbfull_acc))

    # Create batches and get EPost keys
    epost_histories = []
    for batch in [
        fetchaccs[idx : idx + batchsize] for idx in range(0, len(fetchaccs), batchsize)
    ]:
        epost_histories.append(epost_history_with_retries(batch, "nucleotide", retries))
    for history in tqdm(
        epost_histories, desc="Fetching full GenBank records", disable=disabletqdm
    ):
        try:
            records = SeqIO.parse(
                efetch_history_with_retries(
                    history, "nucleotide", "gbwithparts", "text", retries
                ),
                "gb",
            )
            for record in records:
                addedrows.append(add_gbfull(cachepath, record.id, record.format("gb")))
        except NCFPMaxretryException:
            failcount += 1

    return addedrows, (failcount * batchsize)


# Query NCBI singly with each record, to recover nucleotide accessions
def search_nt_ids(records, cachepath, retries, disabletqdm=True):
    """Queries NCBI nucleotide database and populates cache

    records - collection of SeqRecords
    cache   - path to cache
    retries - number of Entrez retries

    If the record's ID is in the cache, the ESearch is not
    performed - the cache is presumed to be up to date.

    The cache gets updated in the link table between the query
    sequence ID, and the nucleotide sequence accession from
    NCBI (seq_nt).

    If the record has a protein query (aa_query is populated in
    the cache's seqdata table), then seq_nt is populated from
    an ELink query; if the record has a nucleotide query
    (nt_query is populated, but aa_query is not), then a direct
    ESearch of NCBI's nucleotide databases is performed.
    """
    addedrows = []  # Holds list of added rows in nt_uid_acc
    noresult = 0  # Count of records with no result
    for record in tqdm(records, desc="Search NT IDs", disable=disabletqdm):
        if not has_ncbi_uid(cachepath, record.id) and has_nt_query(
            cachepath, record.id
        ):  # direct ESearch
            result = esearch_with_retries(
                get_nt_query(cachepath, record.id), "nucleotide", retries
            )
            if len(result["IdList"]):
                addedrows.extend(add_ncbi_uids(cachepath, record.id, result["IdList"]))
            else:
                noresult += 1
        elif not has_ncbi_uid(cachepath, record.id) and has_aa_query(
            cachepath, record.id
        ):  # ELink search
            result = elink_fetch_with_retries(
                record.id, "protein", "protein_nuccore", retries
            )
            try:
                idlist = [lid["Id"] for lid in result[0]["LinkSetDb"][0]["Link"]]
            except IndexError:  # No result returned - possible deleted record
                raise NCFPEFetchException("No link/record returned for %s" % record.id)
            if len(idlist):
                addedrows.extend(add_ncbi_uids(cachepath, record.id, idlist))
            else:
                noresult += 1
    return addedrows, noresult


# Update existing cache nt_uid_acc table with accessions from NCBI
def update_gb_accessions(cachepath, retries, disabletqdm=True):
    """Updates the cache table with GenBank accession for each UID.

    cachepath     - path to cache database
    retries       - number of Entres retries

    For each UID in nt_uid_acc where there is no GenBank accession,
    obtain the GenBank accession and update the row.
    """
    updatedrows = []
    noupdate = 0
    for uid in tqdm(
        get_nt_noacc_uids(cachepath), desc="Fetch UID accessions", disable=disabletqdm
    ):
        result = (
            efetch_with_retries(uid, "nucleotide", "acc", "text", retries)
            .read()
            .strip()
        )
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
        except Exception:
            tries += 1
    raise NCFPMaxretryException(
        "Query ID %s ESearch failed\n%s" % (query_id, last_exception())
    )


def elink_fetch_with_retries(query_id, dbname, linkdbname, maxretries):
    """Entrez ELink fetch for a single query_id.

    query_id        - query term for search
    dbname          - NCBI target database name for primary query
    linkdbname      - NCBI target database name for link query
    maxretries      - maximum number of attempts to make
    """
    tries = 0
    while tries < maxretries:
        try:
            matches = Entrez.read(
                Entrez.elink(dbfrom=dbname, linkname=linkdbname, id=query_id)
            )
            return matches
        except Exception:
            tries += 1
    raise NCFPMaxretryException(
        "Query ID %s ELink failed\n%s" % (query_id, last_exception())
    )


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
            data = Entrez.efetch(
                db=dbname, rettype=rettype, retmode=retmode, id=query_id
            ).read()
            if rettype in ["gb", "gbwithparts"] and retmode == "text":
                if not data.startswith("LOCUS"):
                    raise NCFPEFetchException(
                        "Data does not begin with expected LOCUS string"
                    )
            # Return data string as stream
            return StringIO(data)
        except Exception:
            tries += 1
    raise NCFPMaxretryException(
        "Query ID %s EFetch failed\n%s" % (query_id, last_exception())
    )


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
            history = Entrez.read(Entrez.epost(dbname, id=",".join(qids)))
            return history
        except Exception:
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
            data = Entrez.efetch(
                db=dbname,
                rettype=rettype,
                retmode=retmode,
                webenv=history["WebEnv"],
                query_key=history["QueryKey"],
            ).read()
            if rettype in ["gb", "gbwithparts"] and retmode == "text":
                if not data.startswith("LOCUS"):
                    raise NCFPEFetchException(
                        "Returned data does not begin with string LOCUS"
                    )
            return StringIO(data)
        except Exception:
            tries += 1
    raise NCFPMaxretryException("Failed to recover batch EFetch")


def set_entrez_email(email):
    """Sets Entrez email address."""
    Entrez.email = email
