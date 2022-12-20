#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) University of Strathclyde 2019-2022
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
# Copyright (c) 2019-2022 University of Strathclyde
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
from tqdm.auto import tqdm

import ncbi_cds_from_protein

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
        match = re.search(
            re_uniprot_head, record.description
        )  # test for UniProt header

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

    u_service = UniProt()

    for record in tqdm(
        records, desc="1/5 Process input sequences", disable=disabletqdm
    ):
        seqtype = guess_seqtype(record)
        if seqtype == "UniParc":
            logger.warning(
                "Record %s looks like a UniParc cluster (skipping)", record.id
            )
            skipped.append(record)
            continue
        elif seqtype == "UniProt":
            match = re.search(re_uniprot_gn, record.description)
            pstring = None  # None if we don't need to look up GeneID
            if match is None:  #  No GN field
                logger.warning(
                    "Uniprot record %s has no GN field (skipping)", record.id
                )
                skipped.append(record)
                continue
            logger.debug(
                "Uniprot record %s has GN field: %s", record.id, match.group(0)
            )
            # The UniProt API was updated in June 2022, requiring a change
            # to the returned field for the cross-reference to EMBL
            # Get the UniProt ID from the accession and use this to query the API
            query_acc = record.id.split("|")[1]
            logger.debug(
                "Querying UniProt with %s to match xref_embl",
                query_acc,
            )
            # Use the UniProt record ID as the query to the EMBL ID, retrieving
            # the EMBL record accession, and the ORF gene name
            result = u_service.search(
                query_acc, columns="xref_embl")  # type: ignore
            qstring = result.split("\n")[1].strip()[:-1]
            result = u_service.search(
                query_acc, columns="gene_orf")  # type: ignore
            pstring = result.split("\n")[1].strip()
            # print(qstring, pstring)
            if qstring == "":
                logger.warning(
                    "Uniprot record %s has no EMBL cross-reference",
                    record.id,
                )
                # If there is no EMBL cross-reference, we can try the GeneID as
                # a fallback. This refers to a gene database record at NCBI, which
                # then must be queried to get the nucleotide cross-reference.
                logger.debug("Trying xref_geneid as last resort")
                gresult = u_service.search(
                    query_acc, columns="xref_geneid")  # type: ignore
                gstring = gresult.split("\n")[1].strip()[:-1]
                if gstring == "":
                    logger.warning(
                        "Uniprot record %s has no GeneID cross-reference (skipping)",
                        record.id,
                    )
                    continue
                # Fetch the GeneID entry and extract nucleotide information from it.
                logger.debug(
                    "Finding nucleotide entry from GeneID %s", gstring)
                handle = ncbi_cds_from_protein.entrez.efetch_with_retries(
                    gstring, "gene", "acc", "text", 10
                )  # NOTE: hard-coded retry count
                acc = [
                    _.strip().split()[1]
                    for _ in handle.readlines()
                    if _.strip().startswith("Annotation")
                ]
                if len(acc) == 0:
                    logger.warning(
                        "No nucleotide entry found for %s (skipping)", gstring
                    )
                    continue
                qstring = acc[0]  # This is our new query string
                logger.debug("Found nucleotide entry %s", qstring)
                # We also need the xref_refseq entry to get a protein ID for
                # the query.
                logger.debug("Finding protein entry from GeneID %s", gstring)
                presult = u_service.search(
                    gstring, columns="xref_refseq")  # type: ignore
                pstring = presult.split("\n")[1].strip()[:-1]
                if pstring == "":
                    logger.warning(
                        "Could not identify RefSeq protein ID for %s (skipping)",
                        gstring,
                    )
                    continue
                logger.debug("Found RefSeq protein IDs %s", pstring)
            logger.debug("Recovered NCBI database accession: %s", qstring)
            logger.debug("Accession has ORF gene name: %s", pstring)
            # UniProt can return multiple UIDs separated by semicolons. Sometimes the same
            # UID is repeated. However, the current cache schema uses the accession as primary
            # key in the same table as the query IDs.
            # TODO: Update schema to allow multiple queries per record
            for qid in qstring.split(";"):
                logger.debug(
                    "Adding record %s to cache with query %s", record.id, qid)
                # Uniprot sequences are added to cache as (accession, NULL, nt_query)
                try:
                    logger.debug(
                        "Adding %s as nt ID (protein ID: %s)", qid, pstring)
                    add_input_sequence(cachepath, record.id, pstring, qid)
                except sqlite3.IntegrityError:  # Sequence exists
                    logger.warning(
                        "Additional query terms found for %s: %s (not used)",
                        record.id,
                        qid,
                    )
                    continue
        elif seqtype == "NCBI":
            # NCBI sequences are added to cache as (accession, aa_query, NULL)
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
def extract_feature_by_locus_tag(record, tag, ftype="CDS"):
    """Returns the gene feature with passed tag from passed seqrecord.

    record      - Biopython SeqRecord
    tag         - locus tag to search for
    ftype       - feature types to search
    """
    logger = logging.getLogger(__name__)

    for feature in [ftr for ftr in record.features if ftr.type == ftype]:
        try:
            if tag in feature.qualifiers["locus_tag"]:
                logger.debug("Found %s in locus_tag", tag)
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
    logger = logging.getLogger(__name__)

    for feature in [ftr for ftr in record.features if ftr.type == ftype]:
        try:
            if tag in feature.qualifiers["protein_id"]:
                logger.debug("Found %s in protein_id", tag)
                return feature
        except KeyError:
            continue
    return None


# Extract a gene feature by gene_id
def extract_feature_by_gene_id(record, tag, ftype="CDS"):
    """Returns the gene feature with passed tag from passed seqrecord.

    record      - Biopython SeqRecord
    tag         - gene ID to search for
    ftype       - feature types to search
    """
    logger = logging.getLogger(__name__)

    for feature in [ftr for ftr in record.features if ftr.type == ftype]:
        try:
            if tag in feature.qualifiers["gene"]:
                logger.debug("Found %s in gene", tag)
                return feature
        except KeyError:
            continue
    return None


# Extract the coding sequence from a feature
def extract_feature_cds(feature, record, stockholm, args):
    """Returns SeqRecord with CDS and translation of GenBank record feature.

    :param feature:  SeqFeature object
    :param record:  SeqRecord object
    :param stockholm:  Tuple with (start, end) aa positions
    :param args:  CLI script arguments
    """
    logger = logging.getLogger(__name__)
    logger.debug("Extracting CDS from feature %s", feature)

    # Extract nucleotide coding sequence from GenBank record
    ntseq = feature.extract(record.seq)

    # Account for offset start codon, if necessary (not usually an issue)
    if "codon_start" in feature.qualifiers:
        startpos = int(feature.qualifiers["codon_start"][0]) - 1
    else:
        startpos = 0
    ntseq = ntseq[startpos:]

    # If Stockholm (start, end) headers were provided, unpack tuple and trim sequences
    if len(stockholm) != 0:
        start, end = stockholm
        logger.debug(
            "Trimming CDS to Stockholm coordinates: %d..%d", start, end)
        ntseq = ntseq[(start - 1) * 3: (end * 3)]

    # Has a region been requested that is outside the CDS?
    # NOTE: This can happen if a Stockholm region is provided, but
    # the 'wrong' CDS is selected. A route out of this may be to
    # iterate over all matching CDS and choose one (the first?)
    # that provides the correct conceptual translation.
    if not len(ntseq):
        logger.warning(
            "Requested region %d..%d is outside CDS, skipping", start, end)
        return None, None

    # Generate conceptual translation from extracted nucleotide sequence
    aaseq = ntseq.translate()
    if aaseq[-1] == "*":
        aaseq = aaseq[:-1]

    # Create SeqRecords of CDS and conceptual translation
    field = "protein_id" if (
        "locus_tag" not in feature.qualifiers or args.use_protein_ids) else "locus_tag"
    logger.debug("Using %s as ID", feature.qualifiers[field][0])
    ntrecord = SeqRecord(
        seq=ntseq,
        description="coding sequence",
        id=feature.qualifiers[field][0],
    )
    aarecord = SeqRecord(
        seq=aaseq,
        description="conceptual translation",
        id=feature.qualifiers[field][0],
    )
    # if (args.use_protein_ids) or ("locus_tag" not in feature.qualifiers):
    #     ntrecord = SeqRecord(
    #         seq=ntseq,
    #         description="coding sequence",
    #         id=feature.qualifiers["protein_id"][0],
    #     )
    #     aarecord = SeqRecord(
    #         seq=aaseq,
    #         description="conceptual translation",
    #         id=feature.qualifiers["protein_id"][0],
    #     )
    # else:
    #     ntrecord = SeqRecord(
    #         seq=ntseq,
    #         description="coding sequence",
    #         id=feature.qualifiers["locus_tag"][0],
    #     )
    #     aarecord = SeqRecord(
    #         seq=aaseq,
    #         description="conceptual translation",
    #         id=feature.qualifiers["locus_tag"][0],
    #     )

    return ntrecord, aarecord


def strip_stockholm_from_seqid(seqid):
    """Strip Stockholm header from seqid.

    seqid        - seqid to strip
    """
    return seqid.split("/")[0]
