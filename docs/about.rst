.. _ncfp-about:

==============
About ``ncfp``
==============

``ncfp`` is a program and accompanying Python package (``ncbi_cds_from_protein``) that, given a set of
protein sequences with appropriate identifier values, retrieves corresponding coding nucleotide
sequences from the ``NCBI nucleotide`` databases. This is useful, for instance, to help backthread coding
sequences onto protein alignments for phylogenetic analyses.

------------------
How ``ncfp`` works
------------------

1. The ``ncfp`` program accepts input FASTA sequence files describing protein sequences. These describe a set
of queries that will be used to obtain corresponding nucleotide coding sequences, if possible.

Input sequence formats
    Please see :ref:`input-sequence-formats`.

2. An `SQLite`_ cache database is created. This will hold information about the query sequences, and about
the data downloaded from ``NCBI`` using the input sequence data as queries. The cache enables recovery from
interrupted jobs, and reuse of data without needing to transfer across the network again by interrogating
the SQLite database directly.

Cache directory
    By default, the cache is created in the hidden subdirectory ``.ncfp_cache``, but any location can be specified
    using the ``-d`` or ``--cachedir`` arguments.

Cache filename
    By default, the cache has a filestem reflecting the date and time that ``ncfp`` is run, but 
    this can be changed using the ``-c`` or ``--cachestem`` arguments.

3. The `Biopython`_ ``Entrez`` library is used to make a connection to the ``NCBI`` sequence
databases. Using this connection, the program identifies ``nucleotide`` database coding sequence entries
that are related to each input protein sequence. The relationship is determined on the basis of either
the NCBI accession, or the UniProt gene identifier, depending on the input file.

Batched downloads
    By default, ``ncfp`` makes queries and downloads sequences in batches of 100. The batch size can be
    controlled using the arguments ``-b`` or ``--batchsize``

Retry attempts
    Sometimes network connections are flaky, so by default ``ncfp`` will try each request up to 10 times. The
    number of retries can be controlled with the ``-r`` or ``--retries`` arguments.

4. If the results of each ``NCBI`` query are not already present in the cache, they are downloaded and
recorded in the cache as header information. Some specific data are extracted, (sequence length, taxonomy, etc.)

5. The shortest available complete coding sequence [#f1]_ that recapitulates each input protein sequence
is identified. If the sequence is not already present in the cache, it is downloaded.

6. Pairs of protein and corresponding coding sequences are written to two files: one for nucleotide sequences,
and one for protein sequences. Sequences are written to each file in teh same order, so they can be used for
backtranslation with a tool such as `T-Coffee`_. If any proteins could not be matched to their coding
sequence at ``NCBI``, they are written to a third file.

Output filenames
    The filestem for the paired protein and coding sequence files is always suffixed by ``_aa`` or ``_nt``,
    depending on the type of sequence being written. The filestem is ``ncfp`` by default, but this can be
    controlled with the ``--filestem`` argument.

Skipped sequences filename
    By default, the protein sequences for which a coding sequence could not be found are written to the
    ``skipped.fas`` file. An alternative path can be provided with the ``--skippedfile`` argument.



.. _Biopython: http://biopython.org/
.. _SQLite: https://www.sqlite.org/
.. _T-Coffee: http://www.tcoffee.org/Projects/tcoffee/

.. rubric:: Footnotes

.. [#f1] We require the complete coding sequence, but if we can use a shorter sequence than the complete genome, we do to save bandwidth.