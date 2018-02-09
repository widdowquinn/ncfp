.. _ncfp-examples:


========
Examples
========

----------------------------------------
NCBI format input sequences - no introns
----------------------------------------

The command below [#f1]_ identifies coding sequences from NCBI format
[#f2]_ input for nine proteins that do not contain introns. 

.. code-block:: bash

    ncfp tests/test_input/sequences/input_ncbi.fasta \
        tests/test_output/ncbi dev@null.com -v


-------------------------------------------
UniProt format input sequences - no introns
-------------------------------------------

The command below [#f1]_ identifies coding sequences from UniProt
format [#f3]_ input for ten proteins that do not contain introns. The
``-u`` or ``--uniprot`` argument is required to specify that the input
sequences are UniProt format, otherwise an error is thrown.

.. code-block:: bash

    ncfp -u tests/test_input/sequences/input_uniprot.fasta \
        tests/test_output/uniprot dev@null.com -v


----------------------------------------------
UniProt/Stockholm input sequences - no introns
----------------------------------------------

The command below [#f1]_ identifies coding sequences from UniProt
format [#f3]_ input for 57 amino acid sequences specifying regions
of a protein in Stockholm notation [#f4]_. The ``-u`` or ``--uniprot``
argument is required to specify that the input sequences are UniProt
format, and the ``-s`` or ``--stockholm`` arguments are required to
tell ``ncfp`` to parse the region locations.

.. code-block:: bash

    ncfp -us tests/test_input/sequences/input_uniprot_stockholm.fasta \
        tests/test_output/uniprot_stockholm dev@null.com -v


----------------------------------------------------
Human sequences - isoforms and intron/exon structure
----------------------------------------------------

The command below [#f1]_ identifies coding sequences from NCBI
format [#f3]_ input for four human proteins with intron/exon structure,
including three isoforms of the same protein from the same locus
(GPR137: ``NP_001164351.1``, ``NP_001164352.1``, and ``XP_005274161.1``).

.. code-block:: bash

    ncfp tests/test_input/sequences/human.fasta \
        tests/test_output/human dev@null.com -v

-------
Logging
-------

Verbose output can be written persistently to a logfile using the
``-l`` or ``--logfile`` argument and specifying the path to which
the logfile should be written. An example is given in the command below.

.. code-block:: bash

    ncfp tests/test_input/sequences/human.fasta \
        tests/test_output/logging dev@null.com \
        -l tests/test_output/logging/human.log


-----------------------------
Specifying the cache location
-----------------------------

By default a new cache database is created every time that ``ncfp`` is
run, in the ``.ncfp_cache`` hidden subdirectory. The default cache
database filename is ``ncfpcache_YYYY-MM-DD-HH-MM-SS.sqlite3``,
indicating the time that the command was run. This location and
naming convention can be overridden with the ``-d``/``--cachedir`` and
``-c``/``--cachestem`` arguments, as in the command below.

.. code-block:: bash

    ncfp tests/test_input/sequences/human.fasta \
        tests/test_output/caches dev@null.com \
        -d tests/test_output/caches \
        -c ncfp_cache

-------------------------
Reusing an existing cache
-------------------------

To avoid unnecessary bandwidth/``NCBI`` queries, an existing cache
database can be used. The location of the cache is specified with the
``-d``/``--cachedir`` and ``-c``/``--cachestem`` arguments, and the
``--keepcache`` option must be specified. If the specified location
does not contain a cache database, one is created. For example:

.. code-block:: bash

    ncfp tests/test_input/sequences/human.fasta \
        tests/test_output/caches dev@null.com \
        -d tests/test_output/caches \
        -c ncfp_cache

will create a cache at ``tests/test_output/caches/ncfp_cache.sqlite3``,
and

.. code-block:: bash

    ncfp tests/test_input/sequences/human.fasta \
        tests/test_output/caches dev@null.com \
        -d tests/test_output/caches \
        -c ncfp_cache \
        --filestem cached \
        --keepcache

will reuse the cachefile without making new queries at ``NCBI``, and
write the output to ``cached_aa.fasta`` and ``cached_nt.fasta`` [#f5]_.



.. rubric:: Footnotes

.. [#f1] The ``-v`` option shows verbose output in ``STDOUT``.
.. [#f2] The sequence identifier in the FASTA header is a valid NCBI protein accession.
.. [#f3] The sequence description in the FASTA header contains a valid ``GN=<accession>`` gene identifier.
.. [#f4] The sequence identifier in the FASTA header ends with a Stockholm format region definition, e.g. ``/47-134``.
.. [#f5] The ``--filestem`` argument changes the filestem of the output nucleotide and amino acid sequence files.
