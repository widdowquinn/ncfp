.. _ncfp-quickstart:

================
QuickStart Guide
================

------------
Installation
------------

At the command-line, use ``git`` to clone the current version of the ``ncfp`` repository:

.. code-block:: bash

    git clone git@github.com:widdowquinn/ncfp.git

Change to the newly-created ``ncfp`` subdirectory:

.. code-block:: bash

    cd ncfp

Install the package and program, using the ``setup.py`` file:

.. code-block:: bash

    python setup.py install

(other installation methods can be found at :ref:`ncfp-installation`)

------------------------
``ncfp`` Program Example
------------------------

To see the options available for the ``ncfp`` program, use the ``-h``
(help) option:

.. code-block:: bash

    ncfp -h

In the ``ncfp/tests/test_input/sequences`` subdirectory there is a file
called ``input_ncbi.fasta``. This contains a number of protein sequences in
FASTA format. The identifier for each sequence in this file is a valid NCBI
sequence identifier.

Using ``ncfp``, to obtain a corresponding nucleotide coding sequence for
each protein, issue the following command (substituting your own email
address, where indicated):


.. code-block:: bash

    ncfp tests/test_input/sequences/input_ncbi.fasta \
         example_output \
         my.name@my.domain

You should see progress bars appear for processing of the input protein,
sequences, searching those sequences against the remote ``NCBI`` databases,
then retrieving the corresponding sequence identifiers, GenBank headers and
finally the full GenBank records.

On completion, a list of the recovered sequences will be presented,
and the directory ``example_output`` will be created, with the following
contents:

.. code-block:: bash

    $ tree example_output/
    example_output/
    ├── ncfp_aa.fasta
    └── ncfp_nt.fasta
