.. _ncfp-testing:

=======
Testing
=======

``ncfp`` tests are implemented in the `Nose`_ framework [#f1]_.

To run all tests locally, please issue the command:

.. code-block:: bash

    nosetests -v

-----------------------------------
Obtaining test coverage information
-----------------------------------

The `Nose`_ framework integrates with the `coverage.py`_ module, and an account
of the extent of test coverage can be obtained by running the following command:

.. code-block:: bash

    nosetests -v --with-coverage \
        --cover-package=ncbi_cds_from_protein \
        --cover-html



.. _coverage.py: https://coverage.readthedocs.io/en/coverage-4.5/
.. _Nose: https://nose.readthedocs.io/en/latest/

.. rubric:: Footnotes

.. [#f1] We are aware that ``nosetests`` is in maintenance mode, and a move to ``py.test`` is planned.