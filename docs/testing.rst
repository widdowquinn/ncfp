.. _ncfp-testing:

=======
Testing
=======

``ncfp`` tests use the `pytest`_ framework [#f1]_.

To run all tests locally, please issue the following command from the root of the repository:

.. code-block:: bash

    pytest -v

-----------------------------------
Obtaining test coverage information
-----------------------------------

The `pytest`_ framework integrates with the `coverage.py`_ module, and an account
of the extent of test coverage can be obtained by running the following command:

.. code-block:: bash

    pytest -v --with-coverage \
        --cover-package=ncbi_cds_from_protein \
        --cover-html

----------------------
Continuous Integration
----------------------

Continuous Integration for ``ncfp`` is provided by `CircleCI`_.

.. _CircleCI: https://app.circleci.com/pipelines/github/widdowquinn/ncfp
.. _coverage.py: https://coverage.readthedocs.io/en/coverage-4.5/
.. _pytest: https://docs.pytest.org/en/latest/

.. rubric:: Footnotes
