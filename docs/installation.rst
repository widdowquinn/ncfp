.. _ncfp-installation:

============
Installation
============

``ncfp`` can be installed in any of several ways.

-------------
Using ``pip``
-------------

The most recent release of ``ncfp`` is available at the `PyPI`_ warehouse, and can be installed using ``pip``:

.. code-block:: bash

    pip install ncfp


------------------
Using ``bioconda``
------------------

``ncfp`` is available through the `bioconda`_ channel of the `conda`_ package management system. To install
``ncfp``, you will need `Anaconda`_ or `miniconda`_, and to set up the ``bioconda`` channel. Then you can use
the ``conda install ncfp`` command to install the package:

.. code-block:: bash

    conda install -c bioconda ncfp

Alternatively:

.. code-block:: bash

    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda install ncfp


---------------------------------
From source (most recent release)
---------------------------------

``ncfp`` releases are available at the `GitHub releases page`_. To install from source:

1. Download the `.tar.gz` or `.zip` file containing the package (click the link, use ``curl`` or ``wget``).
2. Uncompress the archive.
3. Change to the newly-expanded directory.
4. Install Python module requirements with ``pip``.
5. Install ``ncfp`` with ``setup.py``

.. code-block:: bash

    $ wget https://github.com/widdowquinn/ncfp/archive/v0.2.0.tar.gz
    $ tar -zxvf v0.2.0.tar.gz
    $ cd ncfp
    $ pip install -r requirements.txt
    $ python setup.py install

---------------------------
From source (bleeding edge)
---------------------------

To get the very latest development version of ``ncfp``, you can clone the repository from the `GitHub project page`_

1. Clone the repository from ``git@github.com:widdowquinn/ncfp.git``
2. Change to the repository root directory.
3. Install Python module requirements with ``pip``.
4. Install ``ncfp`` with ``setup.py``

.. code-block:: bash

    $ git clone git@github.com:widdowquinn/ncfp.git
    $ cd ncfp
    $ pip install -r requirements.txt
    $ pip install -e .  # for development installation
    $ python setup.py install  # for installation as a package/script


.. _conda: https://conda.io/
.. _bioconda: https://bioconda.github.io/
.. _miniconda: https://conda.io/miniconda.html
.. _Anaconda: https://www.anaconda.com/
.. _GitHub project page: https://github.com/widdowquinn/ncfp
.. _GitHub releases page: https://github.com/widdowquinn/ncfp/releases
.. _PyPI: https://pypi.python.org/pypi/ncfp