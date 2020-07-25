#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) University of Strathclyde 2019-2020
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
# Copyright (c) 2019-2020 University of Strathclyde
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
"""Distribution setup."""

import sys
import re
import setuptools

from pathlib import Path

# try using distribute or setuptools or distutils.
try:
    import distribute_setup

    distribute_setup.use_setuptools()
except ImportError:
    pass

# Get long description from README.md
with Path("README.md").open("r") as dfh:
    long_description = dfh.read()  # pyline: disable=C0103

# parse version from package/module without importing or evaluating the code
with Path("ncbi_cds_from_protein/__init__.py").open("r") as ifh:
    for line in ifh:
        print(line)
        match = re.search(r'^__version__ = "(?P<version>[^"]+)"$', line)
        if match:
            version = match.group("version")
            break

if sys.version_info <= (3, 6):
    sys.stderr.write("ERROR: ncfp requires Python 3.7 or above (exiting).\n")
    sys.exit(1)

setuptools.setup(
    name="ncfp",
    version=version,
    author="Leighton Pritchard",
    author_email="leighton.pritchard@hutton.ac.uk",
    description=" ".join(
        [
            "ncfp is a script and module that facilitates",
            "recovery of nucleotide sequences from NCBI",
            "encoding a specified set of input protein sequences",
        ]
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    keywords="genome bioinformatics sequence",
    platforms="Posix; MacOS X",
    url="http://widdowquinn.github.io/ncfp/",  # project home page
    download_url="https://github.com/widdowquinn/ncfp/releases",
    packages=["ncbi_cds_from_protein", "ncbi_cds_from_protein.scripts"],
    package_data={},
    include_package_date=True,
    install_requires=["biopython", "bioservices", "tqdm"],
    python_requires="~=3.5",
    entry_points={
        "console_scripts": ["ncfp = ncbi_cds_from_protein.scripts.ncfp:run_main"]
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
