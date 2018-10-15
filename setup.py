# try using distribute or setuptools or distutils.
try:
    import distribute_setup

    distribute_setup.use_setuptools()
except ImportError:
    pass

import setuptools

import os
import sys
import re

# Get long description from README.md
with open("README.md", "r") as dfh:
    long_description = dfh.read()

# parse version from package/module without importing or evaluating the code
vfname = os.path.join("ncbi_cds_from_protein", "__init__.py")
with open(vfname, "r") as fh:
    for line in fh:
        m = re.search(r"^__version__ = '(?P<version>[^']+)'$", line)
        if m:
            version = m.group("version")
            break

if sys.version_info < (3, 5):
    sys.stderr.write("ERROR: ncfp requires Python 3.5 " + "or above...exiting.\n")
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
    install_requires=["biopython", "tqdm"],
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
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
