# try using distribute or setuptools or distutils.
try:
    import distribute_setup
    distribute_setup.use_setuptools()
except ImportError:
    pass

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import os
import sys
import re

# parse version from package/module without importing or evaluating the code
with open('ncbi_cds_from_protein/__init__.py') as fh:
    for line in fh:
        m = re.search(r"^__version__ = '(?P<version>[^']+)'$", line)
        if m:
            version = m.group('version')
            break

if sys.version_info <= (3, 0):
    sys.stderr.write("ERROR: pyani requires Python 3 " +
                     "or above...exiting.\n")
    sys.exit(1)

setup(
    name="ncbi_cds_from_protein",
    version=version,
    author="Leighton Pritchard",
    author_email="leighton.pritchard@hutton.ac.uk",
    description=' '.join(["ncfp is a script and module that facilitates",
                          "recovery of nucleotide sequecnes from NCBI",
                          "given a set of input protein sequences"]),
    license="MIT",
    keywords="genome bioinformatics sequence",
    platforms="Posix; MacOS X",
    url="http://widdowquinn.github.io/ncbi_cds_from_protein/",  # project home page
    download_url="https://github.com/widdowquinn/ncbi_cds_from_protein/releases",
    scripts=[os.path.join('bin', 'ncfp')],
    packages=['ncbi_cds_from_protein'],
    package_data={},
    include_package_date=True,
    install_requires=['biopython', ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
