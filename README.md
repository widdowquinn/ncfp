# README.md - `ncfp`

This repository contains code for a script that identifies and writes the corresponding nucleotide sequences for each protein in an input multiple sequence file to be used, for example, in backthreading coding sequences onto protein alignments for phylogenetic analyses. `ncfp` uses the NCBI accession or UniProt gene name (as appropriate) to identify source nucleotide sequences in the NCBI databases, download them, and write them to a file.

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/53ecf92293ae4a97820fde81d1bd947e)](https://app.codacy.com/app/widdowquinn/ncfp?utm_source=github.com&utm_medium=referral&utm_content=widdowquinn/ncfp&utm_campaign=Badge_Grade_Dashboard)
[![ncfp TravisCI build status](https://api.travis-ci.org/widdowquinn/ncfp.svg?branch=master)](https://travis-ci.org/widdowquinn/ncfp/branches)
[![ncfp codecov.io coverage](https://img.shields.io/codecov/c/github/widdowquinn/ncfp/master.svg)](https://codecov.io/github/widdowquinn/ncfp)
[![ncfp documentation](https://readthedocs.org/projects/ncfp/badge/?version=latest)](https://ncfp.readthedocs.io/en/latest/?badge=latest)

[![ncfp PyPi version](https://img.shields.io/pypi/v/ncfp.svg "PyPi version")](https://pypi.python.org/pypi/ncfp)
[![ncfp licence](https://img.shields.io/pypi/l/ncfp.svg "PyPi licence")](https://github.com/widdowquinn/ncfp/blob/master/LICENSE)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/ncfp/badges/version.svg)](https://anaconda.org/bioconda/ncfp)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/ncfp/badges/latest_release_date.svg)](https://anaconda.org/bioconda/ncfp)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/ncfp/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/ncfp/badges/downloads.svg)](https://anaconda.org/bioconda/ncfp)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/ncfp/badges/platforms.svg)](https://anaconda.org/bioconda/ncfp)


## Quickstart: `ncfp` at the command-line

Providing an input file of protein sequences as `<INPUT>.fasta`, and writing output to the directory `<OUTPUT>`, while specifying a user email to NCBI of `<EMAIL>` will generate two files: `<OUTPUT>/ncfp_aa.fasta` and `<OUTPUT>/ncfp_nt.fasta`.

```bash
ncfp <INPUT>.fasta <OUTPUT> <EMAIL>
```

The file `<OUTPUT>/ncfp_aa.fasta` contains the amino acid sequences for all input proteins for which a corresponding nucleotide coding sequence could be identified, in FASTA format.

The file `<OUTPUT>/ncfp_nt.fasta` contains nucleotide coding sequences, where they could be found, for all the input proteins, in FASTA format.

Any input protein sequences for which a corresponding nucleotide sequence cannot be recovered, for any reason, are placed in the file `<OUTPUT>/skipped.fas`.

To find out more about what `ncfp` can do, try

```bash
ncfp --help
```

at the command-line

## Documentation

For more detailed information about `ncfp` as a program, or using the underlying `ncbi_cds_from_protein` Python module, please see the most recent documentation at <https://ncfp.readthedocs.io/en/latest/>

## License

Unless otherwise indicated, all code is licensed under the MIT license and subject to the following agreement:

    (c) The James Hutton Institute 2017-2018
    Author: Leighton Pritchard

    Contact: leighton.pritchard@hutton.ac.uk

    Address:
    Leighton Pritchard,
    Information and Computational Sciences,
    James Hutton Institute,
    Errol Road,
    Invergowrie,
    Dundee,
    DD2 5DA,
    Scotland,
    UK

The MIT License

Copyright (c) 2017-2018 The James Hutton Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
