#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions to interact with Entrez

(c) The James Hutton Institute 2017
Author: Leighton Pritchard

Contact: leighton.pritchard@hutton.ac.uk
Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD6 9LH,
Scotland,
UK

The MIT License

Copyright (c) 2017 The James Hutton Institute

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
"""

from Bio import Entrez


# EXCEPTIONS
# ==========
# General exception for scripts
class NCFPELinkException(Exception):
    """Exception for ELink qeries."""

    def __init__(self, msg="Error in ncfp ELink query"):
        """Instantiate class."""
        Exception.__init__(self, msg)



def elink_fetch_with_retries(record, dbfrom, linkname, retries):
    """Returns ELink search result.

    Allows the maximum number of retries specified.
    """
    tries = 0
    while tries < retries:
        try:
            link = Entrez.elink(dbfrom=dbfrom,
                                linkname=linkname,
                                id=record.query)
            matches = Entrez.read(link)
            return matches
        except:
            tries += 1
    raise NCFPELinkException("Could not complete ELink query")
