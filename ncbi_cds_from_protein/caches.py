#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions for handling data caches

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

import os

from collections import namedtuple


# Paths cache files for script
Cachepaths = namedtuple("Cachepaths", "elink acc gb gbfull")


# Initialise caches
def initialise_caches(cachedir, cachestem):
    """Initialise caches for downloading.

    cachedir     - path to cache directory
    cachestem    - suffix for caches to identify a run
    """
    cachepaths = Cachepaths(os.path.join(cachedir,
                                         'elink_{}'.format(cachestem)),
                            os.path.join(cachedir,
                                         'gb_{}'.format(cachestem)),
                            os.path.join(cachedir,
                                         'gbfull_{}'.format(cachestem)),
                            os.path.join(cachedir,
                                         'acc_{}'.format(cachestem)))
    for path in cachepaths:
        check_and_create_cache(path)
    return cachepaths


# Create cache if it doesn't exist
def check_and_create_cache(path):
    """If no file at passed path, creates one."""
    if os.path.isfile(path):
        return
    os.makedirs(os.path.split(path)[0], exist_ok=True)
    with open(path, 'w') as cfh:
        cfh.write("#ncfp cache: %s" % path)
