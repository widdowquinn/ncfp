# -*- coding: utf-8 -*-
# (c) The University of Strathclude 2019-2020
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute of Pharmaceutical and Biomedical Sciences
# The University of Strathclyde
# 161 Cathedral Street
# Glasgow
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# (c) The University of Strathclude 2019-2020
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
"""Support for testing"""

import copy

from argparse import Namespace


def check_files(outdirname, tgtdirname, fnames):
    """Test whether output and target files are identical in output and target directories"""
    for fname in fnames:
        with (outdirname / fname).open() as fh1:
            with (tgtdirname / fname).open() as fh2:
                assert fh1.read() == fh2.read()


def modify_namespace(namespace: Namespace, **kwargs):
    """Update arguments in a passed Namespace.

    :param namespace:       argparse.Namespace object
    :param kwargs:          keyworded arguments

    The expected usage pattern is, for a command-line application with many
    or complex arguments, to define a base argparse.Namespace object, then
    change only a few arguments, specific to a test. This function takes
    a base namespace and a dictionary of argument: value pairs, and
    returns the modified namespace.
    """
    new_namespace = copy.deepcopy(namespace)
    for argname, argval in kwargs.items():
        setattr(new_namespace, argname, argval)

    return new_namespace
