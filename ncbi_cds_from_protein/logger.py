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
"""Code providing a script logger"""

import logging
import logging.config
import re
import sys

from argparse import Namespace
from pathlib import Path


class NoColorFormatter(logging.Formatter):

    """Log formatter that strips terminal colour escape codes from the log message."""

    ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")

    def format(self, record):
        """Return logger message with terminal escapes removed."""
        return "[%s] [%s]: %s" % (
            record.levelname,
            record.name,
            re.sub(self.ANSI_RE, "", record.msg % record.args),
        )


def config_logger(args: Namespace):
    """Configure package/script-level logging

    :param args:  CLI arguments

    Instantiates a logger for the script, and adds basic info.
    """
    # Default package-level logger
    logger = logging.getLogger(__package__)
    logger.setLevel(logging.WARNING)

    # Create and add STDERR handler
    err_formatter = logging.Formatter("[%(levelname)s] [%(name)s]: %(message)s")
    err_handler = logging.StreamHandler(sys.stderr)
    if args is not None and args.verbose:
        err_handler.setLevel(logging.INFO)
    elif args is not None and args.debug:
        err_handler.setLevel(logging.DEBUG)
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)

    # If args.logfile is provided, add a FileHandler for it
    if args.logfile is not None:
        logdir = args.logfile.parents[0]
        try:
            if not logdir == Path.cwd():
                logdir.mkdir(exist_ok=True, parents=True)
        except OSError:
            logger.error(
                "Could not create log directory %s (exiting)", logdir, exc_info=True
            )
            raise SystemExit(1)

        # create handler
        logformatter = NoColorFormatter()
        loghandler = logging.FileHandler(args.logfile, mode="w", encoding="utf8")
        if args.debug:
            loghandler.setLevel(logging.DEBUG)
        else:
            loghandler.setLevel(logging.INFO)
        loghandler.setFormatter(logformatter)
        logger.addHandler(loghandler)
