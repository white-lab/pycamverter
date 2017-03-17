"""
Main module for running pycamv from the commandline.
"""

from __future__ import absolute_import, division

import multiprocessing
import sys

from . import main

if __name__ == "__main__":
    multiprocessing.freeze_support()

    main.main(sys.argv[1:])
