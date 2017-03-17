"""
Main module for running pycamv from the commandline.
"""

from __future__ import absolute_import, division

from . import multi

import multiprocessing
import sys

from . import main

if __name__ == "__main__":
    multiprocessing.freeze_support()

    main.main(sys.args[1:])
