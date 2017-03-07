"""
Main module for running pycamv from the commandline.
"""

from __future__ import absolute_import, division

import sys

from . import main

if __name__ == "__main__":
    main.main(sys.args[1:])
