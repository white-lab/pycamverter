
"""
Main module for running pycamv from the commandline.
"""

import sys

from . import main

if __name__ == "__main__":
    main.main(sys.args[1:])
