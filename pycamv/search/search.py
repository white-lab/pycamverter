"""
Provides functionality for interacting with MASCOT data.
"""

from __future__ import absolute_import, division

import logging
import os

from . import discoverer

LOGGER = logging.getLogger("pycamv.search")

BACKENDS = {
    ".msf": discoverer.read_discoverer_msf,
}


def read_search_file(path, scan_list=None, score=None):
    """
    Parse a search input file.

    Parameters
    ----------
    path : str
        Path to search input file.
    scan_list : list of int, optional
    score : int, optional

    Returns
    -------
    fixed_mods : list of str
    var_mods : list of str
    out : list of :class:`PeptideQuery<pycamv.pep_query.PeptideQuery>`
    """
    ext = os.path.splitext(path)[1]
    backend = BACKENDS.get(ext)

    if backend is None:
        raise Exception("Unable to open search file: {}".format(path))

    LOGGER.debug("Using {} backend for {}".format(backend.__name__, ext))

    return backend(path, scan_list=scan_list, score=score)
