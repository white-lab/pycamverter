"""
Provides functionality for interacting with MASCOT data.
"""

from __future__ import absolute_import, division

import logging
import os

from . import discoverer, mascot


LOGGER = logging.getLogger("pycamv.search")

BACKENDS = {
    ".xml": mascot.read_mascot_xml,
    ".msf": discoverer.read_discoverer_msf,
}


def read_search_file(path):
    """
    Parse a search input file.

    Parameters
    ----------
    path : str
        Path to search input file.

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

    return backend(path)
