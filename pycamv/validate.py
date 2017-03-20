"""
This module provides functionality for interacting with CAMV.

Currently limited to importing and outputing scan lists.
"""

# Built-ins
from __future__ import absolute_import, division

from collections import OrderedDict
import logging
import multiprocessing
import os
import tempfile
import shutil

from . import compare, export, fragments, gen_sequences, scans, search
from .scan_list import load_scan_list
from .utils import LenGen


LOGGER = logging.getLogger("pycamv.validate")


def _remap_pst(pep_mods):
    return [
        (
            count,
            mod,
            letters + (
                ("Y",)
                if (mod == "Phospho" and set(letters) == set(["S", "T"])) else
                ()
            ),
        )
        for count, mod, letters in pep_mods
    ]


def _map_seq(kv):
    pep_query = kv
    return (
        pep_query,
        tuple(
            seq
            for index, seq in zip(
                range(10),
                gen_sequences.gen_possible_seq(
                    pep_query.pep_seq,
                    pep_query.pep_mods,
                )
            )
        ),
    )


def _map_frag_compare(kv):
    pep_query, scan_query, sequence, scan = kv

    frag_ions = fragments.fragment_ions(
        sequence, pep_query.pep_exp_z,
        c13_num=scan_query.c13_num,
    )

    LOGGER.debug("{} - {} ions".format(pep_query.pep_seq, len(frag_ions)))

    peaks = compare.compare_spectra(
        scan, frag_ions,
        tol=compare.COLLISION_TOLS[scan_query.collision_type],
    )

    del scan

    return (pep_query, tuple(sequence)), peaks


def validate_spectra(
    search_path,
    raw_paths,
    scans_path=None,
    scan_list=None,
    out_path=None,
    cpu_count=None,
):
    """
    Generate CAMV web page for validating spectra.

    Parameters
    ----------
    search_path : str
    raw_paths : list of str
    scans_path : str, optional
    scan_list : list of int, optional
    out_path : str, optional
    cpu_count : int, optional
    """
    if cpu_count is None:
        try:
            cpu_count = multiprocessing.cpu_count()
        except NotImplementedError:
            cpu_count = 2

    LOGGER.debug("Validating data set with {} cpus".format(cpu_count))

    if scan_list is None:
        scan_list = []

    if scans_path is not None:
        scan_list += load_scan_list(scans_path)

    if out_path is None:
        out_path = os.path.splitext(search_path)[0] + ".camv.gz"

    # Read peptide search file
    fixed_mods, var_mods, pep_queries = search.read_search_file(search_path)

    LOGGER.info("Found info for {} peptide queries".format(len(pep_queries)))

    # Optionally filter queries using a scan list
    if scan_list:
        pep_queries = [
            query
            for query in pep_queries
            if query.scan in scan_list
        ]

    # Remap pST -> pSTY
    LOGGER.info("Remapping pST -> pSTY")
    for pep_query in pep_queries:
        pep_query.pep_var_mods = _remap_pst(pep_query.pep_var_mods)
        pep_query.pep_fixed_mods = _remap_pst(pep_query.pep_fixed_mods)

    # Get scan data from RAW file
    required_raws = set(query.basename for query in pep_queries)
    base_raw_paths = [os.path.basename(path) for path in raw_paths]

    for base_raw in required_raws:
        if base_raw not in base_raw_paths:
            raise Exception(
                "Unable to find {} in input RAW files: {}"
                .format(base_raw, base_raw_paths)
            )

    # Generate sequences
    LOGGER.info(
        "Generating all possible sequence-modification combinations."
    )
    pool = multiprocessing.Pool(processes=cpu_count)
    sequence_mapping = dict(
        pool.imap_unordered(_map_seq, pep_queries)
    )

    total_num_seq = sum(len(i) for i in sequence_mapping.values())

    LOGGER.info(
        "Generating fragment ions for {} possible sequences."
        .format(
            total_num_seq,
        )
    )

    out_dir = tempfile.mkdtemp()

    LOGGER.info("Getting scan data.")

    scan_queries, ms_two_data, ms_data = scans.get_scan_data(
        raw_paths, pep_queries, out_dir,
    )

    LOGGER.info("Found data for {} scans".format(len(scan_queries)))

    scan_mapping = OrderedDict(
        (pep_query, scan_query)
        for pep_query, scan_query in zip(pep_queries, scan_queries)
    )

    LOGGER.info(
        "Comparing predicted to actual peaks for {} spectra.".format(
            len(scan_mapping),
        )
    )

    peak_hits = dict(
        pool.imap_unordered(
            _map_frag_compare,
            LenGen(
                (
                    (
                        pep_query,
                        scan_mapping[pep_query],
                        sequence,
                        ms_two_data[pep_query.basename][pep_query.scan]
                        .deRef(),
                    )
                    for pep_query, sequences in sequence_mapping.items()
                    for sequence in sequences
                ),
                total_num_seq,
            ),
            cpu_count * 4,
        )
    )

    del pool

    # XXX: Determine SILAC precursor masses?

    # XXX: Remove precursor contaminated scans from validation list?

    # Check each assignment to each scan

    # Output data
    export.export_to_camv(
        out_path,
        peak_hits, scan_mapping,
        ms_data, ms_two_data,
    )

    LOGGER.info("Removing directory of temporary files.")
    del ms_data
    del ms_two_data
    shutil.rmtree(out_dir)
