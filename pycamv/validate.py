"""
This module provides functionality for interacting with CAMV.

Currently limited to importing and outputing scan lists.
"""

# Built-ins
from __future__ import absolute_import, division

from collections import OrderedDict
from functools import partial
import logging
import multiprocessing
import os
import tempfile
import shutil

from . import (
    camv_mat, compare, export, fragments, gen_sequences, ms_labels, scans,
    search,
)
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


def _map_seq(kv, limit_comb=False):
    lowpriority()

    pep_query = kv

    gen = gen_sequences.gen_possible_seq(
        pep_query.pep_seq,
        pep_query.pep_mods,
    )

    if limit_comb:
        return (
            pep_query,
            tuple(
                seq
                for index, seq in zip(
                    range(gen_sequences.MAX_NUM_COMB),
                    gen,
                )
            )
        )
    else:
        return (
            pep_query,
            tuple(gen),
        )


def _to_str(seq):
    return "".join(
        letter.lower()
        if [i for i in mods if i not in ms_labels.LABEL_NAMES]
        else letter.upper()
        for letter, mods in seq[1:-1]
    )


def _close_scans(ms_datas, out_dir=None):
    LOGGER.info("Removing directory of temporary files.")

    for ms_data in ms_datas:
        for raw in ms_data.values():
            raw.info['fileObject'].close()
            raw.seeker.close()

    if out_dir:
        shutil.rmtree(out_dir)


def _get_window_coverage(pep_query, scan_query, precursor_win):
    precursor_win = [
        mz
        for mz, _ in precursor_win
        if mz >= scan_query.isolation_mz - scan_query.window_offset[0] and
        mz <= scan_query.isolation_mz + scan_query.window_offset[1]
    ]
    return (
        max([
            c13
            for c13 in range(
                1, 1 + round(scan_query.window_offset[1] * pep_query.pep_exp_z)
            )
            if min([
                1e6 * abs(
                    scan_query.isolation_mz +
                    fragments.DELTA_C13 * c13 / pep_query.pep_exp_z -
                    mz
                ) / mz
                for mz in precursor_win
            ] + [compare.MS_TOL]) < compare.MS_TOL
        ] + [0])
        if scan_query.window_offset
        else 0
    )


def _map_frag_compare(kv):
    lowpriority()

    try:
        (
            pep_query, scan_query, sequence,
            ms_two_data, ms_data,
            validation_data, auto_maybe,
        ) = kv

        ms_scan = ms_data[scan_query.basename][scan_query.precursor_scan]

        if ms_scan["id"] != scan_query.precursor_scan:
            LOGGER.warning(
                "{} #{}: Precursor scan id different from expected: {} != {}"
                .format(
                    pep_query.query,
                    pep_query.scan,
                    ms_scan["id"],
                    scan_query.precursor_scan,
                )
            )

        precursor_win = scans.get_precursor_peak_window(
            scan_query, ms_scan
        )

        # Get the max C13 peak found within the isolation window
        window_coverage = _get_window_coverage(
            pep_query, scan_query, precursor_win,
        )

        frag_ions = fragments.fragment_ions(
            sequence, pep_query.pep_exp_z,
            c13_num=scan_query.c13_num + window_coverage,
        )

        # Compare MS^2 data with predicted fragment ions
        ms_two_scan = ms_two_data[pep_query.basename][pep_query.scan]

        if ms_two_scan["id"] != scan_query.scan:
            LOGGER.warning(
                "{} #{}: MS^2 scan id different from expected: {} != {}"
                .format(
                    pep_query.query,
                    pep_query.scan,
                    ms_two_scan["id"],
                    scan_query.scan,
                )
            )

        peaks = compare.compare_spectra(
            ms_two_scan, frag_ions,
            tol=compare.COLLISION_TOLS[scan_query.collision_type],
        )

        # XXX: Unlabeled runs?
        quant_scan = (
            ms_two_data[pep_query.basename][pep_query.quant_scan]
            if pep_query.scan != pep_query.quant_scan else
            ms_two_scan
        ) if pep_query.quant_scan is not None else None

        label_win = scans.get_label_peak_window(pep_query, quant_scan)

        choice = validation_data.get((pep_query.scan, _to_str(sequence)), None)

        if not choice and auto_maybe:
            if (
                pep_query.rank_pos is not None and
                pep_query.rank_pos.get(1, None) == set(
                    (pos, mod)
                    for pos, (_, mods) in enumerate(sequence[1:-1])
                    for mod in mods
                )
            ):
                choice = "maybe"

        return (
            pep_query, tuple(sequence), choice,
            peaks, precursor_win, label_win,
        )
    except:
        _close_scans([ms_data, ms_two_data])
        raise


def fill_map_frag_compare(
    sequence_mapping, scan_mapping,
    ms_two_data, ms_data,
    queue, cpu_count,
    validation_data, auto_maybe,
):
    pool = multiprocessing.Pool(
        processes=cpu_count - 1,
    )

    total_num_seq = sum(len(i) for i in sequence_mapping.values())

    try:
        peak_hits = pool.imap_unordered(
            func=_map_frag_compare,
            iterable=LenGen(
                gen=(
                    (
                        pep_query,
                        scan_mapping[pep_query],
                        sequence,
                        ms_two_data,
                        ms_data,
                        validation_data,
                        auto_maybe,
                    )
                    for pep_query, sequences in sequence_mapping.items()
                    for sequence in sequences
                ),
                len=total_num_seq,
            ),
        )

        for item in peak_hits:
            queue.put(item)

        pool.close()
    except:
        pool.terminate()
        raise
    finally:
        pool.join()
        _close_scans([ms_data, ms_two_data])


# Taken from https://stackoverflow.com/questions/1023038/
def lowpriority():
    """ Set the priority of the process to below-normal."""

    import sys
    try:
        sys.getwindowsversion()
    except AttributeError:
        isWindows = False
    else:
        isWindows = True

    if isWindows:
        # Based on:
        #   "Recipe 496767: Set Process Priority In Windows" on ActiveState
        #   http://code.activestate.com/recipes/496767/
        import win32api
        import win32process
        import win32con

        pid = win32api.GetCurrentProcessId()
        handle = win32api.OpenProcess(win32con.PROCESS_ALL_ACCESS, True, pid)
        win32process.SetPriorityClass(
            handle, win32process.BELOW_NORMAL_PRIORITY_CLASS,
        )
    else:
        import os

        os.nice(1)


def validate_spectra(
    search_path,
    raw_paths=None,
    scans_path=None,
    scan_list=None,
    mat_sessions=None,
    out_path=None,
    cpu_count=None,
    reprocess=False,
    auto_maybe=False,
):
    """
    Generate CAMV web page for validating spectra.

    Parameters
    ----------
    search_path : str
    raw_paths : list of str
    scans_path : str, optional
    scan_list : list of int, optional
    mat_sessions : list of str, optional
    out_path : str, optional
    cpu_count : int, optional
    reprocess : bool, optional
    auto_maybe : bool, optional
    """
    if cpu_count is None:
        try:
            cpu_count = multiprocessing.cpu_count()
        except NotImplementedError:
            cpu_count = 2

    LOGGER.debug("Validating data set with {} cpus".format(cpu_count))

    if raw_paths is None:
        raw_paths = []

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
    missing = []

    for base_raw in required_raws:
        if base_raw in base_raw_paths:
            continue

        for base_dir in [
            os.path.dirname(search_path),
            os.path.join(os.path.dirname(search_path), "..", "MS RAW")
        ]:
            local_raw_path = os.path.join(base_dir, base_raw)
            if os.path.exists(local_raw_path):
                raw_paths.append(local_raw_path)
                break
        else:
            missing.append(base_raw)

    if missing:
        raise Exception(
            "Unable to find {} in input RAW files: {}"
            .format(missing, base_raw_paths)
        )

    # Generate sequences
    LOGGER.info(
        "Generating all possible sequence-modification combinations."
    )
    pool = multiprocessing.Pool(
        processes=cpu_count,
    )
    try:
        sequence_mapping = dict(
            pool.imap_unordered(
                partial(_map_seq, limit_comb=not reprocess),
                pep_queries,
            )
        )
        pool.close()
    except:
        pool.terminate()
        raise
    finally:
        pool.join()

    total_num_seq = sum(len(i) for i in sequence_mapping.values())

    LOGGER.info(
        "Generating fragment ions for {} possible sequences."
        .format(
            total_num_seq,
        )
    )

    # Import validation data from CAMV-Matlab
    validation_data = {}

    if mat_sessions:
        LOGGER.info("Loading validation data from {}".format(mat_sessions))

        validation_data = camv_mat.load_mat_validation(mat_sessions)

        LOGGER.debug(
            "Loaded {} sequence validations from MATLAB sessions"
            .format(len(validation_data))
        )

    # Import raw scan data
    out_dir = tempfile.mkdtemp()

    LOGGER.info("Getting scan data.")

    scan_queries, ms_two_data, ms_data = scans.get_scan_data(
        raw_paths, pep_queries, out_dir,
    )

    try:
        LOGGER.info("Found data for {} scans".format(len(scan_queries)))

        scan_mapping = OrderedDict(
            (pep_query, scan_query)
            for pep_query, scan_query in zip(pep_queries, scan_queries)
        )

        # Generate fragments and assign peaks to fragments
        LOGGER.info(
            (
                "Comparing predicted to actual peaks for {} spectra "
                " ({} pep-scan combinations)."
            ).format(
                len(scan_mapping),
                total_num_seq,
            )
        )

        queue = multiprocessing.Queue()
        process = multiprocessing.Process(
            target=fill_map_frag_compare,
            args=(
                sequence_mapping,
                scan_mapping,
                ms_two_data,
                ms_data,
                queue,
                cpu_count,
                validation_data,
                auto_maybe,
            ),
        )
        process.start()

        # XXX: Determine SILAC precursor masses?

        # XXX: Remove precursor contaminated scans from validation list?

        # Check each assignment to each scan

        # Output data
        export.export_to_sql(
            os.path.splitext(out_path)[0] + ".db",
            queue, scan_mapping,
            search_path, raw_paths,
            total_num_seq=total_num_seq,
            reprocess=reprocess,
        )
    except:
        process.terminate()
        raise
    finally:
        process.join()
        _close_scans([ms_data, ms_two_data], out_dir=out_dir)

    LOGGER.info(
        "Exported {} total peptide-scan combinations"
        .format(total_num_seq)
    )
