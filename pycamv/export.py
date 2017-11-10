"""
Provides functionality for exporting processed data to .camv files (JSON
format).
"""

from __future__ import absolute_import, division

import errno
import logging
import os
import sqlite3
from time import time

from . import sql

LOGGER = logging.getLogger("pycamv.export")

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError


def _extract_pep_seq(sequence):
    return "".join(
        letter
        for letter, _ in sequence
        if letter not in ["N-term", "C-term", "C=O"]
    )


def _extract_mods(sequence):
    return tuple(
        tuple(mods)
        for _, mods in sequence
    )


def _pep_mod_name(pep_seq, mods):
    return "".join(
        letter.lower() if mods else letter.upper()
        for letter, mods in zip(pep_seq, mods[1:-1])
    )


def export_to_sql(
    out_path, queue, scan_mapping, search_path, raw_paths,
    reprocess=False, total_num_seq=None,
):
    assert os.path.splitext(out_path)[1] in sql.DB_EXTS

    if not reprocess:
        try:
            os.remove(out_path)
        except FileNotFoundError as e:
            if type(e) == IOError and e.errno != errno.EEXIST:
                raise e

    db = sqlite3.connect(out_path, isolation_level="EXCLUSIVE")
    cursor = db.cursor()

    sql.create_tables(cursor)
    sql.run_migrations(cursor)

    sql.insert_path_data(cursor, search_path, raw_paths)

    total = time()

    # frag_map = defaultdict(list)
    index = 1
    while index <= total_num_seq:
        pep_query, seq, choice, peaks, precursor_win, label_win = queue.get()

        scan_query = scan_mapping[pep_query]

        LOGGER.debug(
            "Exporting: {} - {}{} {}".format(
                pep_query.scan,
                _pep_mod_name(
                    _extract_pep_seq(seq),
                    _extract_mods(seq),
                ),
                "- {}".format(choice) if choice else "",
                "({} / {})".format(index, total_num_seq)
                if total_num_seq else
                "(# {})".format(index),
            )
        )

        # Protein
        protein_ids = sql.insert_protein(cursor, pep_query)
        protein_set_id = sql.insert_protein_set(cursor, pep_query)

        # Peptide
        peptide_id = sql.insert_peptide(cursor, pep_query, protein_set_id)
        sql.insert_pep_prot(
            cursor,
            peptide_id, protein_ids,
            pep_query.pep_offsets,
        )

        # Modification data
        mod_state_id = sql.insert_mod_state(cursor, pep_query, peptide_id)
        ptm_id = sql.insert_ptm(cursor, pep_query, seq, mod_state_id)

        # Scan header data
        quant_mz_id = sql.insert_quant_mz(cursor, pep_query)
        file_id = sql.insert_file(cursor, pep_query)
        scan_id = sql.insert_scans(
            cursor, pep_query, scan_query,
            quant_mz_id,
            file_id,
            reprocessed=reprocess,
        )

        # Scan peak data
        sql.insert_peaks(
            cursor, peaks, scan_id,
        )
        sql.insert_precursor_peaks(
            cursor, scan_query, precursor_win, scan_id,
        )
        sql.insert_quant_peaks(
            cursor, pep_query, label_win, scan_id,
        )

        # PTM - Scan Mapping
        scan_ptm_id = sql.insert_scan_ptms(
            cursor, pep_query, scan_id, ptm_id,
            choice=choice,
        )
        sql.insert_fragments(cursor, peaks, scan_ptm_id)

        db.commit()
        LOGGER.debug(
            "done - avg: {:.3f} sec".format((time() - total) / (index + 1))
        )

        index += 1

    LOGGER.debug(
        "total: {:.3f} min ({:.3f} sec / peptide)"
        .format((time() - total) / 60, (time() - total) / (index + 1))
    )

    db.close()
