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
    index = 0
    while index < total_num_seq:
        query, seq, choice, peaks, precursor_win, label_win = queue.get()

        LOGGER.debug(
            "Exporting: {}{} - {} - {}{}".format(
                index,
                " / {}".format(total_num_seq) if total_num_seq else "",
                query.scan,
                _pep_mod_name(
                    _extract_pep_seq(seq),
                    _extract_mods(seq),
                ),
                "- {}".format(choice) if choice else "",
            )
        )

        scan_query = scan_mapping[query]

        # Peptide sequence / modification data
        protein_ids = sql.insert_protein(cursor, query)
        protein_set_id = sql.insert_protein_set(cursor, query)
        peptide_id = sql.insert_peptide(cursor, query, protein_set_id)
        sql.insert_pep_prot(cursor, peptide_id, protein_ids)
        mod_state_id = sql.insert_mod_state(cursor, query, peptide_id)
        ptm_id = sql.insert_ptm(cursor, query, seq, mod_state_id)

        # Scan data
        quant_mz_id = sql.insert_quant_mz(cursor, query)
        file_id = sql.insert_file(cursor, query)
        scan_id = sql.insert_scans(
            cursor, query, scan_query,
            quant_mz_id,
            file_id,
            reprocessed=reprocess,
        )

        sql.insert_peaks(
            cursor, peaks, scan_id,
        )
        sql.insert_precursor_peaks(
            cursor, scan_query, precursor_win, scan_id,
        )
        sql.insert_quant_peaks(
            cursor, query, label_win, scan_id,
        )

        # PTM - Scan Mapping
        scan_ptm_id = sql.insert_scan_ptms(
            cursor, query, scan_id, ptm_id,
            choice=choice,
        )
        sql.insert_fragments(cursor, peaks, scan_ptm_id)

        # cursor.execute("COMMIT TRANSACTION")
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
