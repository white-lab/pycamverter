
import logging
import sqlite3

from . import ms_labels, regexes, utils


LOGGER = logging.getLogger("pycamv.sql")

DB_EXTS = [".db", ".sql"]

CAMV_SCHEMA = """
PRAGMA foreign_keys = ON;
PRAGMA journal_mode = OFF;
PRAGMA synchronous = OFF;
PRAGMA temp_store = MEMORY;

-- Individual protein names (i.e. Src)
CREATE TABLE proteins
(
    protein_id              integer primary key autoincrement not null,
    protein_name            text,
    UNIQUE(protein_name)
);

-- Individual peptide sequences (i.e. IVLEYK)
CREATE TABLE peptides
(
    peptide_id              integer primary key autoincrement not null,
    peptide_seq             text,
    UNIQUE(peptide_seq)
);

-- Mapping between proteins and peptides
CREATE TABLE protein_peptide
(
    prot_pep_id             integer primary key autoincrement not null,
    peptide_id              integer,
    protein_id              integer,
    FOREIGN KEY(peptide_id) REFERENCES peptides(peptide_id),
    FOREIGN KEY(protein_id) REFERENCES proteins(protein_id),
    UNIQUE(peptide_id, protein_id)
);

-- Peptides with unmapped modifications (i.e. +1 pY)
CREATE TABLE mod_states
(
    mod_state_id            integer primary key autoincrement not null,
    peptide_id              integer,
    mod_desc                text,
    FOREIGN KEY(peptide_id) REFERENCES peptides(peptide_id),
    UNIQUE(peptide_id, mod_desc)
);

-- Peptides with modifications exactly positioned (i.e pY114)
CREATE TABLE ptms
(
    ptm_id                  integer primary key autoincrement not null,
    mod_state_id            integer,
    name                    text,
    full_name               text,
    FOREIGN KEY(mod_state_id) REFERENCES mod_states(mod_state_id),
    UNIQUE(mod_state_id, full_name)
);

-- Raw files sourced for each scan
CREATE TABLE files
(
    file_id                 integer primary key autoincrement not null,
    filename                text,
    UNIQUE(filename)
);

-- Blobs of scan data, for full ms2, quantification window, and precursor data
CREATE TABLE scan_data
(
    data_id                 integer primary key autoincrement not null,
    scan_id                 integer,
    data_type               text,
    data_blob               blob,
    FOREIGN KEY(scan_id) REFERENCES scans(scan_id),
    UNIQUE(scan_id, data_type)
);

-- Set of peaks that are used for quantitation
CREATE TABLE quant_mz
(
    quant_mz_id             integer primary key autoincrement not null,
    label_name              text,
    UNIQUE(label_name)
);

-- Individual peak / mz names used for quantitation
CREATE TABLE quant_mz_peaks
(
    quant_mz_peak_id        integer primary key autoincrement not null,
    quant_mz_id             integer not null,
    mz                      real not null,
    peak_name               text not null,
    FOREIGN KEY(quant_mz_id) REFERENCES quant_mz(quant_mz_id),
    UNIQUE(quant_mz_id, mz, peak_name)
);

CREATE TABLE scans
(
    scan_id                 integer primary key autoincrement not null,
    scan_num                integer not null,
    charge                  integer not null,
    collision_type          text not null,
    precursor_mz            real not null,
    isolation_window_lower  real,
    isolation_window_upper  real,
    quant_mz_id             integer,
    c13_num                 integer,
    file_id                 integer,
    FOREIGN KEY(quant_mz_id) REFERENCES quant_mz(quant_mz_id),
    FOREIGN KEY(file_id) REFERENCES files(file_id),
    UNIQUE(scan_num, file_id)
);

CREATE TABLE scan_ptms
(
    scan_ptm_id             integer primary key autoincrement not null,
    scan_id                 integer not null,
    ptm_id                  integer not null,
    choice                  text,
    mascot_score            real,
    FOREIGN KEY(scan_id) REFERENCES scans(scan_id)
    FOREIGN KEY(ptm_id) REFERENCES ptms(ptm_id)
);

CREATE TABLE fragments
(
    fragment_id             integer primary key autoincrement not null,
    peak_id                 integer not null,
    scan_ptm_id             integer not null,
    name                    text not null,
    display_name            text not null,
    mz                      real not null,
    best                    integer,
    ion_type                text,
    ion_pos                 integer,
    FOREIGN KEY(scan_ptm_id) REFERENCES scan_ptms(scan_ptm_id)
    -- UNIQUE(ptm_id, name)
);
"""


def _insert_or_update_row(cursor, table, id, data, unique_on=None):
    if unique_on is None:
        unique_on = list(data.keys())

    row_id = None

    for row in cursor.execute(
        """
        SELECT ({})
        FROM {}
        WHERE {}
        """.format(
            id,
            table,
            " and ".join("{}=(?)".format(i) for i in unique_on)
        ),
        [data[i] for i in unique_on],
    ):
        row_id = row[0]

    if row_id is None:
        cursor.execute(
            """
            INSERT OR IGNORE INTO {}
            ({})
            VALUES ({})
            """.format(
                table,
                ", ".join(data.keys()),
                ",".join("?" for i in data.keys()),
            ),
            list(data.values()),
        )
        row_id = cursor.lastrowid

    return row_id


def insert_protein(cursor, query):
    return [
        _insert_or_update_row(
            cursor, "proteins", "protein_id",
            {
                "protein_name": prot_desc,
            },
        )
        for prot_desc in query.prot_descs
    ]


def insert_peptide(cursor, query):
    return _insert_or_update_row(
        cursor, "peptides", "peptide_id",
        {
            "peptide_seq": query.pep_seq,
        },
    )


def insert_pep_prot(cursor, pep_id, prot_ids):
    return [
        _insert_or_update_row(
            cursor, "protein_peptide", "prot_pep_id",
            {
                "peptide_id": pep_id,
                "protein_id": prot_id,
            },
        )
        for prot_id in prot_ids
    ]


def _get_labels_mz(query):
    def _within(val, bounds):
        return val >= bounds[0] and val <= bounds[1]
    return [
        (mz, name)
        for mod in set(query.get_label_mods)
        for mz, name in zip(
            ms_labels.LABEL_MASSES.get(mod, []),
            ms_labels.LABEL_NAMES.get(mod, []),
        )
        if _within(mz, ms_labels.LABEL_MZ_WINDOW.get(mod, [mz, mz]))
    ]


def insert_quant_mz(cursor, query):
    quant_mz_id = _insert_or_update_row(
        cursor, "quant_mz", "quant_mz_id",
        {
            "label_name": ";".join(sorted(set(query.get_label_mods))),
        },
    )

    cursor.executemany(
        """
        INSERT OR IGNORE INTO quant_mz_peaks
        (
            quant_mz_id,
            mz,
            peak_name
        )
        VALUES (?, ?, ?)
        """,
        [
            (
                quant_mz_id,
                mz,
                name
            )
            for mz, name in _get_labels_mz(query)
        ],
    )

    return quant_mz_id


def insert_file(cursor, query):
    return _insert_or_update_row(
        cursor, "files", "file_id",
        {
            "filename": query.filename,
        },
    )


def _peaks_to_blob(peaks):
    return sqlite3.Binary(
        ";".join(
            "{},{}".format(
                i.mz if hasattr(i, "mz") else i[0],
                i.intensity if hasattr(i, "intensity") else i[1],
            )
            for i in peaks
        ).encode("utf-8")
    )


def insert_peaks(cursor, peaks, scan_id):
    return _insert_or_update_row(
        cursor, "scan_data", "data_id",
        {
            "scan_id": scan_id,
            "data_type": "ms2",
            "data_blob": _peaks_to_blob(peaks),
        },
        ["scan_id", "data_type"],
    )


def insert_precursor_peaks(cursor, scan_query, precursor_win, scan_id):
    return _insert_or_update_row(
        cursor, "scan_data", "data_id",
        {
            "scan_id": scan_id,
            "data_type": "precursor",
            "data_blob": _peaks_to_blob(precursor_win),
        },
        ["scan_id", "data_type"],
    )


def insert_quant_peaks(cursor, query, label_win, scan_id):
    return _insert_or_update_row(
        cursor, "scan_data", "data_id",
        {
            "scan_id": scan_id,
            "data_type": "quant",
            "data_blob": _peaks_to_blob(label_win),
        },
        ["scan_id", "data_type"],
    )


def insert_scans(
    cursor, query, scan_query,
    quant_mz_id, file_id,
):
    return _insert_or_update_row(
        cursor, "scans", "scan_id",
        {
            "scan_num": query.scan,
            "charge": query.pep_exp_z,
            "collision_type": scan_query.collision_type,
            "precursor_mz": query.pep_exp_mz,
            "isolation_window_lower": scan_query.window_offset[0],
            "isolation_window_upper": scan_query.window_offset[1],
            "c13_num": scan_query.c13_num,
            "quant_mz_id": quant_mz_id,
            "file_id": file_id,
        },
        ["scan_num", "file_id"],
    )


def insert_scan_ptms(cursor, query, scan_id, ptm_id):
    return _insert_or_update_row(
        cursor, "scan_ptms", "scan_ptm_id",
        {
            "scan_id": scan_id,
            "ptm_id": ptm_id,
            "mascot_score": query.pep_score,
        },
        ["scan_id", "ptm_id"],
    )


def _get_mods_description(pep_seq, mod_state):
    return ("+ " if mod_state else "") + " - ".join(
        "{} {}{}".format(count, mod[0].lower(), "".join(letters))
        for count, mod, letters in mod_state
    )


def insert_mod_state(cursor, query, peptide_id):
    return _insert_or_update_row(
        cursor, "mod_states", "mod_state_id",
        {
            "peptide_id": peptide_id,
            "mod_desc": _get_mods_description(
                query.pep_seq,
                tuple(query.pep_var_mods)
            ),
        },
    )


def _pep_mod_name(pep_seq, mods):
    return "".join(
        letter.lower() if mods else letter.upper()
        for letter, mods in zip(pep_seq, mods[1:-1])
    )


def _pep_mod_full_name(pep_seq, mods):
    return "-".join(
        letter + (
            " ({})".format(",".join(mods))
            if mods else
            ""
        )
        for letter, mods in zip(
            ["N-term"] + list(pep_seq) + ["C-term"],
            mods,
        )
    )


def _extract_mods(sequence):
    return tuple(
        tuple(mods)
        for _, mods in sequence
    )


def insert_ptm(cursor, query, seq, mod_state_id):
    mods = _extract_mods(seq)
    return _insert_or_update_row(
        cursor, "ptms", "ptm_id",
        {
            "mod_state_id": mod_state_id,
            "name": _pep_mod_name(query.pep_seq, mods),
            "full_name": _pep_mod_full_name(query.pep_seq, mods),
        },
        ["mod_state_id", "full_name"],
    )


def _ion_type_pos(name):
    name_match = regexes.RE_BY_ION_POS.match(name)
    if name_match:
        ion_pos = int(name_match.group(2))
        ion_type = "b" if name_match.group(1) in "abc" else "y"
    else:
        ion_type, ion_pos = None, None

    return ion_type, ion_pos


def insert_fragments(cursor, peaks, scan_ptm_id):
    gen = (
        (
            scan_ptm_id,
            peak_index,
            name,
            utils.rewrite_ion_name(name),
            mz,
            name == peak_hit.name,
        ) + _ion_type_pos(name)
        for peak_index, peak_hit in enumerate(peaks)
        if peak_hit.match_list
        for name, (mz, _) in peak_hit.match_list.items()
    )
    cursor.execute("BEGIN")
    for tup in gen:
        cursor.execute(
            """
            INSERT OR IGNORE INTO fragments
            (
                scan_ptm_id,
                peak_id,
                name,
                display_name,
                mz,
                best,
                ion_type,
                ion_pos
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            tup,
        )
    cursor.connection.commit()
