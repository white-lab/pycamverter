
import logging
import sqlite3

from . import gen_sequences, migrations, ms_labels, regexes, utils, version


LOGGER = logging.getLogger("pycamv.sql")

DB_EXTS = [".db", ".sql"]
DATA_VERSION = "1.1.0"

CAMV_SCHEMA = """
PRAGMA foreign_keys = ON;
PRAGMA journal_mode = OFF;
PRAGMA synchronous = OFF;
PRAGMA temp_store = MEMORY;

-- Individual protein names (i.e. Src)
CREATE TABLE IF NOT EXISTS proteins
(
    protein_id              integer primary key autoincrement not null,
    protein_name            text,
    protein_accession       text,
    protein_uniprot         text,
    full_sequence           text,
    UNIQUE(protein_accession)
);

-- Pre-processed set of protein names (i.e. Cdk1 / Cdk2 / Cdk3)
CREATE TABLE IF NOT EXISTS protein_sets
(
    protein_set_id          integer primary key autoincrement not null,
    protein_set_name        text,
    protein_set_accession   text,
    protein_set_uniprot     text,
    UNIQUE(protein_set_accession)
);

-- Individual peptide sequences (i.e. IVLEYK)
CREATE TABLE IF NOT EXISTS peptides
(
    peptide_id              integer primary key autoincrement not null,
    protein_set_id          integer not null,
    protein_set_offsets     text,
    peptide_seq             text,
    FOREIGN KEY(protein_set_id) REFERENCES protein_sets(protein_set_id)
    UNIQUE(peptide_seq)
);

-- Mapping between proteins and peptides
CREATE TABLE IF NOT EXISTS protein_peptide
(
    prot_pep_id             integer primary key autoincrement not null,
    peptide_id              integer not null,
    protein_id              integer not null,
    peptide_offset          integer,
    FOREIGN KEY(peptide_id) REFERENCES peptides(peptide_id),
    FOREIGN KEY(protein_id) REFERENCES proteins(protein_id),
    UNIQUE(peptide_id, protein_id)
);

-- Peptides with unmapped modifications (i.e. +1 pY)
CREATE TABLE IF NOT EXISTS mod_states
(
    mod_state_id            integer primary key autoincrement not null,
    peptide_id              integer not null,
    mod_desc                text,
    num_comb                integer,
    FOREIGN KEY(peptide_id) REFERENCES peptides(peptide_id),
    UNIQUE(peptide_id, mod_desc)
);

-- Peptides with modifications exactly positioned (i.e pY114)
CREATE TABLE IF NOT EXISTS ptms
(
    ptm_id                  integer primary key autoincrement not null,
    mod_state_id            integer not null,
    name                    text,
    full_name               text,
    FOREIGN KEY(mod_state_id) REFERENCES mod_states(mod_state_id),
    UNIQUE(mod_state_id, full_name)
);

-- Raw files sourced for each scan
CREATE TABLE IF NOT EXISTS files
(
    file_id                 integer primary key autoincrement not null,
    filename                text,
    UNIQUE(filename)
);

-- Blobs of scan data, for full ms2, quantification window, and precursor data
CREATE TABLE IF NOT EXISTS scan_data
(
    data_id                 integer primary key autoincrement not null,
    scan_id                 integer not null,
    data_type               text,
    data_blob               blob,
    FOREIGN KEY(scan_id) REFERENCES scans(scan_id),
    UNIQUE(scan_id, data_type)
);

-- Set of peaks that are used for quantitation
CREATE TABLE IF NOT EXISTS quant_mz
(
    quant_mz_id             integer primary key autoincrement not null,
    label_name              text,
    UNIQUE(label_name)
);

-- Individual peak / mz names used for quantitation
CREATE TABLE IF NOT EXISTS quant_mz_peaks
(
    quant_mz_peak_id        integer primary key autoincrement not null,
    quant_mz_id             integer not null,
    mz                      real not null,
    peak_name               text not null,
    FOREIGN KEY(quant_mz_id) REFERENCES quant_mz(quant_mz_id),
    UNIQUE(quant_mz_id, mz, peak_name)
);

CREATE TABLE IF NOT EXISTS scans
(
    scan_id                 integer primary key autoincrement not null,
    scan_num                integer not null,
    charge                  integer not null,
    pep_exp_mz              integer not null,
    collision_type          text not null,
    precursor_mz            real not null,
    isolation_window_lower  real,
    isolation_window_upper  real,
    quant_mz_id             integer,
    c13_num                 integer,
    file_id                 integer,
    truncated               integer,
    FOREIGN KEY(quant_mz_id) REFERENCES quant_mz(quant_mz_id),
    FOREIGN KEY(file_id) REFERENCES files(file_id),
    UNIQUE(scan_num, file_id)
);

CREATE TABLE IF NOT EXISTS scan_ptms
(
    scan_ptm_id             integer primary key autoincrement not null,
    scan_id                 integer not null,
    ptm_id                  integer not null,
    choice                  text,
    mascot_score            real,
    FOREIGN KEY(scan_id) REFERENCES scans(scan_id),
    FOREIGN KEY(ptm_id) REFERENCES ptms(ptm_id),
    UNIQUE(scan_id, ptm_id)
);

CREATE TABLE IF NOT EXISTS fragments
(
    fragment_id             integer primary key autoincrement not null,
    peak_id                 integer not null,
    scan_ptm_id             integer not null,
    name                    text not null,
    display_name            text not null,
    mz                      real not null,
    intensity               real not null,
    best                    integer,
    ion_type                text,
    ion_pos                 integer,
    FOREIGN KEY(scan_ptm_id) REFERENCES scan_ptms(scan_ptm_id)
    UNIQUE(peak_id, scan_ptm_id, name)
);
CREATE INDEX IF NOT EXISTS fragments_idx ON fragments(scan_ptm_id);
CREATE INDEX IF NOT EXISTS fragments_idx ON fragments(peak_id, scan_ptm_id);

CREATE TABLE IF NOT EXISTS camv_meta
(
    key                     text not null,
    val                     text,
    UNIQUE(key, val)
);
"""


def create_tables(cursor):
    cursor.executescript(CAMV_SCHEMA)
    insert_camv_meta(cursor)
    cursor.connection.commit()


def run_migrations(cursor):
    rows = cursor.execute(
        """
        SELECT camv_meta.key, camv_meta.val
        FROM camv_meta
        WHERE camv_meta.key IN (?)
        """,
        [
            "camvDataVersion",
        ],
    )

    camv_data_version = None

    for key, val in rows:
        if key == "camvDataVersion":
            camv_data_version = val
        else:
            raise Exception(
                "Unexpected key in camv_meta: {} - {}".format(key, val)
            )

    if camv_data_version is None:
        raise Exception(
            "Unable to determine CAMV data version of existing database",
        )

    if camv_data_version != DATA_VERSION:
        migrations.run_migrations(camv_data_version, DATA_VERSION, cursor)


def _insert_or_update_row(
    cursor, table, id, data,
    unique_on=None, update=False,
):
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
    elif update:
        keys = [i for i in data.keys() if i not in unique_on]
        cursor.execute(
            """
            UPDATE {}
            SET {}
            WHERE {}
            """.format(
                table,
                ", ".join("{}=?".format(i) for i in keys),
                " and ".join("{}=?".format(i) for i in unique_on),
            ),
            [data[i] for i in keys] + [data[i] for i in unique_on],
        )

    return row_id


def insert_protein(cursor, query):
    return [
        _insert_or_update_row(
            cursor, "proteins", "protein_id",
            {
                "protein_name": prot_desc,
                "protein_accession": access,
                "protein_uniprot": uniprot,
                "full_sequence": seq,
            },
            unique_on=["protein_accession"],
            update=True,
        )
        for prot_desc, access, uniprot, seq in zip(
            query.prot_descs, query.accessions,
            query.uniprot_accessions, query.full_seqs,
        )
    ]


def insert_protein_set(cursor, query):
    return _insert_or_update_row(
        cursor, "protein_sets", "protein_set_id",
        {
            "protein_set_name": " / ".join(
                i[0]
                for i in sorted(
                    zip(query.prot_descs, query.accessions),
                    key=lambda x: x[1]
                )
            ),
            "protein_set_accession": " / ".join(
                sorted(query.accessions)
            ),
            "protein_set_uniprot": ";".join(
                i[0]
                for i in sorted(
                    zip(query.uniprot_accessions, query.accessions),
                    key=lambda x: x[1],
                )
            ),
        },
        unique_on=["protein_set_accession"],
        update=True,
    )


def insert_peptide(cursor, query, prot_set_id):
    return _insert_or_update_row(
        cursor, "peptides", "peptide_id",
        {
            "peptide_seq": query.pep_seq,
            "protein_set_id": prot_set_id,
            "protein_set_offsets": ";".join(
                str(i[0])
                for i in sorted(
                    zip(query.pep_offsets, query.accessions),
                    key=lambda x: x[1],
                )
            ),
        },
        unique_on=["peptide_seq", "protein_set_id"],
        update=True,
    )


def insert_pep_prot(cursor, pep_id, prot_ids, prot_offsets):
    return [
        _insert_or_update_row(
            cursor, "protein_peptide", "prot_pep_id",
            {
                "peptide_id": pep_id,
                "protein_id": prot_id,
                "peptide_offset": offset,
            },
            unique_on=[
                "peptide_id",
                "protein_id",
            ],
            update=True,
        )
        for prot_id, offset in zip(prot_ids, prot_offsets)
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
        unique_on=["scan_id", "data_type"],
    )


def insert_precursor_peaks(cursor, scan_query, precursor_win, scan_id):
    return _insert_or_update_row(
        cursor, "scan_data", "data_id",
        {
            "scan_id": scan_id,
            "data_type": "precursor",
            "data_blob": _peaks_to_blob(precursor_win),
        },
        unique_on=["scan_id", "data_type"],
    )


def insert_quant_peaks(cursor, query, label_win, scan_id):
    return _insert_or_update_row(
        cursor, "scan_data", "data_id",
        {
            "scan_id": scan_id,
            "data_type": "quant",
            "data_blob": _peaks_to_blob(label_win),
        },
        unique_on=["scan_id", "data_type"],
    )


def insert_scans(
    cursor, query, scan_query,
    quant_mz_id, file_id,
    reprocessed=False, c13_num=None,
):
    return _insert_or_update_row(
        cursor, "scans", "scan_id",
        {
            "scan_num": query.scan,
            "charge": query.pep_exp_z,
            "pep_exp_mz": query.pep_exp_mz,
            "collision_type": scan_query.collision_type,
            "precursor_mz": scan_query.isolation_mz,
            "isolation_window_lower": scan_query.window_offset[0],
            "isolation_window_upper": scan_query.window_offset[1],
            "c13_num": c13_num or scan_query.c13_num,
            "quant_mz_id": quant_mz_id,
            "file_id": file_id,
            "truncated": (
                not reprocessed and query.num_comb > gen_sequences.MAX_NUM_COMB
            ),
        },
        unique_on=["scan_num", "file_id"],
        update=reprocessed,
    )


def insert_scan_ptms(cursor, query, scan_id, ptm_id, choice=None):
    return _insert_or_update_row(
        cursor, "scan_ptms", "scan_ptm_id",
        {
            "scan_id": scan_id,
            "ptm_id": ptm_id,
            "mascot_score": query.pep_score,
            "choice": choice,
        },
        unique_on=["scan_id", "ptm_id"],
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
                tuple(query.pep_var_mods),
            ),
            "num_comb": query.num_comb,
        },
        unique_on=[
            "peptide_id",
            "mod_desc",
        ]
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
        unique_on=["mod_state_id", "full_name"],
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
            peak_hit.intensity,
            name == peak_hit.name,
        ) + _ion_type_pos(name)
        for peak_index, peak_hit in enumerate(peaks)
        if peak_hit.match_list
        for name, (mz, _) in peak_hit.match_list.items()
    )
    cursor.executemany(
        """
        INSERT OR IGNORE INTO fragments
        (
            scan_ptm_id,
            peak_id,
            name,
            display_name,
            mz,
            intensity,
            best,
            ion_type,
            ion_pos
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        gen,
    )


def insert_camv_meta(cursor):
    cursor.executemany(
        """
        INSERT INTO camv_meta
        (
            key,
            val
        ) SELECT (?), (?)
        WHERE NOT EXISTS (SELECT 1 FROM camv_meta WHERE key=(?))
        """,
        [
            (i[0], i[1], i[0])
            for i in {
                "pycamverterVersion": version.__version__,
                "camvDataVersion": DATA_VERSION,
            }.items()
        ],
    )


def insert_path_data(cursor, search_path, raw_paths):
    cursor.executemany(
        """
        INSERT OR IGNORE INTO camv_meta
        (
            key,
            val
        ) VALUES (?, ?)
        """,
        [
            ("search_path", search_path),
        ] + [
            ("raw_path", raw_path)
            for raw_path in raw_paths
        ],
    )
