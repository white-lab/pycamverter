
import logging

from . import version

LOGGER = logging.getLogger("pycamv.migrations")


def run_migrations(from_version, to_version, cursor):
    if from_version == "1.0.0":
        if to_version == "1.1.0":
            return migrate_1_0_0_to_1_1_0(cursor)

    raise NotImplementedError(
        "Unable to migrate CAMV data from {} to {}".format(
            from_version, to_version,
        )
    )


def migrate_1_0_0_to_1_1_0(cursor):
    LOGGER.info("Migrating data from 1.0.0 to 1.1.0")

    cursor.executescript("""
    ALTER TABLE proteins ADD COLUMN protein_uniprot text;
    ALTER TABLE proteins ADD COLUMN full_sequence   text;

    ALTER TABLE protein_peptide ADD COLUMN peptide_offset integer;

    ALTER TABLE protein_sets ADD COLUMN protein_set_uniprot text;

    ALTER TABLE peptides ADD COLUMN protein_set_offsets text;
    """)
    cursor.connection.commit()
    cursor.executemany(
        """
        UPDATE camv_meta
        SET val=(?)
        WHERE key=(?)
        """,
        [
            (i[1], i[0])
            for i in {
                "pycamverterVersion": version.__version__,
                "camvDataVersion": "1.1.0",
            }.items()
        ]
    )
    cursor.connection.commit()
