"""
Provides functionality for interacting with ProteomeDiscoverer data.
"""

from __future__ import absolute_import, division

from collections import Counter, defaultdict
from datetime import datetime
import logging
import os
import sqlite3
import re
import xml.etree.ElementTree as ET

from . import masses, regexes, pep_query


LOGGER = logging.getLogger("pycamv.discoverer")
RE_PSP = re.compile("(\w+)\((\d+)\): ([\d\.]+)")


def _parse_letters(letters):
    """
    Turns a string of residue letters (i.e. "STY") into a list.

    Includes special provisions for N- and C-term modifications.

    Parameters
    ----------
    letters : str

    Returns
    -------
    tuple of str
    """
    if letters in ["N-term", "C-term"]:
        return (letters,)

    return tuple(letters)


def _get_fixed_var_mods_pd21(conn):
    query = conn.cursor().execute(
        """
        SELECT
        Workflows.WorkflowXML
        FROM
        Workflows
        """
    )

    fixed_mods = []
    var_mods = []

    for xml, in query:
        root = ET.fromstring(xml)
        params = root.findall(
            "WorkflowTree/WorkflowNode/"
            "ProcessingNodeParameters/ProcessingNodeParameter"
        )
        for param in params:
            if not param.get("IsValueSet", False) or not param.text:
                continue

            if param.get("Name", "").startswith("DynMod_"):
                if param.text.strip():
                    var_mods.append(param.text.strip())
            elif param.get("Name", "").startswith("StaticMod_"):
                if param.text.strip():
                    fixed_mods.append(param.text.strip())

    return fixed_mods, var_mods


def _get_fixed_var_mods(conn):
    try:
        query = conn.cursor().execute(
            """
            SELECT
            ProcessingNodeParameters.ParameterName,
            ProcessingNodeParameters.ParameterValue
            FROM
            ProcessingNodeParameters
            """,
        )
    except sqlite3.OperationalError:
        return _get_fixed_var_mods_pd21(conn)

    fixed_mods = []
    var_mods = []

    for name, val in query:
        if name.startswith("StaticMod_"):
            if val.strip():
                fixed_mods.append(val.strip())
        elif name.startswith("DynMod_"):
            if val.strip():
                var_mods.append(val.strip())

    return fixed_mods, var_mods


def _find_mod(abbrev, letter, pot_mods):
    for pot_mod in pot_mods:
        if pot_mod[0] == abbrev and letter in pot_mod[1]:
            return pot_mod

    return None


def _count_mods(mod_list):
    return [
        (count, abbrev, _parse_letters(letters))
        for (abbrev, letters), count in Counter(mod_list).items()
    ]


def _get_pep_mods(conn, pep_id, pep_seq, var_mods, fixed_mods):
    aa_mods = conn.cursor().execute(
        """
        SELECT
        Peptides.SearchEngineRank,
        AminoAcidModifications.ModificationName,
        PeptidesAminoAcidModifications.Position
        FROM Peptides
        JOIN PeptidesAminoAcidModifications
        JOIN AminoAcidModifications
        ON PeptidesAminoAcidModifications.AminoAcidModificationID=
        AminoAcidModifications.AminoAcidModificationID
        WHERE
        Peptides.PeptideID=:pepID AND
        PeptidesAminoAcidModifications.PeptideID=:pepID
        """,
        {
            "pepID": pep_id,
        }
    )

    pep_var_mods = []
    pep_fixed_mods = []
    rank_pos = defaultdict(set)

    for rank, abbrev, pos in sorted(
        aa_mods,
        key=lambda x: 0 if x[1].startswith("Mapping") else 1,
    ):
        letter = pep_seq[pos]
        rank_pos[rank].add((pos, abbrev))

        mod = _find_mod(abbrev, letter, var_mods)

        if mod:
            pep_var_mods.append(mod)
            continue

        mod = _find_mod(abbrev, letter, fixed_mods)

        if mod:
            pep_fixed_mods.append(mod)
            continue

        # Weird ProteomeDiscoverer setting for proteins with unknown amino
        # acids. We can remap them to their correct amino acid:
        if (
            letter == "X" and
            (letter, abbrev) in masses.MODIFICATIONS and
            abbrev.startswith("Mapping")
        ):
            pep_seq = pep_seq[:pos] + abbrev[-1] + pep_seq[pos + 1:]
            continue

        raise Exception(
            "Unexpected modification: {} {}".format(letter, abbrev)
        )

    term_mods = conn.cursor().execute(
        """
        SELECT
        AminoAcidModifications.ModificationName,
        AminoAcidModifications.PositionType
        FROM PeptidesTerminalModifications
        JOIN AminoAcidModifications
        ON PeptidesTerminalModifications.TerminalModificationID=
        AminoAcidModifications.AminoAcidModificationID
        WHERE
        PeptidesTerminalModifications.PeptideID=:pepID
        """,
        {
            "pepID": pep_id,
        }
    )

    for abbrev, pos_type in term_mods:
        letter = "N-term" if pos_type == 1 else "C-term"
        mod = _find_mod(abbrev, letter, var_mods)

        if mod:
            pep_var_mods.append(mod)
            continue

        mod = _find_mod(abbrev, letter, fixed_mods)

        if mod:
            pep_fixed_mods.append(mod)
            continue

        raise Exception(
            "Unexpected modification: {} {}".format(letter, abbrev)
        )

    pep_var_mods = _count_mods(pep_var_mods)
    pep_fixed_mods = _count_mods(pep_fixed_mods)

    return pep_seq, pep_var_mods, pep_fixed_mods, rank_pos


def _get_prot_info(conn, peptide_id):
    proteins = conn.cursor().execute(
        """
        SELECT
        ProteinAnnotations.Description,
        Proteins.Sequence
        FROM PeptidesProteins
        JOIN ProteinAnnotations
        ON ProteinAnnotations.ProteinID=PeptidesProteins.ProteinID
        JOIN Proteins
        ON Proteins.ProteinID=PeptidesProteins.ProteinID
        WHERE
        PeptidesProteins.PeptideID=:pepID
        """,
        {
            "pepID": peptide_id,
        }
    )

    (
        accessions,
        descriptions,
        uniprot_accessions,
        full_seqs,
    ) = [], [], [], []

    for full_prot_desc, sequence in proteins:
        uniprot_accession, accession, prot_desc = \
            regexes.RE_DISCOVERER_DESCRIPTION.match(
                full_prot_desc,
            ).group(1, 2, 3)

        if not accession:
            raise Exception(
                "Unable to find accession ID for {}".format(full_prot_desc)
            )

        accessions.append(accession)
        descriptions.append(prot_desc)
        uniprot_accessions.append(uniprot_accession)
        full_seqs.append(sequence)

    return (
        accessions,
        descriptions,
        uniprot_accessions,
        full_seqs,
    )


def get_quant_scan(conn, scan):
    query = conn.cursor().execute(
        """
        SELECT
        quant_scans.FirstScan
        FROM ReporterIonQuanResultsSearchSpectra mapping
        JOIN SpectrumHeaders quant_scans
        ON mapping.SpectrumID=quant_scans.SpectrumID
        JOIN SpectrumHeaders search_scans
        ON search_scans.SpectrumID=mapping.SearchSpectrumID
        WHERE search_scans.FirstScan=:scanId
        """,
        {
            "scanId": scan,
        }
    )
    for quant_scan, in query:
        return quant_scan


def _get_peptide_queries(conn, fixed_mods, var_mods):
    out = []
    # index = 0
    # scan_used = {}
    fixed_mods = [
        regexes.RE_DYN_MODS.match(i).group(3, 4)
        for i in fixed_mods
    ]
    var_mods = [
        regexes.RE_DYN_MODS.match(i).group(3, 4)
        for i in var_mods
    ]

    try:
        # ProteomeDiscoverer == 1.4
        query = conn.cursor().execute(
            """
            SELECT
            Peptides.PeptideID,
            Peptides.Sequence,
            PeptideScores.ScoreValue,
            SpectrumHeaders.FirstScan,
            Masspeaks.Mass,
            Masspeaks.Charge,
            SpectrumHeaders.Charge,
            FileInfos.FileName
            FROM Peptides
            JOIN PeptideScores
            ON PeptideScores.PeptideID=Peptides.PeptideID
            JOIN PeptidesProteins
            ON PeptidesProteins.PeptideID=Peptides.PeptideID
            JOIN ProteinAnnotations
            ON ProteinAnnotations.ProteinID=PeptidesProteins.ProteinID
            JOIN SpectrumHeaders
            ON SpectrumHeaders.SpectrumID=Peptides.SpectrumID
            JOIN FileInfos
            ON FileInfos.FileID=MassPeaks.FileID
            JOIN Masspeaks
            ON Masspeaks.MassPeakID=SpectrumHeaders.MassPeakID
            WHERE
            Peptides.SearchEngineRank=1
            """
        )
    except sqlite3.OperationalError:
        # ProteomeDiscoverer == 2.1
        query = conn.cursor().execute(
            """
            SELECT
            Peptides.PeptideID,
            Peptides.Sequence,
            PeptideScores.ScoreValue,
            SpectrumHeaders.FirstScan,
            Masspeaks.Mass,
            Masspeaks.Charge,
            SpectrumHeaders.Charge,
            WorkflowInputFiles.FileName
            FROM Peptides
            JOIN PeptideScores
            ON PeptideScores.PeptideID=Peptides.PeptideID
            JOIN PeptidesProteins
            ON PeptidesProteins.PeptideID=Peptides.PeptideID
            JOIN ProteinAnnotations
            ON ProteinAnnotations.ProteinID=PeptidesProteins.ProteinID
            JOIN SpectrumHeaders
            ON SpectrumHeaders.SpectrumID=Peptides.SpectrumID
            JOIN WorkflowInputFiles
            ON WorkflowInputFiles.FileID=MassPeaks.FileID
            JOIN Masspeaks
            ON Masspeaks.MassPeakID=SpectrumHeaders.MassPeakID
            WHERE
            Peptides.SearchEngineRank=1
            """
        )

    psp_vals = _get_phosphors_psp_vals(conn.cursor())
    changed_peptides = 0
    ambiguous_peptides = 0

    for (
        pep_id, pep_seq, pep_score,
        scan, exp_mz, mass_z, exp_z, filename,
    ) in query:
        if exp_z < 1:
            LOGGER.warning(
                "Charge: {}, {} for {} (Scan {})"
                .format(exp_z, mass_z, pep_seq, scan)
            )

        pep_seq, pep_var_mods, pep_fixed_mods, rank_pos = _get_pep_mods(
            conn, pep_id, pep_seq, var_mods, fixed_mods,
        )

        quant_scan = get_quant_scan(conn, scan)

        (
            accessions,
            descriptions,
            uniprot_accessions,
            full_seqs,
        ) = _get_prot_info(conn, pep_id)

        if pep_id in psp_vals:
            rank_pos, reassigned, ambiguous = _reassign_rank(
                pep_fixed_mods + pep_var_mods,
                rank_pos,
                psp_vals[pep_id],
            )

            if reassigned:
                changed_peptides += 1

            if ambiguous:
                ambiguous_peptides += 1

        out.append(
            pep_query.PeptideQuery(
                accessions=accessions,
                prot_descs=descriptions,
                uniprot_accessions=uniprot_accessions,
                full_seqs=full_seqs,
                query=pep_id,
                filename=filename,
                pep_score=pep_score,
                pep_exp_mz=exp_mz,
                pep_exp_z=exp_z,
                pep_seq=pep_seq,
                pep_var_mods=pep_var_mods,
                pep_fixed_mods=pep_fixed_mods,
                scan=scan,
                rank_pos=rank_pos,
                quant_scan=quant_scan,
            )
        )

    out = sorted(out, key=lambda x: x.scan)

    LOGGER.info(
        "Reassigned {} top modification assignments using phosphoRS"
        .format(changed_peptides)
    )

    return out


def _is_pmod(mod):
    pos, mod_type = mod
    return (
        mod_type in ["Phospho"]
    )


def _sort_mods(mods):
    return tuple(
        sorted(
            mods,
            key=lambda x: (x[0], x[1]),
        )
    )


def _reassign_rank(mods, rank_pos, psp_val):
    reassigned = False
    ambiguous = False

    # phophoRS example format: "T(4): 99.6; S(6): 0.4; S(10): 0.0"
    # Error messages include: "Too many isoforms"
    psp_val = [
        RE_PSP.match(i.strip())
        for i in psp_val.split(";")
    ]
    psp_val = [
        i.groups()
        for i in psp_val
        if i
    ]
    psp_val = [
        (i[0], int(i[1]), float(i[2]))
        for i in psp_val
    ]

    if 1 not in rank_pos:
        return rank_pos, False, False

    o_mods = [i for i in rank_pos[1] if not _is_pmod(i)]
    p_mods = [i for i in rank_pos[1] if _is_pmod(i)]
    psp_val_f = [i for i in psp_val if i[2] > 50]

    if len(p_mods) != len(psp_val_f):
        LOGGER.debug(
            "Not enough info to assign phophosite: {}".format(psp_val)
        )
        ambiguous = True
        rank_pos = None
    elif set(i[0] + 1 for i in p_mods) != set(i[1] for i in psp_val_f):
        p_mods = [
            (i[1] - 1, "Phospho")
            for i in psp_val_f
        ]
        reassigned = True

        rank_pos = _sort_mods(o_mods + p_mods)

    return rank_pos, reassigned, ambiguous


def _get_phosphors_psp_vals(cursor):
    fields = cursor.execute(
        """
        SELECT
        CustomDataFields.FieldID,
        CustomDataFields.DisplayName
        FROM CustomDataFields
        """,
    )
    field_ids = [
        field_id
        for field_id, name in fields
        if name in ["phosphoRS Site Probabilities"]
    ]

    if not field_ids:
        return {}

    psp_vals = cursor.execute(
        """
        SELECT
        CustomDataPeptides.PeptideID,
        CustomDataPeptides.FieldValue
        FROM CustomDataPeptides
        WHERE CustomDataPeptides.FieldID IN ({})
        """.format(
            ", ".join("?" * len(field_ids))
        ),
        field_ids,
    )
    return {
        pep_id: psp_val
        for pep_id, psp_val in psp_vals
    }


def read_discoverer_msf(msf_path):
    """
    Parse a ProteomeDiscoverer MSF file.

    Parameters
    ----------
    msf_path : str
        Path to MSF file.

    Returns
    -------
    fixed_mods : list of str
    var_mods : list of str
    out : list of :class:`PeptideQuery<pycamv.pep_query.PeptideQuery>`
    """
    LOGGER.info(
        "Loading ProteomeDiscoverer peptides from \"{}\"".format(
            os.path.basename(msf_path),
        )
    )
    start = datetime.now()

    with sqlite3.connect(msf_path) as conn:
        fixed_mods, var_mods = _get_fixed_var_mods(conn)
        out = _get_peptide_queries(conn, fixed_mods, var_mods)

    LOGGER.info(
        "Loaded {} peptides in {} hr:min:sec"
        .format(
            len(out),
            str(datetime.now() - start).split('.')[0],
        )
    )

    return fixed_mods, var_mods, out
