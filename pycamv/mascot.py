"""
Provides functionality for interacting with MASCOT data.
"""

from __future__ import absolute_import, division

import logging
import xml.etree.ElementTree as ET

from . import pep_query, regexes


MASCOT_NS = {
    "mascot":
        "http://www.matrixscience.com/xmlns/schema/mascot_search_results_2",
}

LOGGER = logging.getLogger("pycamv.mascot")


def _count_instances(pep_seq, letters):
    return sum(
        (["N-term"] + list(pep_seq) + ["C-term"]).count(letter)
        for letter in letters
    )


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


def _get_fixed_var_mods(root):
    fixed_mods = [
        i.text
        for i in root.findall(
            "mascot:fixed_mods/mascot:modification/mascot:name",
            MASCOT_NS,
        )
    ]
    var_mods = [
        i.text
        for i in root.findall(
            "mascot:variable_mods/mascot:modification/mascot:name",
            MASCOT_NS,
        )
    ]
    return fixed_mods, var_mods


def _parse_mascot_2_4_1(root):
    filename = root.find("mascot:header/mascot:FILENAME", MASCOT_NS).text
    filename = filename.split(":", 1)[1].strip()
    fixed_mods, var_mods = _get_fixed_var_mods(root)

    scan_used = {}
    index = 0
    out = []

    for hit in root.findall("mascot:hits/mascot:hit", MASCOT_NS):
        accession = hit.find("mascot:protein", MASCOT_NS).get("accession")
        prot_desc = hit.find("mascot:protein/mascot:prot_desc", MASCOT_NS).text

        if "FGR" in prot_desc.upper() or "FGR" in accession.upper():
            LOGGER.info(accession, prot_desc)

        match = regexes.RE_MASCOT_DESCRIPTION.match(prot_desc)

        if match:
            prot_desc = match.group(1)

        for peptide in hit.findall("mascot:protein/mascot:peptide", MASCOT_NS):
            query = int(peptide.get("query"))
            rank = int(peptide.get("rank"))
            pep_score = float(peptide.find("mascot:pep_score", MASCOT_NS).text)
            exp_mz = float(peptide.find("mascot:pep_exp_mz", MASCOT_NS).text)
            exp_z = int(peptide.find("mascot:pep_exp_z", MASCOT_NS).text)
            pep_seq = peptide.find("mascot:pep_seq", MASCOT_NS).text

            pep_var_mods = peptide.find("mascot:pep_var_mod", MASCOT_NS).text

            if fixed_mods:
                pep_fixed_mods = [
                    regexes.RE_DYN_MODS.match(mod.strip()).group(3, 4)
                    for mod in fixed_mods
                ]
                pep_fixed_mods = [
                    (name, _parse_letters(letters))
                    for name, letters in pep_fixed_mods
                ]
                pep_fixed_mods = [
                    (_count_instances(pep_seq, letters), name, letters)
                    for name, letters in pep_fixed_mods
                ]
                pep_fixed_mods = [
                    (count, name, letters)
                    for count, name, letters in pep_fixed_mods
                    if count > 0
                ]
            else:
                pep_fixed_mods = []

            if pep_var_mods:
                # i.e. "2 Phospho (STY)""
                pep_var_mods = [
                    regexes.RE_DYN_MODS.match(mod.strip()).group(2, 3, 4)
                    for mod in pep_var_mods.split(";")
                ]
                pep_var_mods = [
                    (int(count) if count else 1, name, _parse_letters(letters))
                    for count, name, letters in pep_var_mods
                ]
            else:
                pep_var_mods = []

            scan = int(
                regexes.RE_SCAN_NUM.search(
                    peptide.find("mascot:pep_scan_title", MASCOT_NS).text
                ).group(2)
            )

            # scan_used
            if scan in scan_used:
                index_match, rank_match = scan_used[scan]

                if rank >= rank_match:
                    continue
                else:
                    del out[index_match]

            scan_used[scan] = (index, rank)
            out.append(
                pep_query.PeptideQuery(
                    accessions=[accession],
                    prot_descs=[prot_desc],
                    query=query,
                    filename=filename,
                    pep_score=pep_score,
                    pep_exp_mz=exp_mz,
                    pep_exp_z=exp_z,
                    pep_seq=pep_seq,
                    pep_var_mods=pep_var_mods,
                    pep_fixed_mods=pep_fixed_mods,
                    scan=scan,
                    quant_scan=scan,
                )
            )
            index += 1

    out = sorted(out, key=lambda x: x.scan)

    return fixed_mods, var_mods, out


def read_mascot_xml(xml_path):
    """
    Parse a MASCOT XML file.

    Parameters
    ----------
    xml_path : str
        Path to XML file.

    Returns
    -------
    fixed_mods : list of str
    var_mods : list of str
    out : list of :class:`PeptideQuery<pycamv.pep_query.PeptideQuery>`
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    ms_version = root.find("mascot:header/mascot:MascotVer", MASCOT_NS).text
    ms_version = tuple(int(i) for i in ms_version.split("."))

    if ms_version >= (2, 4, 1):
        return _parse_mascot_2_4_1(root)
#     elif ms_version == "2.4.0":
#         return _parse_mascot_2_4_0(root)
#     elif ms_version == "2.3.02":
#         return _parse_mascot_2_3_02(root)
#     elif ms_version == "2.1.03":
#         return _parse_mascot_2_1_03(root)
    else:
        raise Exception("Unsupported Mascot Version: {}".format(ms_version))
