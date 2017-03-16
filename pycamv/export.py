"""
Provides functionality for exporting processed data to .camv files (JSON
format).
"""

from __future__ import absolute_import, division

from collections import OrderedDict
import gzip
import logging
import json

from .utils import DefaultOrderedDict
from . import ms_labels, regexes, version


LOGGER = logging.getLogger("pycamv.export")

SUPERSCRIPT_UNICODE_START = ord(u"\u2070")
SUBSCRIPT_UNICODE_START = ord(u'\u2080')
SCRIPT_MAPPING = {
    str(i): i
    for i in range(10)
}
SCRIPT_MAPPING["+"] = 10
SCRIPT_MAPPING["-"] = 11
SCRIPT_MAPPING["("] = 12
SCRIPT_MAPPING[")"] = 13


def _rewrite_ion_name(name):
    m = regexes.RE_B_Y_IONS.match(name)

    if m:
        name = "".join(m.group(1, 2))

    ret = ""
    sup, sub = False, False
    paren = 0

    for char in name:
        if char in "^_{}":
            if char == "^":
                sup, sub = True, False
            elif char == "_":
                sup, sub = False, True
            elif char == "}":
                sup, sub = False, False
                paren -= 1
            elif char == "{":
                paren += 1
            continue

        if sup:
            if char == "1":
                ret += u"\u00B9"
            elif char == "2":
                ret += u"\u00B2"
            elif char == "3":
                ret += u"\u00B3"
            else:
                ret += chr(SUPERSCRIPT_UNICODE_START + SCRIPT_MAPPING[char])
        elif sub:
            ret += chr(SUBSCRIPT_UNICODE_START + SCRIPT_MAPPING[char])
        else:
            ret += char

        if sup or sub:
            if not paren:
                sup, sub = False, False

    return ret


def _peaks_to_dict(peaks):
    return [
        OrderedDict([
            ("mz", mz),
            ("into", i),
        ])
        for mz, i in peaks
    ]


def _join_seq_mods(seq, mods):
    return tuple(zip(["N-term"] + list(seq) + ["C-term"], mods))


def _extract_pep_seq(sequence):
    return "".join(
        letter
        for letter, _ in sequence
        if letter not in ["N-term", "C-term"]
    )


def _extract_mods(sequence):
    return tuple(
        tuple(mods)
        for _, mods in sequence
    )


def _get_mods_description(pep_seq, mod_state):
    return ("+ " if mod_state else "") + " - ".join(
        "{} {}{}".format(count, mod[0].lower(), "".join(letters))
        for count, mod, letters in mod_state
    )


def _mod_positions(mods):
    return [
        index
        for index, mod in enumerate(mods)
        if mod
    ]


def _pep_mod_name(pep_seq, mods):
    return "".join(
        letter.lower() if mods else letter.upper()
        for letter, mods in zip(pep_seq, mods[1:-1])
    )


def _get_labels_mz(query):
    def _within(val, bounds):
        return val >= bounds[0] and val <= bounds[1]
    return [
        mz
        for mod in set(query.get_label_mods)
        for mz in ms_labels.LABEL_MASSES.get(mod, [])
        if _within(mz, ms_labels.LABEL_MZ_WINDOW.get(mod, [mz, mz]))
    ]


def export_to_camv(
    out_path, peak_hits, scan_mapping,
    precursor_windows, label_windows,
):
    """
    Export data to .camv file.

    Parameters
    ----------
    out_path : str
    peak_hits
    scan_mapping
    precursor_windows
    label_windows

    Returns
    -------
    data : dict
        Dictionary of data written to file.

    Example
    -------
    >>> from pycamv import export, validate
    >> name = "2016-06-27-CKX13-pY-NTA-elute-pre67-col77"
    >>> options, peak_hits, scan_mapping, precursor_windows, label_windows = \
    ...     validate.validate_spectra(name)
    >>> json_out = export.export_to_camv(
    ...     "{}.camv".format(name),
    ...     peak_hits, scan_mapping, precursor_windows, label_windows,
    ... )
    """
    ###
    # Mappings between proteins, peptides, modifications, queries, and peaks
    ###
    # Mapping for protein -> peptides
    prot_dict = DefaultOrderedDict(set)

    for query, _ in peak_hits.keys():
        prot_dict[query.prot_name].add(query.pep_seq)

    pep_dict = DefaultOrderedDict(set)

    for query, _ in peak_hits.keys():
        # query.pep_var_mods => "mod_state" / "+1 pSTY"
        pep_dict[query.prot_name, query.pep_seq].add(tuple(query.pep_var_mods))

    mod_states_dict = DefaultOrderedDict(set)

    for query, seq in peak_hits.keys():
        mod_states_dict[
            query.prot_name, query.pep_seq, tuple(query.pep_var_mods)
        ].add(_extract_mods(seq))

    mod_query_dict = DefaultOrderedDict(set)

    for query, seq in peak_hits.keys():
        mod_query_dict[
            query.prot_name, query.pep_seq, tuple(query.pep_var_mods)
        ].add(query)

    # Mapping for modifications -> queries
    exact_mods_dict = DefaultOrderedDict(list)

    for query, seq in peak_hits.keys():
        exact_mods_dict[
            query.prot_name, query.pep_seq, _extract_mods(seq)
        ].append(query)

    # Mapping back from queries -> sequence + modifications
    query_dict = DefaultOrderedDict(list)

    for query, seq in peak_hits.keys():
        # _extract_mods(seq) => "exact_mods" / "pY100, oxM102"
        query_dict[query].append(
            (
                query.prot_name, query.pep_seq,
                tuple(query.pep_var_mods), _extract_mods(seq),
            )
        )

    ###
    # Pre-calculate IDs for later reference
    ###
    # Protein IDs
    prot_index = {
        prot_name: index
        for index, prot_name in enumerate(sorted(prot_dict.keys()))
    }

    # Peptide IDs
    pep_index = {
        (prot_name, pep_seq, mod_state): index
        for prot_name, peptides in prot_dict.items()
        for index, (pep_seq, mod_state) in enumerate(
            (pep_seq, mod_state)
            for pep_seq in peptides
            for mod_state in pep_dict[prot_name, pep_seq]
        )
    }

    # Peptide Data IDs
    pep_data_index = {
        (prot_name, pep_seq): index
        for index, (prot_name, pep_seq) in enumerate(
            sorted(
                (prot_name, pep_seq)
                for prot_name, peptides in prot_dict.items()
                for pep_seq in peptides
            )
        )
    }

    # Mod State IDs
    mod_state_index = {
        (prot_name, pep_seq, mod_state): index
        for (prot_name, pep_seq), mod_states in pep_dict.items()
        for index, mod_state in enumerate(sorted(mod_states))
    }

    # Modification IDs
    exact_mods_index = {
        (prot_name, pep_seq, mod_state, exact_mods): index
        for (prot_name, pep_seq, mod_state), exact_mods_list
        in mod_states_dict.items()
        for index, exact_mods in enumerate(sorted(exact_mods_list))
    }

    # Scan IDs
    scan_index = {
        query: index
        for (prot_name, pep_seq), mod_states in pep_dict.items()
        for mod_state in mod_states
        for index, query in enumerate(
            sorted(
                mod_query_dict[prot_name, pep_seq, mod_state],
                key=lambda x: x.scan,
            )
        )
    }

    # Peak Match IDs
    match_index = {
        (prot_name, pep_seq, exact_mods, name): index
        for (prot_name, pep_seq, exact_mods), queries
        in exact_mods_dict.items()
        for index, name in enumerate(
            name
            for name in sorted(
                set(
                    name
                    for query in queries
                    for peak_hit in peak_hits[
                        query, _join_seq_mods(pep_seq, exact_mods)
                    ]
                    if peak_hit.match_list
                    for name in peak_hit.match_list.keys()
                )
            )
        )
    }

    ###
    # Individual data parsing functions
    ###
    def _gen_match_data(prot_name, pep_seq, exact_mods):
        """
        Generate a list of all potential peak matches
        """
        visited = set()
        queries = exact_mods_dict[prot_name, pep_seq, exact_mods]

        for query in queries:
            seq_mods = _join_seq_mods(pep_seq, exact_mods)

            for peak_hit in peak_hits[query, seq_mods]:
                if not peak_hit.match_list:
                    continue

                for name, (mz, _) in peak_hit.match_list.items():
                    if name in visited:
                        continue

                    name_match = regexes.RE_BY_ION_POS.match(name)

                    if name_match:
                        ion_pos = int(name_match.group(2))
                        ion_type = "b" if name_match.group(1) in "abc" else "y"
                    else:
                        ion_type, ion_pos = None, None

                    yield OrderedDict([
                        (
                            "ionId",
                            match_index[
                                query.prot_name, query.pep_seq,
                                exact_mods, name,
                            ]
                        ),
                        ("mz", mz),
                        ("name", _rewrite_ion_name(name)),
                        ("ionType", ion_type),
                        ("ionPosition", ion_pos),
                    ])

                    visited.add(name)

    def _get_match_data(prot_name, pep_seq, exact_mods):
        return sorted(
            _gen_match_data(prot_name, pep_seq, exact_mods),
            key=lambda x: x["ionId"],
        )

    def _get_exact_mod_data(prot_name, pep_seq, mod_state):
        return [
            OrderedDict([
                (
                    "exactModsId",
                    exact_mods_index[
                        prot_name, pep_seq, mod_state, exact_mods
                    ],
                ),
                ("positions", _mod_positions(exact_mods)),
                ("name", _pep_mod_name(pep_seq, exact_mods)),
                ("matchData", _get_match_data(prot_name, pep_seq, exact_mods)),
            ])
            for exact_mods in sorted(
                mod_states_dict[prot_name, pep_seq, mod_state],
                key=lambda x: exact_mods_index[
                    prot_name, pep_seq, mod_state, x
                ],
            )
        ]

    def _get_mod_states(prot_name, pep_seq, mod_states):
        return [
            OrderedDict([
                ("modStateId", mod_state_index[prot_name, pep_seq, mod_state]),
                ("modDesc", _get_mods_description(pep_seq, mod_state)),
                (
                    "mods",
                    _get_exact_mod_data(
                        prot_name,
                        pep_seq,
                        mod_state,
                    ),
                ),
            ])
            for mod_state in sorted(
                mod_states,
                key=lambda x: mod_state_index[prot_name, pep_seq, x],
            )
        ]

    def _get_peptide_data():
        """
        Return all information mapping (modified) peptides to their sequence,
        descriptions, ion fragmentation patterns, etc.
        """
        return [
            OrderedDict([
                ("proteinId", index),
                ("proteinName", prot_name),
                ("peptideSequence", pep_seq),
                (
                    "modificationStates",
                    _get_mod_states(
                        prot_name, pep_seq, pep_dict[prot_name, pep_seq],
                    ),
                ),
            ])
            for (prot_name, pep_seq), index in sorted(
                pep_data_index.items(),
                key=lambda x: x[1],
            )
        ]

    def _get_default_choice_data(prot_name, pep_seq, mod_state):
        return [
            OrderedDict([
                (
                    "modsId",
                    exact_mods_index[
                        prot_name, pep_seq, mod_state, exact_mods
                    ]),
                ("state", None),  # null
            ])
            for exact_mods in sorted(
                mod_states_dict[prot_name, pep_seq, mod_state],
                key=lambda x: exact_mods_index[
                    prot_name, pep_seq, mod_state, x
                ],
            )
        ]

    def _get_matches(prot_name, peak_index, query):
        return [
            OrderedDict([
                (
                    "modsId",
                    exact_mods_index[prot_name, seq, mod_state, exact_mods],
                ),
                (
                    "matchId",
                    match_index.get(
                        (
                            prot_name,
                            seq,
                            exact_mods,
                            peak_hits[
                                query,
                                _join_seq_mods(seq, exact_mods),
                            ][peak_index].name,
                        ),
                        None,
                    ),
                ),
            ])
            for _, seq, mod_state, exact_mods in sorted(
                query_dict[query],
                key=lambda x: exact_mods_index[x],
            )
        ]

    def _get_scan_assignments(query, seq):
        exact_mods = query_dict[query][0][-1]

        return [
            OrderedDict([
                ("mz", peak_hit.mz),
                ("into", peak_hit.intensity),
                (
                    "matchInfo",
                    _get_matches(query.prot_name, peak_index, query)
                ),
            ])
            for peak_index, peak_hit in sorted(
                enumerate(peak_hits[query, _join_seq_mods(seq, exact_mods)]),
                key=lambda x: x[1].mz,
            )
        ]

    def _get_scans(prot_name, pep_seq, mod_state):
        """
        Return information on individual scans, including peaks, precursor
        ions, and peptide modification assignments.
        """
        return [
            OrderedDict([
                ("scanId", scan_index[query]),
                ("fileName", query.filename),
                ("scanNumber", query.scan),
                ("chargeState", query.pep_exp_z),
                (
                    "precursorScanData",
                    _peaks_to_dict(precursor_windows[query]),
                ),
                ("precursorMz", query.pep_exp_mz),
                (
                    "precursorIsolationWindow",
                    scan_mapping[query].window_offset,
                ),
                ("collisionType", scan_mapping[query].collision_type),
                ("c13Num", scan_mapping[query].c13_num),
                ("quantScanData", _peaks_to_dict(label_windows[query])),
                ("quantMz", _get_labels_mz(query)),
                (
                    "choiceData",
                    _get_default_choice_data(prot_name, pep_seq, mod_state),
                ),
                (
                    "scanData",
                    _get_scan_assignments(query, pep_seq),
                ),
            ])
            for query in sorted(
                mod_query_dict[prot_name, pep_seq, mod_state],
                key=lambda x: scan_index[x],
            )
        ]

    def _get_peptide_scan_data(prot_name, peptides):
        """
        Map peptides to their data IDs, scans, and candidate modification
        states.
        """
        return [
            OrderedDict([
                ("peptideId", pep_index[prot_name, pep_seq, mod_state]),
                ("peptideDataId", pep_data_index[prot_name, pep_seq]),
                (
                    "modificationStateId",
                    mod_state_index[prot_name, pep_seq, mod_state],
                ),
                ("scans", _get_scans(prot_name, pep_seq, mod_state)),
            ])
            for _, pep_seq, mod_state in sorted(
                (
                    (prot_name, pep_seq, mod_state)
                    for pep_seq in peptides
                    for mod_state in pep_dict[prot_name, pep_seq]
                ),
                key=lambda x: pep_index[x],
            )
        ]

    def _get_scan_data():
        """
        Return all information mapping proteins / peptides to their scans and
        a list of candidate modification patterns.
        """
        return [
            OrderedDict([
                ("proteinId", prot_index[prot_name]),
                ("proteinName", prot_name),
                ("peptides", _get_peptide_scan_data(prot_name, peptides)),
            ])
            for prot_name, peptides in sorted(
                prot_dict.items(),
                key=lambda x: prot_index[x[0]],
            )
        ]

    data = OrderedDict([
        ("pycamverterVersion", version.__version__),
        ("peptideData", _get_peptide_data()),
        ("scanData", _get_scan_data()),
    ])

    LOGGER.info("Exporting CAMV data to {}".format(out_path))

    if out_path.endswith(".gz"):
        with gzip.GzipFile(out_path, 'w') as f:
            f.write(json.dumps(data, indent=None).encode("utf-8"))
    else:
        with open(out_path, "w") as f:
            json.dump(data, f, indent=None)

    return data
