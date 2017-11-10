
from __future__ import absolute_import, division

import os

from .utils import nCr, fuzzy_find
from . import ms_labels


class PeptideQuery:
    """
    Attributes
    ----------
    accessions : list of str
    prot_descs : list of str
    uniprot_accessions : list of str
    full_seqs : list of str
    pep_offsets : list of int
    query : int
    filename: str
    pep_score : float
    pep_exp_mz : float
    pep_exp_z : int
    pep_seq : str
    pep_var_mods : list of tuple of (int, str, tuple of str)
    pep_fixed_mods : list of tuple of (int, str, tuple of str)
    scan : int
    quant_scan : int or None
    num_comb : int
    rank_pos : set of tuple of int, str or None
    """
    def __init__(
        self,
        query,
        filename,
        pep_score,
        pep_exp_mz, pep_exp_z,
        pep_seq,
        pep_var_mods, pep_fixed_mods, scan,
        rank_pos=None,
        quant_scan=None,
        accessions=None,
        prot_descs=None,
        uniprot_accessions=None,
        full_seqs=None,
        pep_offsets=None,
    ):
        assert _check_mods(pep_var_mods)
        assert _check_mods(pep_fixed_mods)

        self.prot_descs = prot_descs or []
        self.accessions = accessions or []
        self.uniprot_accessions = uniprot_accessions or []
        self.full_seqs = full_seqs or []
        self.pep_offsets = pep_offsets or _get_offsets(pep_seq, full_seqs)

        self.query = query
        self.filename = filename
        self.pep_score = pep_score
        self.pep_exp_mz = pep_exp_mz
        self.pep_exp_z = pep_exp_z
        self.pep_seq = pep_seq.upper()
        self.pep_var_mods = pep_var_mods
        self.pep_fixed_mods = pep_fixed_mods
        self.scan = scan
        self.quant_scan = quant_scan
        self.num_comb = self._calc_num_comb()
        self.rank_pos = rank_pos

    def _unique_tuple(self):
        return (
            tuple(self.accessions),
            self.query,
            self.pep_seq,
            self.scan,
        )

    def __hash__(self):
        return hash(self._unique_tuple())

    def __eq__(self, other):
        if not isinstance(other, PeptideQuery):
            raise TypeError(other)
        return self._unique_tuple() == other._unique_tuple()

    @property
    def basename(self):
        return os.path.basename(self.filename)

    @property
    def prot_name(self):
        return " / ".join(sorted(self.prot_descs))

    @property
    def pep_mods(self):
        return self.pep_var_mods + self.pep_fixed_mods

    @property
    def get_label_mods(self):
        return [
            mod
            for _, mod, letters in self.pep_mods
            if mod in ms_labels.LABEL_NAMES and "N-term" in letters
        ]

    def _calc_num_comb(self):
        num_comb = 1

        for count, mod, letters in self.pep_var_mods:
            if mod == "Phospho" and letters == ["S", "T"]:
                letters = ["S", "T", "Y"]

            potential_mod_sites = sum(self.pep_seq.count(i) for i in letters)

            # Subtract sites that will be taken up by another modification
            # (i.e. Oxidation and Dioxidation of M)
            for o_count, o_mod, o_letters in self.pep_var_mods:
                if o_letters == letters and o_mod != mod:
                    potential_mod_sites -= o_count

            num_comb *= nCr(potential_mod_sites, count)

        return num_comb


def _check_mods(mods):
    return all(
        isinstance(count, int) and
        isinstance(abbrev, str) and
        isinstance(letters, tuple) and
        all(isinstance(i, str) for i in letters)
        for count, abbrev, letters in mods
    )


def _get_offsets(pep_seq, full_seqs):
    def _get_rel_pos(seq):
        if not seq:
            return 0, False

        pep_pos = seq.find(pep_seq)
        exact = True

        if pep_pos < 0:
            pep_pos = fuzzy_find(pep_seq, seq)
            exact = False

        return pep_pos, exact

    return [_get_rel_pos(i)[0] for i in full_seqs]
