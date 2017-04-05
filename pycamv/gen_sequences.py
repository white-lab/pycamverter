"""
This module provides functionality for generating all possible combinations of
modifications on a peptide.
"""

from __future__ import absolute_import, division

import itertools


MAX_NUM_COMB = 10


def gen_possible_seq(pep_seq, var_mods):
    """
    Generates sequences

    Extra rules: Acetlyation cannot happen on terminal lysines, as trypsin does
    not cut at these sites.

    Parameters
    ----------
    pep_seq : str
    var_mods : list of tuple of (int, str, list of str)

    Returns
    -------
    generator of list of tuple of (str, list of str)
    """
    def _gen_mods(seq, mods):
        if not mods:
            yield seq
            return

        count, mod, letters = mods[0]

        # XXX: Allow multiple modifications on a peptide
        indices = [
            index
            for index, (letter, mods) in enumerate(seq)
            if letter in letters and not mods
        ]

        if len(indices) < count:
            raise Exception(
                (
                    "Too few sites for modification \"{} {} ({})\" "
                    "in peptide \"{}\""
                ).format(count, mod, "".join(letters), "".join(pep_seq[1:-1]))
            )

        for mod_indices in itertools.combinations(indices, count):
            tmp_seq = [
                (
                    letter,
                    (mods + (mod,)) if index in mod_indices else mods,
                )
                for index, (letter, mods) in enumerate(seq)
            ]

            for new_seq in _gen_mods(tmp_seq, mods[1:]):
                yield new_seq

    tmp_seq = ["N-term"] + list(pep_seq) + ["C-term"]

    # Sort var mods so that we assign more restricted sites first
    # (i.e. pY before pSTY)
    var_mods = sorted(var_mods, key=lambda x: len(x[2]))

    # Ideally we would use py>3.3 syntax: yield from, but we are also
    # supporting py2.7
    for seq in _gen_mods(
        seq=list(zip(tmp_seq, [() for _ in tmp_seq])),
        mods=var_mods,
    ):
        yield seq
