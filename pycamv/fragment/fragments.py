"""
This module provides functionality for calculating the masses of peptide
fragments.
"""

from __future__ import absolute_import, division

from collections import Counter

from . import masses, ms_labels

from .losses import PEPTIDE_LOSSES, INTERNAL_LOSSES, AA_LOSSES, MOD_LOSSES


DELTA_C13 = masses.exact_mass({"C": [-1, 1]})


def _sequence_mass(pep_seq):
    return sum(
        masses.AMINO_ACIDS[letter] +
        sum(
            masses.MODIFICATIONS[letter, mod]
            for mod in mods
        )
        for letter, mods in pep_seq
    )


def _sequence_name(pep_seq):
    return "".join(
        letter.lower() if mods else letter.upper()
        for letter, mods in pep_seq
        if letter not in ["N-term", "C-term", "C=O"]
    )


def _internal_fragment_ions(
    pep_seq,
    c13_num=0,
    any_losses=None, aa_losses=None, mod_losses=None,
):
    """
    Calculate the m/z of all internal fragmenets of a peptide.

    Parameters
    ----------
    pep_seq : list of tuple of (str, list of str)
        A list of peptide letters and each residue's modification(s).
    c13_num : int, optional
    aa_losses : list of str, optional
        Potential neutral losses for each fragment (i.e. Water, amine, CO).
        List is composed of neutral loss names.
    aa_losses : dict of str, str, optional
    mod_losses : dict of tuple of (str, str), list of str
        Potential neutral losses for modified amino acids (i.e. pY-HPO_3).
        Dictionary should map (letter, modification) to a list of neutral
        loss names.

    Returns
    -------
    generator of str, float
        Ion names and their corresponding m/z's.
    """
    if any_losses is None:
        any_losses = INTERNAL_LOSSES

    if aa_losses is None:
        aa_losses = AA_LOSSES

    if mod_losses is None:
        mod_losses = MOD_LOSSES

    # Calculate the mass of the peptide cut at any two sites between the N-
    # and C-terminus
    for start in range(2, len(pep_seq) - 2):
        for end in range(start + 1, len(pep_seq) - 1):
            # Only add the mass of an N-terminus, cleavage will be between
            # C=O and N-H bond, adding a hydrogen to N-H
            fragment = (
                [("N-term", ())] +
                list(pep_seq[start:end]) +
                [("C=O", ())]
            )

            mass = _sequence_mass(fragment)
            name = _sequence_name(fragment)

            # Also calculate any potential neutral losses from this fragment
            for loss_name, loss_mass in _generate_losses(
                pep_seq=fragment,
                max_depth=1,
                c13_num=c13_num,
                any_losses=any_losses,
                aa_losses=aa_losses,
                mod_losses=mod_losses,
            ):
                yield name + loss_name, mass + loss_mass


def _get_frag_masses(pep_seq):
    return [
        _sequence_mass([pep_seq[index]])
        for index in range(len(pep_seq))
    ]


def _charged_m_zs(name, mass, max_charge):
    for charge in range(1, max_charge + 1):
        yield (
            (
                name +
                (
                    "^{{{:+}}}".format(charge)
                    if charge > 1 else
                    "^{+}"
                )
            ),
            (mass + charge * masses.PROTON) / charge,
        )


def _format_loss(name, count):
    if not name or count == 0:
        return ""

    return "{}{}{}".format(
        "+" if count > 0 else "-",
        "{} ".format(abs(count)) if abs(count) > 1 else "",
        name,
    )


def _generate_c13(c13_num):
    for c13 in range(0, c13_num + 1):
        c13_mass = c13 * DELTA_C13
        c13_name = _format_loss("^{13}C", c13)
        yield c13_name, c13_mass


def _generate_losses(
    pep_seq=None,
    max_depth=2,
    c13_num=0,
    any_losses=None, aa_losses=None, mod_losses=None,
):
    def _generate_loss_combos(
        seq=None,
        losses=None, max_depth=2,
    ):
        # TODO: Check if duplicates are generated here
        if losses is None:
            losses = Counter({"": 0})
            yield losses

        if not seq:
            return

        if max_depth < 1:
            yield losses
            return

        for loss in any_losses:
            new_losses = losses.copy()

            for ion in loss.split("-"):
                new_losses[ion] -= 1

            yield new_losses

            child_losses = _generate_loss_combos(
                seq, new_losses,
                max_depth=max_depth - 1,
            )

            for new_losses in child_losses:
                yield new_losses

        for aa, a_losses in aa_losses.items():
            for index, (letter, mods) in enumerate(seq):
                if aa != letter:
                    continue

                new_seq = seq[:index] + seq[index + 1:]

                for loss in a_losses:
                    new_losses = losses.copy()

                    for ion in loss.split("-"):
                        new_losses[ion] -= 1

                    yield new_losses

                    child_losses = _generate_loss_combos(
                        new_seq, new_losses,
                        max_depth=max_depth - 1,
                    )

                    for new_losses in child_losses:
                        yield new_losses

        for (aa, mod), m_losses in mod_losses.items():
            for index, (letter, mods) in enumerate(seq):
                if aa != letter or (mod and mod not in mods):
                    continue

                new_seq = seq[:index] + seq[index + 1:]

                for loss in m_losses:
                    new_losses = losses.copy()

                    for ion in loss.split("-"):
                        new_losses[ion] -= 1

                    yield new_losses

                    child_losses = _generate_loss_combos(
                        new_seq, new_losses,
                        max_depth=max_depth - 1,
                    )

                    for new_losses in child_losses:
                        yield new_losses

    for loss in _generate_loss_combos(
        seq=pep_seq,
        max_depth=max_depth,
    ):
        loss_mass = sum(
            masses.MASSES[name] * count
            for name, count in loss.items()
            if name
        )
        loss_name = "".join(
            _format_loss(name, count)
            for name, count in sorted(loss.items(), key=lambda x: x[0])
        )

        for c13_name, c13_mass in _generate_c13(c13_num):
            yield loss_name + c13_name, loss_mass + c13_mass


def _b_y_ions(
    pep_seq, frag_masses, fragment_max_charge,
    c13_num=0,
    any_losses=None, aa_losses=None, mod_losses=None,
):
    def _generate_ions(seq, mass, basename):
        # b/y ions +/- losses
        for loss_name, loss_mass in _generate_losses(
            pep_seq=seq,
            c13_num=c13_num,
            any_losses=any_losses,
            aa_losses=aa_losses,
            mod_losses=mod_losses,
        ):
            for name, mz in _charged_m_zs(
                basename + loss_name,
                mass + loss_mass,
                fragment_max_charge,
            ):
                yield name, mz

    base_ions = {}

    for index in range(2, len(pep_seq) - 1):
        # XXX: iTRAQ / TMT y-adducts?
        base_ions["a_{{{}}}".format(index - 1)] = (
            sum(frag_masses[:index]) - masses.MASSES["CO"],
            pep_seq[:index],
        )
        base_ions["b_{{{}}}".format(index - 1)] = (
            sum(frag_masses[:index]),
            pep_seq[:index],
        )

    for index in range(1, len(pep_seq) - 1):
        # y ion: 1 hydrogen added to NH group, one hydrogen on K/R
        base_ions["y_{{{}}}".format(len(pep_seq) - index - 1)] = (
            sum(frag_masses[index:]) + masses.PROTON,
            pep_seq[index:],
        )

    for ion_name, (ion_mass, seq) in base_ions.items():
        for name, mz in _generate_ions(
            seq, ion_mass, ion_name,
        ):
            yield name, mz


def _label_ions(pep_seq):
    label_mods = [
        mod
        for mod in pep_seq[0][1]
        if mod in ms_labels.LABEL_NAMES
    ]

    for mod in label_mods:
        for name, mz in zip(
            ms_labels.LABEL_NAMES[mod],
            ms_labels.LABEL_MASSES[mod],
        ):
            yield name, mz


def _parent_ions(
    pep_seq, frag_masses, parent_max_charge,
    c13_num=0,
    any_losses=None, aa_losses=None, mod_losses=None,
):
    parent_mass = sum(frag_masses) + masses.PROTON

    for loss_name, loss_mass in _generate_losses(
        pep_seq=pep_seq,
        c13_num=c13_num,
        any_losses=any_losses,
        aa_losses=aa_losses,
        mod_losses=mod_losses,
    ):
        for name, mz in _charged_m_zs(
            "MH" + loss_name,
            parent_mass + loss_mass,
            parent_max_charge,
        ):
            yield name, mz


def _py_ions(pep_seq, c13_num=0):
    if any(
        letter == "Y" and "Phospho" in mods
        for letter, mods in pep_seq
    ):
        name = "pY"
        mass = (
            masses.IMMONIUM_IONS["Y"] +
            masses.MODIFICATIONS["Y", "Phospho"]
        )

        for loss_name, loss_mass in _generate_losses(
            c13_num=c13_num,
        ):
            yield name + loss_name, mass + loss_mass


def fragment_ions(
    pep_seq, charge,
    parent_max_charge=None, fragment_max_charge=None,
    c13_num=0,
    any_losses=None, aa_losses=None, mod_losses=None,
):
    """
    Calculate the m/z of all ions generated by a peptide.

    Parameters
    ----------
    pep_seq : list of tuple of (str, list of str)
    charge : int
    parent_max_charge : int, optional
    fragment_max_charge : int, optional
    c13_num : int, optional
    any_losses : list of str, optional
        Potential neutral losses for each fragment (i.e. Water, amine, CO).
        List is composed of neutral loss names.
    aa_losses : dict of str, str, optional
    mod_losses : dict of tuple of (str, str), str, optional
        Potential neutral losses for modified amino acids (i.e. pY-HPO_3).
        Dictionary should map (letter, modification) to a list of neutral
        loss names.

    Returns
    -------
    dict of int, dict of str, float
        Dictionary mapping fragment position to a dictionary mapping ion names
        to ion m/z's.
    """
    assert "N-term" == pep_seq[0][0]
    assert "C-term" == pep_seq[-1][0]

    if parent_max_charge is None:
        parent_max_charge = charge

    if fragment_max_charge is None:
        # TODO: This correct?
        fragment_max_charge = parent_max_charge - 1

    if any_losses is None:
        any_losses = PEPTIDE_LOSSES

    if aa_losses is None:
        aa_losses = AA_LOSSES

    if mod_losses is None:
        mod_losses = MOD_LOSSES

    # First calculate the masses of each residue along the backbone
    frag_masses = _get_frag_masses(pep_seq)

    frag_ions = {}

    # Get b/y (and associated a/c/x/z) ions
    frag_ions.update(
        _b_y_ions(
            pep_seq, frag_masses, fragment_max_charge,
            c13_num=c13_num,
            any_losses=any_losses,
            aa_losses=aa_losses,
            mod_losses=mod_losses,
        )
    )

    # Get parent ions (i.e. MH^{+1})
    frag_ions.update(
        _parent_ions(
            pep_seq, frag_masses, parent_max_charge,
            c13_num=c13_num,
            any_losses=any_losses,
            aa_losses=aa_losses,
            mod_losses=mod_losses,
        )
    )

    # Get TMT / iTRAQ labels
    frag_ions.update(
        _label_ions(pep_seq)
    )

    # Get pY peak
    frag_ions.update(
        _py_ions(
            pep_seq,
            c13_num=c13_num,
        )
    )

    # Get internal fragments
    frag_ions.update(
        _internal_fragment_ions(
            pep_seq,
            c13_num=c13_num,
            aa_losses=aa_losses,
            mod_losses=mod_losses,
        )
    )

    return frag_ions
