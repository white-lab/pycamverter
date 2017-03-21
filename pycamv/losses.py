PEPTIDE_LOSSES = []
INTERNAL_LOSSES = []
AA_LOSSES = {}
MOD_LOSSES = {
    # Losses from modified amino acids
    ("M", "Oxidation"): ["SOCH_4"],
    ("M", "Dioxidation"): ["SO_2CH_4"],
    ("S", "Phospho"): ["H_3PO_4"],
    ("T", "Phospho"): ["H_3PO_4"],
    ("Y", "Phospho"): ["HPO_3", "HPO_3-H_2O"],

    # Losses from unprotected N-/C-/internal fragment termini
    ("N-term", None): ["NH_3"],
    ("C-term", None): ["H_2O"],
    ("C=O", None): ["CO"],

    # Water losses from hydroxyl / carboxyl groups
    ("S", None): ["H_2O"],
    ("T", None): ["H_2O"],
    ("E", None): ["H_2O"],
    ("D", None): ["H_2O"],

    # Amine losses from amine groups
    ("R", None): ["NH_3"],
    ("K", None): ["NH_3"],
    ("Q", None): ["NH_3"],
    ("N", None): ["NH_3"],
}
