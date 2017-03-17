PEPTIDE_LOSSES = [
    "H_2O",
    "NH_3",
]
INTERNAL_LOSSES = [
    "H_2O",
    "NH_3",
    "CO",
]
AA_LOSSES = {}
MOD_LOSSES = {
    ("M", "Oxidation"): ["SOCH_4"],
    ("M", "Dioxidation"): ["SO_2CH_4"],
    ("S", "Phospho"): ["H_3PO_4"],
    ("T", "Phospho"): ["H_3PO_4"],
    ("Y", "Phospho"): ["HPO_3", "HPO_3-H_2O"],
}
