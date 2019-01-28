from unittest import TestCase

from pycamv.fragment import fragments


class FragmentIonsTest(TestCase):
    def test_fragment_ions(self):
        test_peptides = {
            (
                (
                    ("N-term", ()),
                    ("R", ()),
                    ("V", ()),
                    ("D", ()),
                    ("E", ()),
                    ("N", ()),
                    ("N", ()),
                    ("P", ()),
                    ("E", ()),
                    ("Y", ()),
                    ("C-term", ()),
                ),
                2,
            ): {
                "a_{1}^{+}": 129.11347,
                "a_{2}^{+}": 228.18189,
                "b_{1}^{+}": 157.10839,
                "b_{2}^{+}": 256.17680,
                "b_{3}^{+}": 371.20374,
                "b_{4}^{+}": 500.24634,
                "b_{5}^{+}": 614.28926,
                "b_{6}^{+}": 728.33219,
                "b_{6}^{+2}": 364.66973,
                "b_{7}^{+}": 825.38496,
                "b_{8}^{+}": 954.42755,
                "y_{1}^{+}": 182.08117,
                "y_{2}^{+}": 311.12376,
                "y_{3}^{+}": 408.17653,
                "y_{4}^{+}": 522.21945,
                "y_{5}^{+}": 636.26238,
                "y_{6}^{+}": 765.30497,
                "y_{7}^{+}": 880.33192,
                "y_{8}^{+}": 979.40033,
            },
            (
                (
                    ("N-term", ("TMT6plex",)),
                    ("S", ()),
                    ("V", ()),
                    ("Y", ("Phospho",)),
                    ("T", ()),
                    ("E", ()),
                    ("I", ()),
                    ("K", ("TMT6plex",)),
                    ("C-term", ()),
                ),
                2,
            ): {
                "a_{1}^{+}": 289.21,
                "a_{2}^{+}": 388.28,
                "b_{1}^{+}": 317.20,
                "b_{2}^{+}": 416.27,
                "b_{3}^{+}": 659.30,
                "b_{4}^{+}": 760.35,
                "b_{5}^{+}": 889.39,
                "b_{6}^{+}": 1002.47,
                "y_{1}^{+}": 376.28,
                "y_{2}^{+}": 489.36,
                "y_{3}^{+}": 618.40,
                "y_{4}^{+}": 719.45,
                "y_{5}^{+}": 962.48,
                "y_{6}^{+}": 1061.55,
                "y_{7}^{+}": 1148.58,
                "b_{3}-HPO_3^{+}": 579.33,
                "b_{4}-HPO_3^{+}": 680.38,
                "b_{5}-HPO_3^{+}": 809.42,
                "b_{6}-HPO_3^{+}": 922.51,
                "y_{5}-HPO_3^{+}": 882.51,
                "y_{6}-HPO_3^{+}": 981.58,
                "y_{7}-HPO_3^{+}": 1068.61,
                "pY": 216.04,
                "126": 126.13,
                "127": 127.13,
                "128": 128.13,
                "129": 129.13,
                "130": 130.14,
                "131": 131.14,
                "MH^{+2}": 689.37,
                "MH-HPO_3^{+2}": 649.39,
                "MH-H_2O^{+2}": 680.37,
                "MH-HPO_3-H_2O^{+2}": 640.39,
            },
            (
                (
                    ("N-term", ("TMT6plex",)),
                    ("V", ()),
                    ("I", ()),
                    ("Y", ("Phospho",)),
                    ("D", ()),
                    ("F", ()),
                    ("I", ()),
                    ("E", ()),
                    ("K", ("TMT6plex",)),
                    ("C-term", ()),
                ),
                2,
            ): {
                "IyD": 472.15,
                "IyD-CO": 444.16,
                "IyD-NH_3": 455.12,
                "IyD-H_2O": 454.14,
                "IyD-HPO_3": 392.18,
            },
            (
                (
                    ("N-term", ("TMT6plex",)),
                    ("L", ()),
                    ("F", ()),
                    ("V", ()),
                    ("T", ("Phospho",)),
                    ("P", ()),
                    ("P", ()),
                    ("E", ()),
                    ("G", ()),
                    ("S", ()),
                    ("A", ()),
                    ("R", ()),
                    ("C-term", ()),
                ),
                4,
            ):
            {
                "a_{1}^{+}": 315.26,
                "a_{2}^{+}": 462.33,
                "a_{3}^{+}": 561.40,
                "b_{2}^{+}": 490.32,
                "b_{3}^{+}": 589.39,
                "y_{1}^{+}": 175.12,
                "y_{1}^{+}": 175.12,
                "y_{1}^{+}": 175.12,
                "y_{4}^{+}": 390.21,
                "y_{5}^{+}": 519.25,
                "y_{6}^{+}": 616.31,
                "y_{7}^{+}": 713.36,
                "y_{6}^{+2}": 308.65,
                "y_{7}^{+2}": 357.18,
                "a_{4}-H_3PO_4^{+2}": 322.72,
                "a_{5}-H_3PO_4^{+2}": 371.24,
                "b_{4}-H_3PO_4^{+}": 672.43,
                "PE": 227.10,
            },
            (
                (
                    ("N-term", ("TMT6plex",)),
                    ("D", ()),
                    ("E", ()),
                    ("I", ()),
                    ("L", ()),
                    ("P", ()),
                    ("T", ()),
                    ("T", ("Phospho",)),
                    ("P", ()),
                    ("I", ()),
                    ("S", ()),
                    ("E", ()),
                    ("Q", ()),
                    ("K", ("TMT6plex",)),
                    ("C-term", ()),
                ),
                3,
            ): {
                "PT": 199.11,
                "LP": 211.14,
                "PT-H_2O": 181.10,
                "LP-CO": 183.15,  # "LP-28"
                "PIS": 298.18,
                "PTt": 380.12,
            },
            (
                (
                    ("N-term", ("TMT6plex",)),
                    ("I", ()),
                    ("G", ()),
                    ("E", ()),
                    ("G", ()),
                    ("T", ()),
                    ("Y", ("Phospho",)),
                    ("G", ()),
                    ("V", ()),
                    ("V", ()),
                    ("Y", ()),
                    ("K", ("TMT6plex",)),
                    ("C-term", ()),
                ),
                3,
            ): {
                "a_{1}^{+}": 315.26,
                "a_{3}^{+}": 501.32,
                "a_{5}^{+}": 659.39,
                "a_{7}^{+}": 959.44,
                "b_{1}^{+}": 343.25,
                "b_{2}^{+}": 400.28,
                "b_{3}^{+}": 529.32,
                "b_{4}^{+}": 586.34,
                "b_{5}^{+}": 687.39,
                "b_{5}-H_2O^{+}": 669.38,
                "b_{6}^{+}": 930.42,
                "b_{7}^{+}": 987.44,
                "b_{8}^{+}": 1086.51,
                "y_{1}^{+}": 376.28,
                "y_{2}^{+}": 539.34,
                "y_{3}^{+}": 638.41,
                "y_{4}^{+}": 737.48,
                "y_{5}^{+}": 794.50,
                "pY": 216.04,
                "126": 126.13,
                "127": 127.13,
                "128": 128.13,
                "129": 129.13,
                "130": 130.14,
                "131": 131.14,
                "MH^{+3}": 575.31,
                "MH-HPO_3^{+3}": 548.66,
            }
        }

        for (pep_seq, charge), hits in test_peptides.items():
            frag_ions = fragments.fragment_ions(
                pep_seq, charge,
            )

            for name, mz in hits.items():
                self.assertIn(name, frag_ions)
                self.assertLess(abs(frag_ions[name] - mz), 0.01)

    def test_no_by_losses(self):
        frag_ions = fragments.fragment_ions(
            (
                ("N-term", ("TMT6plex",)),
                ("S", ()),
                ("V", ()),
                ("Y", ("Phospho",)),
                ("T", ()),
                ("E", ()),
                ("I", ()),
                ("K", ("TMT6plex",)),
                ("C-term", ()),
            ),
            2,
        )

        for ion in [
            "b_{1}-HPO_3^{+}",
            "b_{2}-HPO_3^{+}",
            "y_{1}-HPO_3^{+}",
            "y_{2}-HPO_3^{+}",
            "y_{3}-HPO_3^{+}",
            "y_{4}-HPO_3^{+}",
        ]:
            self.assertNotIn(ion, frag_ions)

    def test_c13(self):
        frag_ions = fragments.fragment_ions(
            (
                ("N-term", ("TMT6plex",)),
                ("S", ()),
                ("V", ()),
                ("Y", ("Phospho",)),
                ("T", ()),
                ("E", ()),
                ("I", ()),
                ("K", ("TMT6plex",)),
                ("C-term", ()),
            ),
            2,
        )

        self.assertNotIn("MH+^{13}C^{+}", frag_ions)
        self.assertNotIn("MH+2 ^{13}C^{+}", frag_ions)

        frag_ions = fragments.fragment_ions(
            (
                ("N-term", ("TMT6plex",)),
                ("S", ()),
                ("V", ()),
                ("Y", ("Phospho",)),
                ("T", ()),
                ("E", ()),
                ("I", ()),
                ("K", ("TMT6plex",)),
                ("C-term", ()),
            ),
            2,
            c13_num=1,
        )

        self.assertIn("MH+^{13}C^{+}", frag_ions)
        self.assertNotIn("MH+2 ^{13}C^{+}", frag_ions)

        frag_ions = fragments.fragment_ions(
            (
                ("N-term", ("TMT6plex",)),
                ("S", ()),
                ("V", ()),
                ("Y", ("Phospho",)),
                ("T", ()),
                ("E", ()),
                ("I", ()),
                ("K", ("TMT6plex",)),
                ("C-term", ()),
            ),
            2,
            c13_num=2,
        )

        self.assertIn("MH+^{13}C^{+}", frag_ions)
        self.assertIn("MH+2 ^{13}C^{+}", frag_ions)
